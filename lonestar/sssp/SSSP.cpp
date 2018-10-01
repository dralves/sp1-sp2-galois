/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/FlatMap.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "llvm/Support/CommandLine.h"

#include "Lonestar/BoilerPlate.h"
#include "Lonestar/BFS_SSSP.h"

#include <iostream>
#include <unordered_set>

namespace cll = llvm::cl;

static const char* name = "Single Source Shortest Path";
static const char* desc =
    "Computes the shortest path from a source node to all nodes in a directed "
    "graph using a modified chaotic iteration algorithm";
static const char* url = "single_source_shortest_path";

static cll::opt<std::string>
    filename(cll::Positional, cll::desc("<input graph>"), cll::Required);

static cll::opt<unsigned int>
    startNode("startNode",
              cll::desc("Node to start search from (default value 0)"),
              cll::init(0));
static cll::opt<unsigned int>
    reportNode("reportNode",
               cll::desc("Node to report distance to(default value 1)"),
               cll::init(1));
static cll::opt<unsigned int>
    stepShift("delta",
              cll::desc("Shift value for the deltastep (default value 13)"),
              cll::init(13));

enum Algo {
  deltaTile = 0,
  deltaStep,
  serDeltaTile,
  serDelta,
  dijkstraTile,
  dijkstra,
  topo,
  topoTile,
  serSP1,
};

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep",    "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo",      "topoTile", "serSP1"};

static cll::opt<Algo>
    algo("algo", cll::desc("Choose an algorithm:"),
         cll::values(clEnumVal(deltaTile, "deltaTile"),
                     clEnumVal(deltaStep, "deltaStep"),
                     clEnumVal(serDeltaTile, "serDeltaTile"),
                     clEnumVal(serDelta, "serDelta"),
                     clEnumVal(dijkstraTile, "dijkstraTile"),
                     clEnumVal(dijkstra, "dijkstra"), clEnumVal(topo, "topo"),
                     clEnumVal(topoTile, "topoTile"),
                     clEnumVal(serSP1, "serSP1"), clEnumValEnd),
         cll::init(deltaTile));

// typedef galois::graphs::LC_InlineEdge_Graph<std::atomic<unsigned int>,
// uint32_t>::with_no_lockable<true>::type::with_numa_alloc<true>::type Graph;
//! [withnumaalloc]
using Graph = galois::graphs::LC_CSR_Graph<std::atomic<uint32_t>, uint32_t>::
    with_no_lockable<true>::type ::with_numa_alloc<true>::type;
//! [withnumaalloc]
typedef Graph::GraphNode GNode;

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

using SSSP                 = BFS_SSSP<Graph, uint32_t, true, EDGE_TILE_SIZE>;
using Dist                 = SSSP::Dist;
using UpdateRequest        = SSSP::UpdateRequest;
using UpdateRequestIndexer = SSSP::UpdateRequestIndexer;
using SrcEdgeTile          = SSSP::SrcEdgeTile;
using SrcEdgeTileMaker     = SSSP::SrcEdgeTileMaker;
using SrcEdgeTilePushWrap  = SSSP::SrcEdgeTilePushWrap;
using ReqPushWrap          = SSSP::ReqPushWrap;
using OutEdgeRangeFn       = SSSP::OutEdgeRangeFn;
using TileRangeFn          = SSSP::TileRangeFn;

template <typename T, typename P, typename R>
void deltaStepAlgo(Graph& graph, GNode source, const P& pushWrap,
                   const R& edgeRange) {

  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> BadWork;
  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> WLEmptyWork;

  namespace gwl = galois::worklists;

  using PSchunk = gwl::PerSocketChunkFIFO<CHUNK_SIZE>;
  using OBIM    = gwl::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

  graph.getData(source) = 0;

  galois::InsertBag<T> initBag;
  pushWrap(initBag, source, 0, "parallel");

  galois::for_each(galois::iterate(initBag),
                   [&](const T& item, auto& ctx) {
                     constexpr galois::MethodFlag flag =
                         galois::MethodFlag::UNPROTECTED;
                     const auto& sdata = graph.getData(item.src, flag);

                     if (sdata < item.dist) {
                       if (TRACK_WORK)
                         WLEmptyWork += 1;
                       return;
                     }

                     for (auto ii : edgeRange(item)) {

                       GNode dst          = graph.getEdgeDst(ii);
                       auto& ddist        = graph.getData(dst, flag);
                       Dist ew            = graph.getEdgeData(ii, flag);
                       const Dist newDist = sdata + ew;

                       while (true) {
                         Dist oldDist = ddist;

                         if (oldDist <= newDist) {
                           break;
                         }

                         if (ddist.compare_exchange_weak(
                                 oldDist, newDist, std::memory_order_relaxed)) {

                           if (TRACK_WORK) {
                             //! [per-thread contribution of self-defined stats]
                             if (oldDist != SSSP::DIST_INFINITY) {
                               BadWork += 1;
                             }
                             //! [per-thread contribution of self-defined stats]
                           }

                           pushWrap(ctx, dst, newDist);
                           break;
                         }
                       }
                     }
                   },
                   galois::wl<OBIM>(UpdateRequestIndexer{stepShift}),
                   galois::no_conflicts(), galois::loopname("SSSP"));

  if (TRACK_WORK) {
    //! [report self-defined stats]
    galois::runtime::reportStat_Single("SSSP", "BadWork", BadWork.reduce());
    //! [report self-defined stats]
    galois::runtime::reportStat_Single("SSSP", "WLEmptyWork",
                                       WLEmptyWork.reduce());
  }
}

template <typename T, typename P, typename R>
void serDeltaAlgo(Graph& graph, const GNode& source, const P& pushWrap,
                  const R& edgeRange) {

  SerialBucketWL<T, UpdateRequestIndexer> wl(UpdateRequestIndexer{stepShift});
  ;
  graph.getData(source) = 0;

  pushWrap(wl, source, 0);

  size_t iter = 0ul;
  while (!wl.empty()) {

    auto& curr = wl.minBucket();

    while (!curr.empty()) {
      ++iter;
      auto item = curr.front();
      curr.pop_front();

      if (graph.getData(item.src) < item.dist) {
        // empty work
        continue;
      }

      for (auto e : edgeRange(item)) {

        GNode dst   = graph.getEdgeDst(e);
        auto& ddata = graph.getData(dst);

        const auto newDist = item.dist + graph.getEdgeData(e);

        if (newDist < ddata) {
          ddata = newDist;
          pushWrap(wl, dst, newDist);
        }
      }
    }

    wl.goToNextBucket();
  }

  if (!wl.allEmpty()) {
    std::abort();
  }
  galois::runtime::reportStat_Single("SSSP-Serial-Delta", "Iterations", iter);
}

template <typename T, typename P, typename R>
void dijkstraAlgo(Graph& graph, const GNode& source, const P& pushWrap,
                  const R& edgeRange) {

  using WL = galois::MinHeap<T>;

  graph.getData(source) = 0;

  WL wl;
  pushWrap(wl, source, 0);

  size_t iter = 0;

  while (!wl.empty()) {
    ++iter;

    T item = wl.pop();

    if (graph.getData(item.src) < item.dist) {
      // empty work
      continue;
    }

    for (auto e : edgeRange(item)) {

      GNode dst   = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst);

      const auto newDist = item.dist + graph.getEdgeData(e);

      if (newDist < ddata) {
        ddata = newDist;
        pushWrap(wl, dst, newDist);
      }
    }
  }

  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Iterations", iter);
}

template<typename T,
         typename Hash = std::hash<T>,
         typename KeyEqual = std::equal_to<T>>
using UnorderedSet = std::unordered_set<T, Hash, KeyEqual, galois::runtime::Pow_2_BlockAllocator<T>>;

template<typename Set, typename Map, typename EdgeDistData, typename GraphDistData>
bool processEdgeSP1(Set& fixed, Set& r_set, Set& q_set, Map& pred, const GNode& k,
                    EdgeDistData& z_k_dist, GraphDistData& k_dist, GraphDistData& z_dist) {
  bool changed = false;
  pred[k] -= 1;
  if (k_dist > z_dist + z_k_dist) {
    k_dist = z_dist + z_k_dist;
    changed = true;
  }
  if (pred[k] == 0) {
    fixed.insert(k);
    r_set.insert(k);
  } else if (changed) {
    q_set.insert(k);
  }
  return changed;
}

template <typename T, typename P, typename R>
void serSP1Algo(Graph& graph, const GNode& source, const P& pushWrap,
                const R& edgeRange) {

  using Heap = galois::MinHeap<T>;
  using Map = galois::flat_map<GNode, int>;

  Heap heap;
  // Se the dist in the graph to 0 and push to the heap.
  graph.getData(source) = 0;
  pushWrap(heap, source, 0);

  Map pred;
  // Fill the pred array. This should be easily parallelizable.
  for (auto vertex : graph) {
    std::cout << vertex << std::endl;
    size_t incoming_edges = 0;
    for (auto edge : edgeRange(vertex)) {
      incoming_edges++;
    }
    pred.insert(std::pair<GNode, int>(vertex, incoming_edges));
  }

  // The set of nodes which have been fixed
  UnorderedSet<GNode> fixed;
  // The set of nodes whose distance has changed (used in the inner loop)
  UnorderedSet<GNode> q_set;
  // The set of nodes which have been fixed but not explored
  UnorderedSet<GNode> r_set;

  size_t outer_iter = 0;
  size_t inner_iter = 0;

  // While the heap is not empty
  while (!heap.empty()) {
    // Get the min element from the heap.
    T j = heap.pop();
    outer_iter++;

    // If the min element is not fixed
    if (fixed.find(j.src) == fixed.end()) {
      // Insert into R
      r_set.insert(j.src);
      // Set the element to fixed
      fixed.insert(j.src);

      // Inner loop
      // Go through the all the elements in R
      while (!r_set.empty()) {
        for (auto it = r_set.begin(); it != r_set.end();) {
          GNode z = *it;
          // Remove one element
          it = r_set.erase(it);
          // Get all the vertices that have edges to z
          for (auto e : edgeRange(z)) {
            inner_iter++;
            GNode k = graph.getEdgeDst(e);
            // If k vertex is not fixed, process the edge
            // between z and k.
            if (fixed.find(k) == fixed.end()) {
              auto& k_dist = graph.getData(k);
              auto& z_dist = graph.getData(z);
              auto z_k_dist = graph.getEdgeData(e);
              // Process the edge. If the distance changed we need to remove
              // the item from the heap if it was there. We'll reinsert below
              // when we go through the q set.
              if (processEdgeSP1(fixed, r_set, q_set, pred, k, z_k_dist, k_dist, z_dist)) {
                T k_item(k, 0);
                heap.remove(k_item);
              }
            }
          }
        }
      }
      for (auto z : q_set) {
        if (fixed.find(z) == fixed.end()) {
          pushWrap(heap, z, graph.getData(z));
        }
      }
      q_set.clear();
    }
  }

  galois::runtime::reportStat_Single("SSSP-serSP1", "Outer loop iterations", outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Inner loop iterations", inner_iter);
}

void topoAlgo(Graph& graph, const GNode& source) {

  galois::LargeArray<Dist> oldDist;
  oldDist.allocateInterleaved(graph.size());

  constexpr Dist INFTY = SSSP::DIST_INFINITY;
  galois::do_all(galois::iterate(0ul, graph.size()),
                 [&](size_t i) { oldDist.constructAt(i, INFTY); },
                 galois::no_stats(), galois::loopname("initDistArray"));

  graph.getData(source) = 0;

  galois::GReduceLogicalOR changed;
  size_t rounds = 0;

  do {

    ++rounds;
    changed.reset();

    galois::do_all(galois::iterate(graph),
                   [&](const GNode& n) {
                     const auto& sdata = graph.getData(n);

                     if (oldDist[n] > sdata) {

                       oldDist[n] = sdata;
                       changed.update(true);

                       for (auto e : graph.edges(n)) {
                         const auto newDist = sdata + graph.getEdgeData(e);
                         auto dst           = graph.getEdgeDst(e);
                         auto& ddata        = graph.getData(dst);
                         galois::atomicMin(ddata, newDist);
                       }
                     }
                   },
                   galois::steal(), galois::loopname("Update"));

  } while (changed.reduce());

  galois::runtime::reportStat_Single("SSSP-topo", "rounds", rounds);
}

void topoTileAlgo(Graph& graph, const GNode& source) {

  galois::InsertBag<SrcEdgeTile> tiles;

  graph.getData(source) = 0;

  galois::do_all(galois::iterate(graph),
                 [&](const GNode& n) {
                   SSSP::pushEdgeTiles(
                       tiles, graph, n,
                       SrcEdgeTileMaker{n, SSSP::DIST_INFINITY});
                 },
                 galois::steal(), galois::loopname("MakeTiles"));

  galois::GReduceLogicalOR changed;
  size_t rounds = 0;

  do {
    ++rounds;
    changed.reset();

    galois::do_all(galois::iterate(tiles),
                   [&](SrcEdgeTile& t) {
                     const auto& sdata = graph.getData(t.src);

                     if (t.dist > sdata) {

                       t.dist = sdata;
                       changed.update(true);

                       for (auto e = t.beg; e != t.end; ++e) {
                         const auto newDist = sdata + graph.getEdgeData(e);
                         auto dst           = graph.getEdgeDst(e);
                         auto& ddata        = graph.getData(dst);
                         galois::atomicMin(ddata, newDist);
                       }
                     }
                   },
                   galois::steal(), galois::loopname("Update"));

  } while (changed.reduce());

  galois::runtime::reportStat_Single("SSSP-topo", "rounds", rounds);
}

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  Graph graph;
  GNode source, report;

  std::cout << "Reading from file: " << filename << std::endl;
  galois::graphs::readGraph(graph, filename);
  std::cout << "Read " << graph.size() << " nodes, " << graph.sizeEdges()
            << " edges" << std::endl;

  if (startNode >= graph.size() || reportNode >= graph.size()) {
    std::cerr << "failed to set report: " << reportNode
              << " or failed to set source: " << startNode << "\n";
    assert(0);
    abort();
  }

  auto it = graph.begin();
  std::advance(it, startNode);
  source = *it;
  it     = graph.begin();
  std::advance(it, reportNode);
  report = *it;

  size_t approxNodeData = graph.size() * 64;
  // size_t approxEdgeData = graph.sizeEdges() * sizeof(typename
  // Graph::edge_data_type) * 2;
  galois::preAlloc(numThreads +
                   approxNodeData / galois::runtime::pagePoolSize());
  galois::reportPageAlloc("MeminfoPre");

  if (algo == deltaStep || algo == deltaTile || algo == serDelta ||
      algo == serDeltaTile) {
    std::cout << "INFO: Using delta-step of " << (1 << stepShift) << "\n";
    std::cout
        << "WARNING: Performance varies considerably due to delta parameter.\n";
    std::cout
        << "WARNING: Do not expect the default to be good for your graph.\n";
  }

  galois::do_all(galois::iterate(graph),
                 [&graph](GNode n) { graph.getData(n) = SSSP::DIST_INFINITY; });

  graph.getData(source) = 0;

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;

  galois::StatTimer Tmain;
  Tmain.start();

  switch (algo) {
  case deltaTile:
    deltaStepAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                               TileRangeFn());
    break;
  case deltaStep:
    deltaStepAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                 OutEdgeRangeFn{graph});
    break;
  case serDeltaTile:
    serDeltaAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                              TileRangeFn());
    break;
  case serDelta:
    serDeltaAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                OutEdgeRangeFn{graph});
    break;
  case dijkstraTile:
    dijkstraAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                              TileRangeFn());
    break;
  case dijkstra:
    dijkstraAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                OutEdgeRangeFn{graph});
    break;
  case topo:
    topoAlgo(graph, source);
    break;
  case topoTile:
    topoTileAlgo(graph, source);
    break;
  case serSP1:
    serSP1Algo<UpdateRequest>(graph, source, ReqPushWrap(),
                              OutEdgeRangeFn{graph});
    break;
  default:
    std::abort();
  }

  Tmain.stop();

  galois::reportPageAlloc("MeminfoPost");

  std::cout << "Node " << reportNode << " has distance "
            << graph.getData(report) << "\n";

  if (!skipVerify) {
    if (SSSP::verify(graph, source)) {
      std::cout << "Verification successful.\n";
    } else {
      GALOIS_DIE("Verification failed");
    }
  }

  return 0;
}
