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
#include "galois/gstl.h"
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
  parSP1,
};

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep",    "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo",      "topoTile", "serSP1", "parSP1"};

static cll::opt<Algo>
    algo("algo", cll::desc("Choose an algorithm:"),
         cll::values(clEnumVal(deltaTile, "deltaTile"),
                     clEnumVal(deltaStep, "deltaStep"),
                     clEnumVal(serDeltaTile, "serDeltaTile"),
                     clEnumVal(serDelta, "serDelta"),
                     clEnumVal(dijkstraTile, "dijkstraTile"),
                     clEnumVal(dijkstra, "dijkstra"), clEnumVal(topo, "topo"),
                     clEnumVal(topoTile, "topoTile"),
                     clEnumVal(serSP1, "serSP1"),
                     clEnumVal(parSP1, "parSP1"), clEnumValEnd),
         cll::init(deltaTile));

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

template<typename GraphTypes>
struct SerNodeData {
  typename GraphTypes::SerGraphType::GraphNode node;
  int pred;
  uint32_t dist;
  bool fixed;
  SerNodeData(): node(-1), pred(0), dist(0), fixed(false) {}
};

template<typename GraphTypes>
struct ParNodeData {
  typename GraphTypes::ParGraphType::GraphNode node;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;
  ParNodeData(): node(-1), pred(0), dist(0), fixed(false) {}
};

struct GraphTypes {
  typedef SerNodeData<GraphTypes> SerNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::SerNodeDataType, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type SerGraphType;
  typedef ParNodeData<GraphTypes> ParNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::ParNodeDataType, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type ParGraphType;
};

template<typename GType>
struct GraphTypeTraits {
  typedef BFS_SSSP<GType, uint32_t, true, EDGE_TILE_SIZE> SSSP;
  typedef typename SSSP::Dist Dist;
  typedef typename SSSP::UpdateRequest UpdateRequest;
  typedef typename SSSP::UpdateRequestIndexer UpdateRequestIndexer;
  typedef typename SSSP::SrcEdgeTile SrcEdgeTile;
  typedef typename SSSP::SrcEdgeTileMaker SrcEdgeTileMaker;
  typedef typename SSSP::SrcEdgeTilePushWrap SrcEdgeTilePushWrap;
  typedef typename SSSP::ReqPushWrap ReqPushWrap;
  typedef typename SSSP::OutEdgeRangeFn OutEdgeRangeFn;
  typedef typename SSSP::TileRangeFn TileRangeFn;
};

template<typename GType>
struct GraphType {
  typedef GType Graph;
  typedef typename GType::GraphNode GNode;
  typedef GraphTypeTraits<GType> Traits;
};

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void deltaStepAlgo(Graph& graph,
                   GNode source,
                   const P& pushWrap,
                   const R& edgeRange) {

  using Dist                 = typename Traits::Dist;
  using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> BadWork;
  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> WLEmptyWork;

  namespace gwl = galois::worklists;

  using PSchunk = gwl::PerSocketChunkFIFO<CHUNK_SIZE>;
  using OBIM    = gwl::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

  graph.getData(source).dist = 0;

  galois::InsertBag<T> initBag;
  pushWrap(initBag, source, 0, "parallel");

  galois::for_each(galois::iterate(initBag),
                   [&](const T& item, auto& ctx) {
                     constexpr galois::MethodFlag flag =
                         galois::MethodFlag::UNPROTECTED;
                     const auto& sdata = graph.getData(item.src, flag);

                     if (sdata.dist < item.dist) {
                       if (TRACK_WORK)
                         WLEmptyWork += 1;
                       return;
                     }

                     for (auto ii : edgeRange(item)) {

                       GNode dst          = graph.getEdgeDst(ii);
                       auto& ddist        = graph.getData(dst, flag);
                       Dist ew            = graph.getEdgeData(ii, flag);
                       const Dist newDist = sdata.dist + ew;

                       while (true) {
                         Dist oldDist = ddist.dist;

                         if (oldDist <= newDist) {
                           break;
                         }

                         if (ddist.dist.compare_exchange_weak(
                                 oldDist, newDist, std::memory_order_relaxed)) {

                           if (TRACK_WORK) {
                             //! [per-thread contribution of self-defined stats]
                             if (oldDist != GType::Traits::SSSP::DIST_INFINITY) {
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

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void serDeltaAlgo(Graph& graph, const GNode& source,
                  const P& pushWrap, const R& edgeRange) {

  using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

  SerialBucketWL<T, UpdateRequestIndexer> wl(UpdateRequestIndexer{stepShift});
  ;
  graph.getData(source).dist = 0;

  pushWrap(wl, source, 0);

  size_t iter = 0ul;
  while (!wl.empty()) {

    auto& curr = wl.minBucket();

    while (!curr.empty()) {
      ++iter;
      auto item = curr.front();
      curr.pop_front();

      if (graph.getData(item.src).dist < item.dist) {
        // empty work
        continue;
      }

      for (auto e : edgeRange(item)) {

        GNode dst   = graph.getEdgeDst(e);
        auto& ddata = graph.getData(dst);

        const auto newDist = item.dist + graph.getEdgeData(e);

        if (newDist < ddata.dist) {
          ddata.dist = newDist;
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

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void dijkstraAlgo(Graph& graph, const GNode& source, const P& pushWrap,
                  const R& edgeRange) {

  using WL = galois::MinHeap<T>;

  graph.getData(source).dist = 0;

  WL wl;
  pushWrap(wl, source, 0);

  size_t iter = 0;
  size_t inner_iter = 0;
  size_t heap_pushes =0;
  double average_heap_size = 0;
  size_t duplicated_items = 0;

  while (!wl.empty()) {
    ++iter;
    average_heap_size += wl.size();

    T item = wl.pop();

    if (graph.getData(item.src).dist < item.dist) {
      duplicated_items++;
      continue;
    }

    for (auto e : edgeRange(item)) {

      inner_iter++;
      GNode dst   = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst);

      const auto newDist = item.dist + graph.getEdgeData(e);

      if (newDist < ddata.dist) {
        ddata.dist = newDist;
        pushWrap(wl, dst, newDist);
        heap_pushes++;
      }
    }
  }

  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Iterations", iter);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Inner iteration", inner_iter);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Average heap size", average_heap_size / iter);
}

template<typename Graph, typename R>
void calc_graph_predecessors(Graph& graph, R& edgeRange) {
  // Fill the pred array
  for (auto vertex : graph) {
    for (auto edge : edgeRange(vertex)) {
      graph.getData(graph.getEdgeDst(edge)).pred++;
    }
  }
}

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void serSP1Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange) {

  using Heap = galois::MinHeap<T>;

  Heap heap;
  pushWrap(heap, source, 0);

  // The set of nodes which have been fixed but not explored
  galois::gstl::Vector<GNode> r_set;

  size_t outer_iter = 0;
  size_t middle_iter = 0;
  size_t additional_nodes_explored = 0;
  size_t inner_iter = 0;
  size_t heap_pushes = 0;
  size_t average_heap_size = 0;
  size_t average_rset_size = 0;
  size_t duplicated_items = 0;
  size_t duplicates_avoided = 0;

  // While the heap is not empty
  while (!heap.empty()) {
    average_heap_size += heap.size();
    // Get the min element from the heap.
    T j = heap.pop();
    outer_iter++;

    auto& j_data = graph.getData(j.src);

    if (j_data.dist < j.dist) {
      duplicated_items++;
      continue;
    }

    // If the min element is not fixed
    if (!j_data.fixed) {
      // Set the element to fixed
      j_data.fixed = true;
      GNode z = j.src;

      // Inner loop, go through the all the elements in R
      while(true) {
        auto& z_data = graph.getData(z);
        // Get all the vertices that have edges from z
        for (auto e : edgeRange(z)) {
          inner_iter++;
          GNode k = graph.getEdgeDst(e);
          auto& k_data = graph.getData(k);
          // If k vertex is not fixed, process the edge between z and k.
          if (!k_data.fixed) {
            auto z_k_dist = graph.getEdgeData(e);
            k_data.pred--;
            if (k_data.pred <= 0) {
              k_data.fixed = true;
              r_set.push_back(k);
            }

            if (k_data.dist > z_data.dist + z_k_dist) {
              k_data.dist = z_data.dist + z_k_dist;
              if (!k_data.fixed) {
                heap_pushes++;
                pushWrap(heap, k_data.node, k_data.dist);
              }
            }
          }
        }

        average_rset_size += r_set.size();
        middle_iter++;

        if (r_set.empty()) break;
        z = r_set.back();
        r_set.pop_back();
        additional_nodes_explored++;
      }
    }
  }

  galois::runtime::reportStat_Single("SSSP-serSP1", "Outer loop iterations", outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Inner loop iterations", inner_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes avoided", duplicates_avoided);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Additional nodes explored", additional_nodes_explored);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Average heap size", (average_heap_size * 1.0) / outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Average rset size", (average_rset_size + 1.0) / middle_iter);
}

template<typename GNode, typename GEdge>
struct OutEdge {
  GNode src;
  GEdge edge;

  // OutEdge(const GNode& _src, const GEdge& _edge) : src(_src), edge(_edge) {}

  // // Move constructor.
  // OutEdge(OutEdge&& other) : src(std::move(other.src)), edge(std::move(other.edge)) {
  // }

  // // Copy constructor.
  // OutEdge(const OutEdge& other) : src(0), edge(0) {
  //   src = other.src;
  //   edge = other.edge;
  // }

  // OutEdge& operator=(const OutEdge& other) {
  //   src = other.src;
  //   edge = other.edge;
  //   return *this;
  // }

  //  OutEdge& operator=(OutEdge&& other)
  // {
  //  src = std::move(other.src);
  //  edge = std::move(other.edge);
  //  return *this;
  // }
};

template<typename GNode,
         typename GEdge,
         typename Container,
         typename R>
inline void push_back_edges_of_node(const GNode& node,
                                    Container& c,
                                    const R& edgeRange) {
  for (auto& edge : edgeRange(node)) {
    std::cout << "Edge iter: " << boost::core::demangle(typeid(edge).name()) << std::endl;
    std::cout << "Edge: " << boost::core::demangle(typeid(*edge).name()) << std::endl;
    std::cout << "Container: " << boost::core::demangle(typeid(c).name()) << std::endl;
    OutEdge<GNode, GEdge> out_edge{node, *edge};
    std::cout << "OutEdge: " << boost::core::demangle(typeid(out_edge).name()) << std::endl;
    std::cout << "Adding to r_set: " << node << " edge: " << *edge << std::endl;
    c.push_back(out_edge);
  }
}

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename GEdge = typename GType::Graph::edge_iterator::value_type,
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits>
void parSP1Algo(Graph& graph,
                const GNode& source,
                const P& pushWrap,
                const R& edgeRange) {

  using Heap = galois::MinHeap<T>;
  using OutEdge = OutEdge<GNode, GEdge>;
  using PSchunk = galois::worklists::PerSocketChunkLIFO<16>;

  Heap heap;
  pushWrap(heap, source, 0);

  // The set of nodes which have been fixed but not explored
  // This is a parallel data structure that is optimized for
  // parallel writes, but requires serial reads.
  galois::InsertBag<GNodeData*> q_set;

  // A step in an unbundled r set which looks at a single edge
  // from a node in the r_set.
  auto r_set_edge_loop = [&](OutEdge& z_edge, auto& ctx) {
    std::cout << "Edge: " << z_edge.edge << std::endl;
    std::cout << "Src Node: " << z_edge.src << std::endl;
    auto& z_data = graph.getData(z_edge.src);
    auto k = graph.getEdgeDst(z_edge.edge);
    std::cout << "Dst Node: " << k << std::endl;
    auto& k_data = graph.getData(k);
    // If k vertex is not fixed, process the edge between z and k.
    if (!k_data.fixed) {
      auto z_k_dist = graph.getEdgeData(z_edge.edge);
      k_data.pred--;
      if (k_data.pred <= 0) {
        k_data.fixed = true;
        push_back_edges_of_node<GNode, GEdge>(k, ctx, edgeRange);
      }

      if (k_data.dist > z_data.dist + z_k_dist) {
        k_data.dist = z_data.dist + z_k_dist;
        if (!k_data.fixed) {
          q_set.push_back(&k_data);
        }
      }
    }
  };

  // While the heap is not empty
  while (!heap.empty()) {
    // Get the min element from the heap.
    T j = heap.pop();
    auto& j_data = graph.getData(j.src);

    std::cout << "Popping from heap: " << j_data.node << " : " << j_data.dist << std::endl;

    if (j_data.dist < j.dist) {
      continue;
    }

    // If the min element is not fixed
    if (!j_data.fixed) {
      // Set the element to fixed
      j_data.fixed = true;

      galois::gstl::Vector<OutEdge> initial;
      push_back_edges_of_node<GNode, GEdge>(j.src, initial, edgeRange);
      galois::for_each(galois::iterate(initial),
                       r_set_edge_loop,
                       galois::wl<PSchunk>(),
                       galois::loopname("sssp_parsp1_inner"));
    }

    for (auto& z : q_set) {
      auto& z_data = *z;
      std::cout << "Pushing to heap: " << z_data.node << " : " << z_data.dist << std::endl;
      pushWrap(heap, z_data.node, z_data.dist);
    }

    q_set.clear();
  }
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void topoAlgo(Graph& graph, const GNode& source) {
  using SSSP                 = typename Traits::SSSP;
  using Dist                 = typename Traits::Dist;

  galois::LargeArray<Dist> oldDist;
  oldDist.allocateInterleaved(graph.size());

  constexpr Dist INFTY = SSSP::DIST_INFINITY;
  galois::do_all(galois::iterate(0ul, graph.size()),
                 [&](size_t i) { oldDist.constructAt(i, INFTY); },
                 galois::no_stats(), galois::loopname("initDistArray"));

  graph.getData(source).dist = 0;

  galois::GReduceLogicalOR changed;
  size_t rounds = 0;

  do {

    ++rounds;
    changed.reset();

    galois::do_all(galois::iterate(graph),
                   [&](const GNode& n) {
                     const auto& sdata = graph.getData(n);

                     if (oldDist[n] > sdata.dist) {

                       oldDist[n] = sdata.dist;
                       changed.update(true);

                       for (auto e : graph.edges(n)) {
                         const auto newDist = sdata.dist + graph.getEdgeData(e);
                         auto dst           = graph.getEdgeDst(e);
                         auto& ddata        = graph.getData(dst);
                         galois::atomicMin(ddata.dist, newDist);
                       }
                     }
                   },
                   galois::steal(), galois::loopname("Update"));

  } while (changed.reduce());

  galois::runtime::reportStat_Single("SSSP-topo", "rounds", rounds);
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void topoTileAlgo(Graph& graph, const GNode& source) {

  using SSSP                 = typename Traits::SSSP;
  using SrcEdgeTile          = typename Traits::SrcEdgeTile;
  using SrcEdgeTileMaker     = typename Traits::SrcEdgeTileMaker;

  galois::InsertBag<SrcEdgeTile> tiles;

  graph.getData(source).dist = 0;

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

                     if (t.dist > sdata.dist) {

                       t.dist = sdata.dist;
                       changed.update(true);

                       for (auto e = t.beg; e != t.end; ++e) {
                         const auto newDist = sdata.dist + graph.getEdgeData(e);
                         auto dst           = graph.getEdgeDst(e);
                         auto& ddata        = graph.getData(dst);
                         galois::atomicMin(ddata.dist, newDist);
                       }
                     }
                   },
                   galois::steal(), galois::loopname("Update"));

  } while (changed.reduce());

  galois::runtime::reportStat_Single("SSSP-topo", "rounds", rounds);
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void init_graph(Graph& graph, GNode& source, GNode& report) {

  using SSSP = typename Traits::SSSP;
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
                 [&graph](GNode n) {
                   auto& data = graph.getData(n);
                   data.dist = SSSP::DIST_INFINITY;
                   data.node = n;
                 });

  graph.getData(source).dist = 0;
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void verify_and_report(Graph& graph, GNode& source, GNode& report) {
  galois::reportPageAlloc("MeminfoPost");

  std::cout << "Node " << reportNode << " has distance "
            << graph.getData(report).dist << "\n";

  if (!skipVerify) {
    if (Traits::SSSP::verify(graph, source)) {
      std::cout << "Verification successful.\n";
    } else {
      GALOIS_DIE("Verification failed");
    }
  }
}

int main(int argc, char** argv) {
  using ParGType = GraphType<GraphTypes::ParGraphType>;
  using SerGType = GraphType<GraphTypes::SerGraphType>;

  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;
  galois::StatTimer Tmain;

  switch (algo) {
    case deltaTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      deltaStepAlgo<ParGType, ParGType::Traits::SrcEdgeTile>(
          graph, source,
          ParGType::Traits::SrcEdgeTilePushWrap{graph},
          ParGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report);
      break;
    }
    case deltaStep: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      deltaStepAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          ParGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report);
      break;
    }
    case serDeltaTile: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);

      Tmain.start();
      serDeltaAlgo<SerGType, SerGType::Traits::SrcEdgeTile>(
          graph, source,
          SerGType::Traits::SrcEdgeTilePushWrap{graph},
          SerGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report);
      break;
    }
    case serDelta: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);

      Tmain.start();
      serDeltaAlgo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          SerGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report);
      break;
    }
    case dijkstraTile: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);

      Tmain.start();
      dijkstraAlgo<SerGType, SerGType::Traits::SrcEdgeTile>(
          graph, source,
          SerGType::Traits::SrcEdgeTilePushWrap{graph},
          SerGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report);
      break;
    }
    case dijkstra: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);

      Tmain.start();
      dijkstraAlgo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          SerGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report);
      break;
    }
    case topo: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      topoAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report);
      break;
    }
    case topoTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      topoTileAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report);
      break;
    }
    case serSP1: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      Tmain.start();
      serSP1Algo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report);
      break;
    }
    case parSP1: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      Tmain.start();
      parSP1Algo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report);
      break;
    }
    default:
      std::abort();
  }
  return 0;
}
