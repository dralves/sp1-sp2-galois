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
  serSP2,
};

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep",    "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo",      "topoTile", "serSP1", "serSP2"};

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
		     clEnumVal(serSP2, "serSP2"), clEnumValEnd),
         cll::init(deltaTile));

struct NodeData;

// typedef galois::graphs::LC_InlineEdge_Graph<std::atomic<unsigned int>,
// uint32_t>::with_no_lockable<true>::type::with_numa_alloc<true>::type Graph;
//! [withnumaalloc]
using InnerGraph = galois::graphs::LC_CSR_Graph<NodeData, uint32_t>::
    with_no_lockable<true>::type ::with_numa_alloc<true>::type;
using Graph = galois::graphs::LC_InOut_Graph<InnerGraph>;
//! [withnumaalloc]

using GNode = Graph::GraphNode;

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

using SSSP                 = BFS_SSSP<Graph, uint32_t, true, EDGE_TILE_SIZE>;

struct NodeData {
  std::atomic<uint32_t> dist;
  int pred;
  unsigned int minWeight; 
  bool fixed;
  bool in_heap;
  uint32_t heap_dist;
  GNode node;

  NodeData(): dist(0), pred(0), minWeight(SSSP::DIST_INFINITY), fixed(false), in_heap(false), heap_dist(0), node(-1) {}
};

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

template <typename T, typename P, typename R>
void dijkstraAlgo(Graph& graph, const GNode& source, const P& pushWrap,
                  const R& edgeRange) {

  using WL = galois::MinHeap<T>;

  graph.getData(source).dist = 0;

  WL wl;
  pushWrap(wl, source, 0);
  
  #ifdef STATS
  size_t iter = 0;
  size_t inner_iter = 0;
  size_t heap_pushes =0;
  double average_heap_size = 0;
  size_t duplicated_items = 0;
  #endif

  while (!wl.empty()) {
    #ifdef STATS
    ++iter;
    average_heap_size += wl.size();
    #endif

    T item = wl.pop();

    if (graph.getData(item.src).dist < item.dist) {
      #ifdef STATS
      duplicated_items++;
      #endif
      continue;
    }

    for (auto e : edgeRange(item)) {
      #ifdef STATS
      inner_iter++;
      #endif
      GNode dst = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst);

      const auto newDist = item.dist + graph.getEdgeData(e);

      if (newDist < ddata.dist) {
        ddata.dist = newDist;
        pushWrap(wl, dst, newDist);
        #ifdef STATS
        heap_pushes++;
	#endif
      }
    }
  }
  

  #ifdef STATS
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Iterations", iter);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Inner iteration", inner_iter);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Average heap size", average_heap_size / iter);
  #endif
}

template<typename R>
void initGraph(Graph& graph, R& edgeRange) {
  galois::gstl::Deque<GNode> pred_set;
  // Fill the pred array
  for (auto vertex : graph) {
    for (auto edge : edgeRange(vertex)) {
      graph.getData(graph.getEdgeDst(edge)).pred++;
      if((graph.getData(graph.getEdgeDst(edge)).minWeight) > graph.getEdgeData(edge))
	graph.getData(graph.getEdgeDst(edge)).minWeight = graph.getEdgeData(edge);
    }
  }
  //add vertices with degree 0 to a queue
  for (auto vertex : graph){
    if(graph.getData(vertex).pred == 0)
      pred_set.push_back(vertex);
  }
  //decrement pred count of neighbours having an incoming edge from a vertex with degree 0
  while(!pred_set.empty()){
     GNode Vertex = pred_set.front();
     pred_set.pop_front();
     for (auto edge : edgeRange(Vertex)){
       graph.getData(graph.getEdgeDst(edge)).pred--;
       if(graph.getData(graph.getEdgeDst(edge)).pred == 0){
         pred_set.push_back(graph.getEdgeDst(edge));
       }
     }
  }
}

template <typename T, typename P, typename R>
void serSP1Algo(Graph& graph, const GNode& source, const P& pushWrap,
                const R& edgeRange) {

  using Heap = galois::MinHeap<T>;

  Heap heap;
  pushWrap(heap, source, 0);

  // The set of nodes which have been fixed but not explored
  galois::gstl::Vector<GNode> r_set;
  galois::gstl::Vector<NodeData*> q_set;

  #ifdef STATS
  size_t outer_iter = 0;
  size_t middle_iter = 0;
  size_t additional_nodes_explored = 0;
  size_t inner_iter = 0;
  size_t heap_pushes = 0;
  size_t average_heap_size = 0;
  size_t average_rset_size = 0;
  size_t duplicated_items = 0;
  size_t duplicates_avoided = 0;
  #endif

  // While the heap is not empty
  while (!heap.empty()) {
    #ifdef STATS
    outer_iter++;
    average_heap_size += heap.size();
    #endif
    // Get the min element from the heap.
    T j = heap.pop();

    auto& j_data = graph.getData(j.src);

    if (j_data.dist < j.dist) {
      #ifdef STATS
      duplicated_items++;
      #endif
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
      	  #ifdef STATS
          inner_iter++;
	  #endif
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
		#ifdef STATS
                heap_pushes++;
		#endif
                pushWrap(heap, k_data.node, k_data.dist);
              }
            }
          }
        }
	#ifdef STATS
        average_rset_size += r_set.size();
        middle_iter++;
	#endif

        if (r_set.empty()) break;
        z = r_set.back();
        r_set.pop_back();
	#ifdef STATS
        additional_nodes_explored++;
	#endif
      }
    }
  }

  #ifdef STATS
  galois::runtime::reportStat_Single("SSSP-serSP1", "Outer loop iterations", outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Inner loop iterations", inner_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Duplicated heap pushes avoided", duplicates_avoided);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Additional nodes explored", additional_nodes_explored);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Average heap size", (average_heap_size * 1.0) / outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP1", "Average rset size", (average_rset_size + 1.0) / middle_iter);
  #endif
}


template <typename T, typename P, typename R>
void serSP2Algo(Graph& graph, const GNode& source, const P& pushWrap,
                const R& edgeRange) {

  using Heap = galois::MinHeap<T>;

  Heap heap;
  pushWrap(heap, source, 0);
  auto d = 0; 

  // The set of nodes which have been fixed but not explored implemented as a queue
  galois::gstl::Deque<GNode> r_set;
  galois::gstl::Vector<NodeData*> q_set;
  
  #ifdef STATS
  size_t outer_iter = 0;
  size_t middle_iter = 0;
  size_t additional_nodes_explored = 0;
  size_t inner_iter = 0;
  size_t heap_pushes = 0;
  size_t average_heap_size = 0;
  size_t average_rset_size = 0;
  size_t duplicated_items = 0;
  size_t duplicates_avoided = 0;
  #endif

  // While the heap is not empty
  while (!heap.empty()) {
    #ifdef STATS
    average_heap_size += heap.size();
    outer_iter++;
    #endif
    // Get the min element from the heap.
    T j = heap.pop();

    auto& j_data = graph.getData(j.src);

    if (j_data.dist < j.dist) {
      #ifdef STATS
      duplicated_items++;
      #endif
      continue;
    }
    
    // If the min element is not fixed
    if (!j_data.fixed) {
      // Set the element to fixed
      j_data.fixed = true;
      GNode z = j.src;
      d = j_data.dist;

      while(!heap.empty()){
	T j_ch = heap.read();
	auto& j_ch_data = graph.getData(j_ch.src);

	if( j_ch_data.dist == j_data.dist && j_ch_data.fixed != true ){
          j_ch_data.fixed = true;
	  GNode z_ch = j_ch.src;
	  r_set.push_back(z_ch);
	  heap.pop();
	  continue;
	}
	break;
      }

      // Inner loop, go through the all the elements in R
      while(true) {
        auto& z_data = graph.getData(z);
        // Get all the vertices that have edges from z
        for (auto e : edgeRange(z)) {
          #ifdef STATS
          inner_iter++;
	  #endif
          GNode k = graph.getEdgeDst(e);
          auto& k_data = graph.getData(k);
          // If k vertex is not fixed, process the edge between z and k.
          if (!k_data.fixed) {

            auto& z_k_dist = graph.getEdgeData(e);
	    k_data.pred--;
	    auto k_dist = z_data.dist + z_k_dist;

            if ((k_data.pred <=0) || k_dist <= (d + k_data.minWeight)){
                k_data.fixed = true;
                r_set.push_back(k);
            }

            if (k_data.dist > z_data.dist + z_k_dist){
                k_data.dist = z_data.dist + z_k_dist;
            	if (!k_data.fixed) {
	          q_set.push_back(&k_data);
            	}
	    }
          }
        }

    	#ifdef STATS
        average_rset_size += r_set.size();
        middle_iter++;
	#endif

        if (r_set.empty()) break;
        z = r_set.front();
        r_set.pop_front();
    	#ifdef STATS
        additional_nodes_explored++;
	#endif
      }

      for (auto& z : q_set) {
        auto& z_data = *z;
    	#ifdef STATS
        heap_pushes++;
	#endif
        pushWrap(heap, z_data.node, z_data.dist);
      }
      q_set.clear();
    }
  }

  #ifdef STATS
  galois::runtime::reportStat_Single("SSSP-serSP2", "Outer loop iterations", outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Inner loop iterations", inner_iter);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Additional nodes explored", additional_nodes_explored);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Average heap size", (average_heap_size * 1.0) / outer_iter);
  galois::runtime::reportStat_Single("SSSP-serSP2", "Average rset size", (average_rset_size + 1.0) / middle_iter);
  #endif

}

void topoAlgo(Graph& graph, const GNode& source) {

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

void topoTileAlgo(Graph& graph, const GNode& source) {

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
                 [&graph](GNode n) {
                   auto& data = graph.getData(n);
                   data.dist = SSSP::DIST_INFINITY;
                   data.node = n;
                 });

  graph.getData(source).dist = 0;

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;
  galois::StatTimer Tmain;

  switch (algo) {
  case deltaTile:
    Tmain.start();
    deltaStepAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                               TileRangeFn());
    break;
  case deltaStep:
    Tmain.start();
    deltaStepAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                 OutEdgeRangeFn{graph});
    break;
  case serDeltaTile:
    Tmain.start();
    serDeltaAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                              TileRangeFn());
    break;
  case serDelta:
    Tmain.start();
    serDeltaAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                OutEdgeRangeFn{graph});
    break;
  case dijkstraTile:
    Tmain.start();
    dijkstraAlgo<SrcEdgeTile>(graph, source, SrcEdgeTilePushWrap{graph},
                              TileRangeFn());
    break;
  case dijkstra:
    Tmain.start();
    dijkstraAlgo<UpdateRequest>(graph, source, ReqPushWrap(),
                                OutEdgeRangeFn{graph});
    break;
  case topo:
    Tmain.start();
    topoAlgo(graph, source);
    break;
  case topoTile:
    Tmain.start();
    topoTileAlgo(graph, source);
    break;
    case serSP1: {
      auto edgeRange = OutEdgeRangeFn{graph};
      initGraph(graph, edgeRange);

      Tmain.start();
      serSP1Algo<UpdateRequest>(graph, source, ReqPushWrap(), edgeRange);
      break;
    }

    case serSP2: {
      auto edgeRange = OutEdgeRangeFn{graph};
      initGraph(graph, edgeRange);

      Tmain.start();
      serSP2Algo<UpdateRequest>(graph, source, ReqPushWrap(), edgeRange);
      break;
    }


  default:
    std::abort();
  }

  Tmain.stop();

  galois::reportPageAlloc("MeminfoPost");

  std::cout << "Node " << reportNode << " has distance "
            << graph.getData(report).dist << "\n";

  if (!skipVerify) {
    if (SSSP::verify(graph, source)) {
      std::cout << "Verification successful.\n";
    } else {
      GALOIS_DIE("Verification failed");
    }
  }

  return 0;
}
