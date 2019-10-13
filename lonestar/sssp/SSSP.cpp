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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

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
#include <chrono>
#include <iomanip>
#include <unordered_set>
#include <fstream>
#include "stack_vector.h"


#ifndef STACK_BITVECTOR_SIZE
#define STACK_BITVECTOR_SIZE 23947347
#endif

namespace cll = llvm::cl;
using namespace std::chrono; 
namespace acc = boost::accumulators;

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
  parSP1,
  parSP2
};

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep",    "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo", "topoTile", "serSP1", "serSP2", "parSP1", "parSP2"};

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
                     clEnumVal(serSP2, "serSP2"),
                     clEnumVal(parSP1, "parSP1"),
                     clEnumVal(parSP2, "parSP2"),
                     clEnumValEnd),
         cll::init(deltaTile));

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

high_resolution_clock::time_point start;

template<typename GraphTypes>
struct SerNodeData {
  typename GraphTypes::SerGraphType::GraphNode node;
  int pred;
  uint32_t dist;
  bool fixed;
  high_resolution_clock::time_point fixed_ts;
  uint32_t min_in_weight;
  SerNodeData(): node(-1), pred(0), dist(0), fixed(false) {}
};

template<typename GraphTypes>
struct ParNodeData {
  typename GraphTypes::ParGraphType::GraphNode node;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;
  high_resolution_clock::time_point fixed_ts;
  uint32_t min_in_weight;
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

std::vector<uint32_t> total_pred;

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
  using OBIM = galois::worklists::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

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

                           ddist.fixed = true;
                           ddist.fixed_ts = high_resolution_clock::now();
                           
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

  while (!wl.empty()) {

    T item = wl.pop();

    if (graph.getData(item.src).dist < item.dist) {
      continue;
    }

    auto& nd = graph.getData(item.src);
    nd.fixed = true;
    nd.fixed_ts = high_resolution_clock::now();

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
}

template<typename Graph, typename R>
void calc_graph_predecessors(Graph& graph, R& edgeRange) {
  // Fill the pred array
  for (auto vertex : graph) {
    for (auto edge : edgeRange(vertex)) {
      auto& k_data = graph.getData(graph.getEdgeDst(edge));
      k_data.pred++;
      total_pred[k_data.node]++;
      k_data.min_in_weight = std::min(k_data.min_in_weight, graph.getEdgeData(edge));
    }
  }
}

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#define LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits,
	  typename Cmp = std::less<typename Traits::Dist>>
void serSP1Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange) {

  using Heap = galois::MinHeap<T>;
  using Dist = typename Traits::Dist;

  bool changed;

  Heap heap;
  pushWrap(heap, source, 0);

  galois::gstl::Vector<GNodeData*> r_set;
  r_set.reserve(100);

  GNodeData* min = nullptr;
  Cmp cmp;

  // While the heap is not empty
  while (!heap.empty() || !r_set.empty()) {

    while (LIKELY(!heap.empty())) {
      T item = heap.top();
      GNodeData* item_data = &graph.getData(item.src);

      if (item_data->fixed || item_data->dist < item.dist) {
        // If we got a min, go do some work.
        heap.pop();
        if (!r_set.empty()) break;
        continue;
      }
      heap.pop();
      min = item_data;
      min->fixed = true;
      min->fixed_ts = high_resolution_clock::now();
      r_set.push_back(min);

      if (!(UNLIKELY(!heap.empty()) && heap.top().dist == min->dist)) break;
    }

    // Inner loop, go through the all the elements in R
    while(r_set.size() > 0) {
      GNodeData* z = r_set.back();
      r_set.pop_back();

      // Get all the vertices that have edges from z
      for (auto e : edgeRange(z->node)) {
        // If k vertex is not fixed, process the edge between z and k.
        auto k = &graph.getData(graph.getEdgeDst(e));
        if (k->fixed) continue;

        Dist z_k_dist = graph.getEdgeData(e);

        changed = false;
        if(z->dist + z_k_dist < k->dist) {
          k->dist = z->dist + z_k_dist;
          changed = true;
          if (k->dist < min->dist) min = k;
        }
        if (--k->pred <= 0) {
          k->fixed = true;
          k->fixed_ts = high_resolution_clock::now();
          r_set.push_back(k);
        } else if (changed) {
          pushWrap(heap,k->node,k->dist);
        }
      }

      if (r_set.empty() && !min->fixed && cmp(min->dist, heap.top().dist)){
        // We're done, but before we break, let's just check whether we have the new min in the q set
        // That is, if the heap is not empty and the current min is higher than the min in the q
        // set no point in pushing back to the heap, where it would have to bubble up.
        min->fixed = true;
        min->fixed_ts = high_resolution_clock::now();
        r_set.push_back(min);
      }
    }
  }
}

// template <typename GType,
//           typename T,
//           typename P,
//           typename R,
//           typename Graph = typename GType::Graph,
//           typename GNode = typename GType::GNode,
//           typename GNodeData = typename GType::Graph::node_data_type,
//           typename Traits = typename GType::Traits,
// 	  typename Cmp = std::less<typename Traits::Dist>>
// void serSP2Algo(Graph& graph, const GNode& source,
//                 const P& pushWrap, const R& edgeRange) {

//   using Heap = galois::MinHeap<T>;
//   using Dist = typename Traits::Dist;

//   bool changed;

//   Heap heap;
//   pushWrap(heap, source, 0);

//   galois::gstl::Vector<GNodeData*> r_set;
//   r_set.reserve(100);

//   GNodeData* min = nullptr;
//   Cmp cmp;

//   // While the heap is not empty
//   while (!heap.empty() || !r_set.empty()) {

//     while (LIKELY(!heap.empty())) {
//       T item = heap.top();
//       GNodeData* item_data = &graph.getData(item.src);

//       if (item_data->fixed || item_data->dist < item.dist) {
//         // If we got a min, go do some work.
//         heap.pop();
//         if (!r_set.empty()) break;
//         continue;
//       }
//       heap.pop();
//       min = item_data;
//       min->fixed = true;
//       min->fixed_ts = high_resolution_clock::now();
//       r_set.push_back(min);

//       if (heap.empty() || heap.top().dist != min->dist) break;
//     }

//     galois::gstl::Vector<GNodeData*> q_set1;
//     galois::gstl::Vector<GNodeData*> q_set2;
//     galois::gstl::Vector<GNodeData*>* q_set_current = &q_set1;
//     galois::gstl::Vector<GNodeData*>* q_set_aux = &q_set2;
//     // Inner loop, go through the all the elements in R
//     while(r_set.size() > 0) {
//       GNodeData* z = r_set.back();
//       r_set.pop_back();

//       // Get all the vertices that have edges from z
//       for (auto e : edgeRange(z->node)) {
//         auto k = &graph.getData(graph.getEdgeDst(e));
//         Dist z_k_dist = graph.getEdgeData(e);
//         // If k vertex is not fixed, process the edge between z and k.
//         if (k->fixed) continue;


//         changed = false;
//         if(z->dist + z_k_dist < k->dist) {
//           //std::cout << "Trying to relax to " <<  (z->dist + z_k_dist) << std::endl;
//           k->dist = z->dist + z_k_dist;
//           changed = true;
//           if (k->dist < min->dist) min = k;
//         }
//         if (--k->pred == 0 || k->dist <= (min->dist + k->min_in_weight)) {         
//           k->fixed = true;
//           k->fixed_ts = high_resolution_clock::now();
//           r_set.push_back(k);
//         } else if (changed) {
//           q_set_current->push_back(k);
//         }
//       }

//       if (r_set.empty()) {        
//         if (!min->fixed && min->dist <= heap.top().dist) {
//           // We're done, but before we break, let's just check whether we have the new min in the q set
//           // That is, if the heap is not empty and the current min is higher than the min in the q
//           // set no point in pushing back to the heap, where it would have to bubble up.
//           min->fixed = true;
//           min->fixed_ts = high_resolution_clock::now();
//           r_set.push_back(min);
//         }
//         for (auto& q : *q_set_current) {
//           if (q->fixed) {
//             r_set.push_back(q);
//           } else {
//             q_set_aux->push_back(q);
//           }
//         }
//         if (r_set.empty()) {
//           for (auto& q: *q_set_aux) {
//             pushWrap(heap, q->node, q->dist);
//           }
//         } else {
//           q_set_current->clear();
//           auto temp =  q_set_current;
//           q_set_current = q_set_aux;
//           q_set_aux = temp;          
//         }
//       }
//     }
//   }
// }

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits,
	  typename Cmp = std::less<typename Traits::Dist>>
void serSP2Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange) {

  using Heap = galois::MinHeap<T>;
  using Dist = typename Traits::Dist;

  bool changed;

  Heap heap;
  pushWrap(heap, source, 0);

  galois::gstl::Vector<GNodeData*> r_set;
  r_set.reserve(100);

  GNodeData* min = nullptr;
  Cmp cmp;

  // While the heap is not empty
  while (!heap.empty() || !r_set.empty()) {

    while (LIKELY(!heap.empty())) {
      T item = heap.top();
      GNodeData* item_data = &graph.getData(item.src);

       heap.pop();
      if (item_data->fixed || item_data->dist < item.dist) {
        // If we got a min, go do some work.
        if (!r_set.empty()) break;
        continue;
      }
      min = item_data;
      min->fixed = true;
      min->fixed_ts = high_resolution_clock::now();
      r_set.push_back(min);

      if (heap.empty() || heap.top().dist != min->dist) break;
    }

    // Inner loop, go through the all the elements in R
    while(r_set.size() > 0) {
      GNodeData* z = r_set.back();
      //std::cout << "Processing node " << z->node << std::endl;
      r_set.pop_back();

      // Get all the vertices that have edges from z
      for (auto e : edgeRange(z->node)) {
        auto k = &graph.getData(graph.getEdgeDst(e));
        Dist z_k_dist = graph.getEdgeData(e);
        // If k vertex is not fixed, process the edge between z and k.
        if (k->fixed) continue;


        changed = false;
        if(z->dist + z_k_dist < k->dist) {
          k->dist = z->dist + z_k_dist;
          changed = true;
          if (k->dist < min->dist) min = k;
        }
        if (--k->pred == 0 || k->dist <= (min->dist + k->min_in_weight)) {         
          k->fixed = true;
          k->fixed_ts = high_resolution_clock::now();
          r_set.push_back(k);
        } else if (changed) {
          pushWrap(heap,k->node,k->dist);
        }
      }

      if (r_set.empty() && !min->fixed && cmp(min->dist, heap.top().dist)){
        // We're done, but before we break, let's just check whether we have the new min in the q set
        // That is, if the heap is not empty and the current min is higher than the min in the q
        // set no point in pushing back to the heap, where it would have to bubble up.
        min->fixed = true;
        min->fixed_ts = high_resolution_clock::now();
        r_set.push_back(min);
      }
    }
  }
}


template<typename GNode,
         typename Graph,
         typename Container,
         typename R>
inline int push_back_edges_of_node(const GNode& node,
                                   Graph& graph,
                                   Container& c,
                                   const R& edgeRange) {
  int edges_pushed = 0;
  for (auto& edge : edgeRange(node)) {
    c.push_back(std::make_pair(node, *edge));
    edges_pushed++;
  }
  return edges_pushed;
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
void parSP1VerticesAlgo(Graph& graph,
                        const GNode& source,
                        const P& pushWrap,
                        const R& edgeRange) {

  using Heap = galois::ThreadSafeMinHeap<T>;
  using WorkItem = T;
  using PSchunk = galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>;
  using Dist = typename Traits::Dist;

  constexpr galois::MethodFlag flag = galois::MethodFlag::UNPROTECTED;

  Heap heap;
  std::atomic<uint32_t> last_heap_pop_dist(0);
  std::atomic<int> r_set_counter(1);
  std::mutex mutex;
  auto& source_data = graph.getData(source);
  source_data.fixed = true;
  source_data.fixed_ts = high_resolution_clock::now();

  galois::for_each(galois::iterate({T{source, 0}}),
                   [&](WorkItem z_item, auto& ctx) {

                     const auto& z_data = graph.getData(z_item.src, flag);
                     //std::cout << "Processing node " << z_data.node << std::endl;

                     for (auto& z_edge : edgeRange(z_item.src)) {
                       auto k = graph.getEdgeDst(z_edge);
                       auto& k_data = graph.getData(k, flag);

                       //std::cout << "Processing edge " << *z_edge << " from " << z_item.src << " to dst : " << k_data.node << " Current dist: " << k_dist << " New dist: " << new_dist << std::endl;
                       if (k_data.fixed) continue;

                       auto z_k_dist = graph.getEdgeData(z_edge, flag);
                       auto new_dist = z_data.dist + z_k_dist;
                       unsigned int k_dist = k_data.dist;;
                       bool changed = false;

                       
                       while (new_dist < k_dist) {
                         //std::cout << "Trying to relax to " << new_dist << std::endl;
                         if (k_data.dist.compare_exchange_weak(k_dist, new_dist, std::memory_order_relaxed)) {
                           changed = true;
                           break;
                         }
                         k_dist = k_data.dist;
                         new_dist = z_data.dist + z_k_dist;                         
                       }                      

                       // If the k vertex is now fixed, push the edges to the r set, otherwise push k
                       // to the heap.
                       if (--k_data.pred == 0) {
                         k_data.fixed = true;
                         k_data.fixed_ts = high_resolution_clock::now();                         
                         //std::cout << "Fixed " << k_data.node << " to " << k_data.dist << " at: " << duration_cast<microseconds>(k_data.fixed_ts - start).count() << std::endl;
                         //std::cout << "SP1 cond " << k_data.pred << " SP2 cond" << (last_heap_pop_dist.load() + k_data.min_in_weight) << std::endl;
                         pushWrap(ctx, k, k_data.dist);
                         ++r_set_counter;
                       } else if (changed){
                         pushWrap(heap, k, k_data.dist);
                       }
                     }

                     
                     if (--r_set_counter == 0) {

                       uint32_t min_dist = 0;
                       uint32_t items_pushed = 0;
                       while(!heap.empty()) {
                         auto j = heap.top();
                         auto& j_data = graph.getData(j.src, flag);

                         if (j_data.fixed || j_data.dist < j.dist) {
                           heap.pop();
                           if (r_set_counter.load() > 0) break;
                           continue;
                         }

                         heap.pop();
                         last_heap_pop_dist = (unsigned int)j_data.dist;
                         min_dist = j_data.dist;
                         j_data.fixed = true;
			 j_data.fixed_ts = high_resolution_clock::now();
                         ++r_set_counter;
                         pushWrap(ctx, j.src, j.dist);
                         items_pushed++;

                         if (!(UNLIKELY(!heap.empty()) && heap.top().dist == min_dist)) break;
                       }
                     }
                   },
                   galois::wl<PSchunk>(),
                   galois::no_conflicts(),
                   galois::loopname("sssp_parsp2_inner"));
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
void parSP2VerticesAlgo(Graph& graph,
                        const GNode& source,
                        const P& pushWrap,
                        const R& edgeRange) {

  using Heap = galois::ThreadSafeMinHeap<T>;
  using WorkItem = T;
  using PSchunk = galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>;
  using Dist = typename Traits::Dist;

  constexpr galois::MethodFlag flag = galois::MethodFlag::UNPROTECTED;

  Heap heap;
  std::atomic<uint32_t> min_dist(0);
  std::atomic<int> r_set_counter(1);
  std::mutex mutex;
  auto& source_data = graph.getData(source);
  source_data.fixed = true;
  source_data.fixed_ts = high_resolution_clock::now();

  galois::for_each(galois::iterate({T{source, 0}}),
                   [&](WorkItem z_item, auto& ctx) {

                     const auto& z_data = graph.getData(z_item.src, flag);
                     //std::cout << "Processing node " << z_data.node << std::endl;

                     for (auto& z_edge : edgeRange(z_item.src)) {
                       auto k = graph.getEdgeDst(z_edge);
                       auto& k_data = graph.getData(k, flag);

                       //std::cout << "Processing edge " << *z_edge << " from " << z_item.src << " to dst : " << k_data.node << " Current dist: " << k_dist << " New dist: " << new_dist << std::endl;
                       if (k_data.fixed) continue;

                       auto z_k_dist = graph.getEdgeData(z_edge, flag);
                       auto new_dist = z_data.dist + z_k_dist;
                       unsigned int k_dist = k_data.dist;;
                       bool changed = false;
                       
                       while (new_dist < k_dist) {
                         //std::cout << "Trying to relax to " << new_dist << std::endl;
                         if (k_data.dist.compare_exchange_weak(k_dist, new_dist, std::memory_order_relaxed)) {
                           changed = true;
                           break;
                         }
                         k_dist = k_data.dist;
                         new_dist = z_data.dist + z_k_dist;                         
                       }                      

                       // If the k vertex is now fixed, push the edges to the r set, otherwise push k
                       // to the heap.
                       --k_data.pred;
                       if (k_dist <= min_dist + k_data.min_in_weight || k_data.pred <= 0) {
                         k_data.fixed = true;
                         k_data.fixed_ts = high_resolution_clock::now();                         
                         //std::cout << "Fixed " << k_data.node << " to " << k_data.dist << " at: " << duration_cast<microseconds>(k_data.fixed_ts - start).count() << std::endl;
                         //std::cout << "SP1 cond " << k_data.pred << " SP2 cond" << (last_heap_pop_dist.load() + k_data.min_in_weight) << std::endl;
                         ++r_set_counter;
                         pushWrap(ctx, k, k_data.dist);
                       } else if (changed){
                         pushWrap(heap, k, k_data.dist);
                       }
                     }

                     
                     if (--r_set_counter == 0) {
                       while(!heap.empty()) {
                         auto j = heap.top();
                         auto& j_data = graph.getData(j.src, flag);
                         heap.pop();
                         if (j_data.fixed || j_data.dist < j.dist) {
                           if (r_set_counter.load() > 0) break;
                           continue;
                         }

                         min_dist = (unsigned int)j_data.dist;
                         j_data.fixed = true;
			 j_data.fixed_ts = high_resolution_clock::now();
                         ++r_set_counter;
                         pushWrap(ctx, j.src, j.dist);

                         if (heap.empty() || heap.top().dist != min_dist) break;
                       }
                     }
                   },
                   galois::wl<PSchunk>(),
                   galois::no_conflicts(),
                   galois::loopname("sssp_parsp2_inner"));
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
          typename GNode = typename GType::GNode>
void reset_graph(Graph& graph, GNode& source){

  using SSSP = typename GType::Traits::SSSP;

  galois::do_all(galois::iterate(graph),
                 [&graph](GNode n) {
                   auto& data = graph.getData(n);
                   data.dist = SSSP::DIST_INFINITY;
		   data.pred = total_pred[data.node];
		   data.fixed = false;
                 });

  graph.getData(source).dist = 0;
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
                   data.min_in_weight = SSSP::DIST_INFINITY;
                   data.node = n;
                 });

  graph.getData(source).dist = 0;
  total_pred.resize(graph.size());
}

struct MatchPathSeparator {
  bool operator()( char ch ) const {
    return ch == '/';
  }
};

std::string basename( std::string const& pathname ) {
    return std::string(
        std::find_if( pathname.rbegin(), pathname.rend(),
                      MatchPathSeparator() ).base(),
        pathname.end() );
}

std::string removeExtension( std::string const& filename ) {
    std::string::const_reverse_iterator
                        pivot
            = std::find( filename.rbegin(), filename.rend(), '.' );
    return pivot == filename.rend()
        ? filename
        : std::string( filename.begin(), pivot.base() - 1 );
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void verify_and_report(Graph& graph, GNode& source, GNode& report, unsigned long exec_time) {
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

  acc::accumulator_set<double, acc::stats<acc::tag::count,
                                          acc::tag::mean,
                                          acc::tag::variance,
                                          acc::tag::max,
                                          acc::tag::min,
                                          acc::tag::median>> stats;

  for (auto& node : graph) {
    auto& node_data = graph.getData(node);
    if (node_data.fixed) stats(duration_cast<microseconds>(node_data.fixed_ts - start).count());
  }

  auto count =  acc::count(stats);
  auto mean =  acc::mean(stats);
  auto median = acc::median(stats);
  auto max =  acc::max(stats);
  auto min = acc::min(stats);
  auto variance =  acc::variance(stats);

  std::cout << std::setprecision(16)
            << "count:    " << count    << '\n'
            << "mean:     " << mean     << '\n'
            << "median:   " << median   << '\n'
            << "max:      " << max      << '\n' 
            << "min:      " << min      << '\n'
            << "variance: " << variance << '\n';

    //Auxiliary code for stat collection to CSV file
  std::string resultfile_name = "results.csv";
  std::string graph_name = removeExtension(basename(filename));
  std::ifstream in_data_stream(resultfile_name);
  bool print_header = false;
  if(in_data_stream.peek() == std::ifstream::traits_type::eof())
    print_header = true;
  in_data_stream.close();
  std::ofstream data_stream;
  data_stream.open(resultfile_name, std::ofstream::out | std::ofstream::app);
  if(print_header)
    data_stream << "Algo;Graph_name;Main_time;Count;Mean;Median;max;min;variance;";
  data_stream << std::endl << ALGO_NAMES[algo] << ';' << graph_name << ';' << exec_time << ';' << count << ';' << mean << ';' << median << ';' << max << ';' << min << ';' << variance; 
  data_stream.close();

}

int main(int argc, char** argv) {
  using ParGType = GraphType<GraphTypes::ParGraphType>;
  using SerGType = GraphType<GraphTypes::SerGraphType>;

  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  skipVerify = true;
  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;
  galois::StatTimer Tinit;
  galois::StatTimer Tmain;
  unsigned long Init_Time, Main_Time;

  switch (algo) {
    case deltaTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();
      
      Tmain.start();
      deltaStepAlgo<ParGType, ParGType::Traits::SrcEdgeTile>(
          graph, source,
          ParGType::Traits::SrcEdgeTilePushWrap{graph},
          ParGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    case deltaStep: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      deltaStepAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          ParGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    case serDeltaTile: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();
      
      Tmain.start();
      serDeltaAlgo<SerGType, SerGType::Traits::SrcEdgeTile>(
          graph, source,
          SerGType::Traits::SrcEdgeTilePushWrap{graph},
          SerGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case serDelta: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();
      
      start = high_resolution_clock::now();

      Tmain.start();
      serDeltaAlgo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          SerGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case dijkstraTile: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      dijkstraAlgo<SerGType, SerGType::Traits::SrcEdgeTile>(
          graph, source,
          SerGType::Traits::SrcEdgeTilePushWrap{graph},
          SerGType::Traits::TileRangeFn());
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case dijkstra: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      dijkstraAlgo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          SerGType::Traits::OutEdgeRangeFn{graph});
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case topo: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      topoAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    case topoTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      topoTileAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    case serSP1: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      serSP1Algo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case serSP2: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();

      Tmain.start();
      serSP2Algo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, Tmain.get());
      break;
    }
    case parSP1: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();      
      
      Tmain.start();
      parSP1VerticesAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    case parSP2: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors(graph, edgeRange);

      start = high_resolution_clock::now();      

      Tmain.start();
      parSP2VerticesAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      //skipVerify = true;
      verify_and_report<ParGType>(graph, source, report, Tmain.get());
      break;
    }
    default:
      std::abort();
  }
  return 0;
}
