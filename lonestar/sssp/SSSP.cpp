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

#include <fstream>
#include <iostream>
#include <unordered_set>


#ifndef NDEBUG
#define D(x) x
#else
#define D(x)
#endif

#define LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)


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
  deltaStepSP1,
  serDeltaTile,
  serDelta,
  dijkstraTile,
  dijkstra,
  topo,
  topoTile,
  serSP1,
  serSP2,
  parSP1,
  parSP2V,
  parSP2E
};

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep", "deltaStepSP1", "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo", "topoTile", "serSP1", "serSP2", "parSP1", "parSP2V", "parSP2E"};

static cll::opt<Algo>
    algo("algo", cll::desc("Choose an algorithm:"),
         cll::values(clEnumVal(deltaTile, "deltaTile"),
                     clEnumVal(deltaStep, "deltaStep"),
                     clEnumVal(deltaStepSP1, "deltaStepSP1"),
                     clEnumVal(serDeltaTile, "serDeltaTile"),
                     clEnumVal(serDelta, "serDelta"),
                     clEnumVal(dijkstraTile, "dijkstraTile"),
                     clEnumVal(dijkstra, "dijkstra"), clEnumVal(topo, "topo"),
                     clEnumVal(topoTile, "topoTile"),
                     clEnumVal(serSP1, "serSP1"),
                     clEnumVal(serSP2, "serSP2"),
                     clEnumVal(parSP1, "parSP1"),
                     clEnumVal(parSP2V, "parSP2V"),
                     clEnumVal(parSP2E, "parSP2E"),
                     clEnumValEnd),
         cll::init(deltaTile));

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

typedef std::chrono::steady_clock clock_type;
typedef std::chrono::time_point<clock_type> timestamp;

static timestamp profiling_start;

static std::atomic<uint64_t> fixed_nodes;
static std::atomic<uint64_t> visited_nodes;
static std::atomic<uint64_t> explored_nodes;

template<typename GraphTypes>
struct SerNodeData {
  typename GraphTypes::SerGraphType::GraphNode node;
  int pred;
  uint32_t dist;
  bool fixed;
  uint32_t min_in_weight;
  SerNodeData(): node(-1), pred(0), dist(0), fixed(false) {}

  inline void visit() {}
  inline void explore() {}
  inline void fix() {}
};

template<typename GraphTypes>
struct ParNodeData {
  typename GraphTypes::ParGraphType::GraphNode node;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;
  uint32_t min_in_weight;
  ParNodeData(): node(-1), pred(0), dist(0), fixed(false) {}

  inline void visit() {}
  inline void explore() {}
  inline void fix() {}
};

template<typename GraphTypes>
struct ProfilingNodeData {
  typename GraphTypes::ParGraphType::GraphNode node;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;
  uint32_t min_in_weight;

  std::atomic<uint64_t> num_visits;
  std::atomic<uint64_t> first_visit_usec;
  std::atomic<uint64_t> first_visit_idx;
  std::atomic<uint64_t> last_visit_usec;
  std::atomic<uint64_t> visits_after_fixed;
  std::atomic<uint64_t> visits_after_explored;
  std::atomic<uint64_t> first_explored_usec;
  std::atomic<uint64_t> first_explored_idx;
  std::atomic<uint64_t> num_explores;
  std::atomic<uint64_t> fixed_usec;
  std::atomic<uint64_t> fixed_idx;


  ProfilingNodeData(): node(-1),
                       pred(0),
                       dist(0),
                       fixed(false),
                       num_visits(0),
                       first_visit_usec(0),
                       first_visit_idx(0),
                       last_visit_usec(0),
                       visits_after_fixed(0),
                       visits_after_explored(0),
                       first_explored_usec(0),
                       first_explored_idx(0),
                       num_explores(0),
                       fixed_usec(0),
                       fixed_idx(0) {}

  inline void explore() {
    num_explores++;
    if (first_explored_usec == 0) {
      timestamp now = clock_type::now();
      first_explored_usec = (now - profiling_start).count();
      first_explored_idx = explored_nodes++;
    }
  }

  inline void fix() {
    timestamp now = clock_type::now();
    fixed = true;
    fixed_usec = (now - profiling_start).count();
    fixed_idx = fixed_nodes++;
  }

  inline void visit() {
    if (num_visits == 0) {
      timestamp now = clock_type::now();
      first_visit_usec = (now - profiling_start).count();
      first_visit_idx = visited_nodes++;
    }
    num_visits++;
    if (num_explores > 0) {
      visits_after_explored++;
    }
    if (fixed) {
      timestamp now = clock_type::now();
      visits_after_fixed++;
      last_visit_usec = (now - profiling_start).count();
    }
  }
};


struct EdgeData {
  EdgeData(uint32_t dist_) : dist(dist_) {}
  uint32_t dist;
};

// Edge data to use for profiling purposes
struct ProfilingEdgeData {
  uint32_t dist;
  ProfilingEdgeData(uint32_t dist_) : dist(dist_) {}
};

// struct GraphTypes {
//   typedef SerNodeData<GraphTypes> SerNodeDataType;
//   typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::SerNodeDataType, uint32_t>::
//   with_no_lockable<true>::type ::with_numa_alloc<true>::type SerGraphType;
//   typedef ParNodeData<GraphTypes> ParNodeDataType;
//   typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::ParNodeDataType, uint32_t>::
//   with_no_lockable<true>::type ::with_numa_alloc<true>::type ParGraphType;
// };


// PROFILING GRAPH TYPES
struct GraphTypes {
  typedef ProfilingNodeData<GraphTypes> SerNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::SerNodeDataType, ProfilingEdgeData, false, false, false, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type SerGraphType;
  typedef ProfilingNodeData<GraphTypes> ParNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::ParNodeDataType, ProfilingEdgeData, false, false, false, uint32_t>::
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


template<typename Graph, typename HeapEntryType, typename NodeDataType>
struct GarbageCollectFixedNodes {
  GarbageCollectFixedNodes(Graph& graph_) : graph(graph_) {}
  Graph& graph;
  inline bool operator()(HeapEntryType& item) const {
    return graph.getData(item.src).fixed;
  }
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
                       Dist ew            = graph.getEdgeData(ii, flag).dist;
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

        const auto newDist = item.dist + graph.getEdgeData(e).dist;

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
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits>
void dijkstraAlgo(Graph& graph, const GNode& source, const P& pushWrap,
                  const R& edgeRange) {

  galois::Timer timer;
  timer.start();

  //using WL = galois::MinHeap<T>;
  using WL = galois::GarbageCollectingMinHeap<T, GarbageCollectFixedNodes<Graph, T, GNodeData>>;

  graph.getData(source).dist = 0;

  profiling_start = clock_type::now();
  explored_nodes = 0;
  fixed_nodes = 0;


  GarbageCollectFixedNodes<Graph, T, GNodeData> gc(graph);
  WL wl(gc);

  pushWrap(wl, source, 0);

  size_t iter = 0;
  size_t inner_iter = 0;
  size_t heap_pushes =0;
  double average_heap_size = 0;
  size_t duplicated_items = 0;

  while (!wl.empty()) {
    ++iter;
    average_heap_size += wl.size();

    T item = wl.top();
    wl.pop();

    auto& item_data = graph.getData(item.src);
    if (item_data.dist < item.dist) {
      duplicated_items++;
      continue;
    }

    item_data.fix();
    item_data.explore();

    for (auto e : edgeRange(item)) {

      inner_iter++;
      GNode dst   = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst);
      ddata.visit();

      const auto newDist = item.dist + graph.getEdgeData(e).dist;

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


// template <typename GType,
//           typename T,
//           typename P,
//           typename R,
//           typename Graph = typename GType::Graph,
//           typename GNode = typename GType::GNode,
//           typename GNodeData = typename GType::Graph::node_data_type,
//           typename Traits = typename GType::Traits>
// void deltaStepAlgoSP1(Graph& graph,
//                       GNode source,
//                       const P& pushWrap,
//                       const R& edgeRange) {

//   using Dist                 = typename Traits::Dist;
//   using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

//   //! [reducible for self-defined stats]
//   galois::GAccumulator<size_t> BadWork;
//   //! [reducible for self-defined stats]
//   galois::GAccumulator<size_t> WLEmptyWork;

//   namespace gwl = galois::worklists;

//   using PSchunk = gwl::PerSocketChunkFIFO<CHUNK_SIZE>;
//   using OBIM = galois::worklists::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

//   auto& sdata = graph.getData(source);
//   sdata.dist = 0;
//   sdata.fixed = true;

//   std::atomic<Dist> min(0);

//   galois::InsertBag<T> initBag;
//   pushWrap(initBag, source, 0, "parallel");

//   galois::for_each(galois::iterate(initBag),
//                    [&](const T& item, auto& ctx) {
//                      constexpr galois::MethodFlag flag =
//                          galois::MethodFlag::UNPROTECTED;
//                      auto& sdata = graph.getData(item.src, flag);

//                      if (sdata.dist < item.dist || (sdata.fixed && sdata.explored)) {
//                        if (TRACK_WORK) WLEmptyWork += 1;
//                        return;
//                      }

//                      thread_local std::vector<GNodeData*> r_set(100);
//                      r_set.push_back(&sdata);

//                      while (!r_set.empty()) {
//                        auto elem = r_set.back();
//                        r_set.pop_back();

//                        for (auto ii : edgeRange(elem)) {

//                          GNode dst = graph.getEdgeDst(ii);
//                          auto& ddist = graph.getData(dst, flag);
//                          Dist ew = graph.getEdgeData(ii, flag);
//                          const Dist newDist = elem.dist + ew;

//                          if (ddist.fixed) continue;

//                          bool changed = false;
//                          while (true) {
//                            Dist oldDist = ddist.dist;

//                            if (oldDist <= newDist) break;

//                            if (ddist.dist.compare_exchange_weak(
//                                    oldDist, newDist, std::memory_order_relaxed)) {
//                              changed = true;
//                              break;
//                            }
//                          }

//                          if (elem.fixed && --ddist.pred <= 0) {
//                            ddist.fixed = true;
//                            if (changed) {
//                              r_set.push_back(ddist);
//                              continue;
//                            }
//                          }

//                          if (changed) {
//                            pushWrap(ctx, dst, newDist);
//                          }
//                        }

//                        if (sdata.fixed) sdata.explored = true;
//                      }

//                    },
//                    galois::wl<OBIM>(UpdateRequestIndexer{stepShift}),
//                    galois::no_conflicts(), galois::loopname("SSSP"));

//   if (TRACK_WORK) {
//     //! [report self-defined stats]
//     galois::runtime::reportStat_Single("SSSP", "BadWork", BadWork.reduce());
//     //! [report self-defined stats]
//     galois::runtime::reportStat_Single("SSSP", "WLEmptyWork",
//                                        WLEmptyWork.reduce());
//   }
// }

namespace detail{
template<typename GType, typename Edge>
struct EdgeDist {
  auto operator()(typename GType::Graph& graph,
                  typename GType::Graph::edge_iterator& edge) {
    return graph.getEdgeData(edge);
  }
};

template<typename GType>
struct EdgeDist<GType, ProfilingEdgeData> {
  auto operator()(typename GType::Graph& graph,
                  typename GType::Graph::edge_iterator& edge) {
    return graph.getEdgeData(edge).dist;
  }
};
}

template<typename GType, typename R>
void calc_graph_predecessors(typename GType::Graph& graph, R& edgeRange) {
  detail::EdgeDist<GType, typename GType::Graph::edge_data_type> edgeDist;
  // Fill the pred array
  for (auto vertex : graph) {
    for (auto edge : edgeRange(vertex)) {
      auto& k_data = graph.getData(graph.getEdgeDst(edge));
      k_data.pred++;
      k_data.min_in_weight = std::min(k_data.min_in_weight, edgeDist(graph, edge));
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
      j_data.fix();
      j_data.explore();
      GNode z = j.src;

      // Inner loop, go through the all the elements in R
      while(true) {
        auto& z_data = graph.getData(z);
        // Get all the vertices that have edges from z
        for (auto e : edgeRange(z)) {
          inner_iter++;
          GNode k = graph.getEdgeDst(e);
          auto& k_data = graph.getData(k);
          k_data.visit();
          // If k vertex is not fixed, process the edge between z and k.
          if (!k_data.fixed) {
            auto z_k_dist = graph.getEdgeData(e).dist;
            k_data.pred--;
            if (k_data.pred <= 0) {
              k_data.fix();
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

template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits>
void serSP2Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange) {

  //using Heap = galois::MinHeap<T>;

  using Heap = galois::GarbageCollectingMinHeap<T, GarbageCollectFixedNodes<Graph, T, GNodeData>>;

  using Dist = typename Traits::Dist;
  // The set of nodes which have been fixed but not explored

  size_t nodes_fixed = 0;

  galois::gstl::Vector<GNodeData*> r_set;
  r_set.reserve(100);
  std::vector<std::pair<GNodeData*, Dist>> q_set;
  q_set.reserve(100);

  profiling_start = clock_type::now();
  explored_nodes = 0;
  fixed_nodes = 0;

  GarbageCollectFixedNodes<Graph, T, GNodeData> gc(graph);

  Heap heap(gc);
  pushWrap(heap, source, 0);

  GNodeData* min = nullptr;

  // While the heap is not empty
  while (LIKELY(nodes_fixed < graph.size()) && (!heap.empty() || !r_set.empty())) {

    if (LIKELY(!heap.empty())) {
      T item = heap.top();
      GNodeData* item_data = &graph.getData(item.src);

      if (item_data->fixed || item_data->dist < item.dist) {
        heap.pop();
        // If we got a min, go do some work.
        if (!r_set.empty()) goto mainloop;
        continue;
      }

      heap.pop();
      min = item_data;
      min->fix();
      r_set.push_back(min);
    }

    if (LIKELY(!heap.empty()) && heap.top().dist == min->dist) {
      continue;
    }

    mainloop:
    // Inner loop, go through the all the elements in R
    while(r_set.size() > 0) {
      GNodeData* z = r_set.back();
      r_set.pop_back();
      z->explore();

      // Get all the vertices that have edges from z
      for (auto e : edgeRange(z->node)) {
        auto k = &graph.getData(graph.getEdgeDst(e));
        k->visit();
        // If k vertex is not fixed, process the edge between z and k.
        if (k->fixed) continue;
        auto z_k_dist = graph.getEdgeData(e).dist;

        bool changed = false;
        if (k->dist > z->dist + z_k_dist) {
          k->dist = z->dist + z_k_dist;
          changed = true;
          if (k->dist < min->dist) min = k;
        }

        if (--k->pred <= 0 || k->dist <= (min->dist + k->min_in_weight)) {
          k->fix();
          r_set.push_back(k);
        } else if (changed) {
          q_set.push_back({k, k->dist});
        }
      }

      if (r_set.empty() && !min->fixed && min->dist <= heap.top().dist) {
        // We're done, but before we break, let's just check whether we have the new min in the q set
        // That is, if the heap is not empty and the current min is higher than the min in the q
        // set no point in pushing back to the heap, where it would have to bubble up.
        min->fix();
        r_set.push_back(min);
      }
    }

    for (auto& item : q_set) {
      if (item.first->dist == item.second && !item.first->fixed) {
        pushWrap(heap, item.first->node, item.first->dist);
      }
    }
    q_set.clear();
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
void parSP1Algo(Graph& graph,
                const GNode& source,
                const P& pushWrap,
                const R& edgeRange) {

  using Heap = galois::ThreadSafeMinHeap<T>;
  using WorkItem = std::pair<GNode, GEdge>;
  using PSchunk = galois::worklists::PerThreadChunkLIFO<CHUNK_SIZE>;
  using Dist = typename Traits::Dist;

  constexpr galois::MethodFlag flag = galois::MethodFlag::UNPROTECTED;

  Heap heap;
  pushWrap(heap, source, 0, "parallel");

  // A step in an unbundled r set which looks at a single edge
  // from a node in the r_set.
  auto r_set_edge_loop = [&](const WorkItem&  z_edge, auto& ctx) {
    auto k = graph.getEdgeDst(z_edge.second);
    auto& k_data = graph.getData(k, flag);

    if (k_data.fixed) return;

    // If k vertex is not fixed, process the edge between z and k.
    k_data.pred--;
    if (k_data.pred <= 0) {
      k_data.fixed = true;
      push_back_edges_of_node(k, graph, ctx, edgeRange);
    }

    // Note that if another thread reduces the distance in z
    // while we're in the while loop the algorithm is still
    // correct (if we change k because z is closer, then if
    // z get even closer it's still valid to change k)
    const auto& z_data = graph.getData(z_edge.first, flag);
    auto z_k_data = graph.getEdgeData(z_edge.second, flag);
    Dist new_dist = z_data.dist + z_k_data.dist;

    while (true) {
      Dist old_dist = k_data.dist;
      if (new_dist > old_dist) break;

      if (k_data.dist.compare_exchange_weak(old_dist, new_dist, std::memory_order_relaxed)) {
        if (!k_data.fixed) pushWrap(heap, k_data.node, k_data.dist);
        break;
      }
    }
  };


  size_t outer_iter = 0;
  size_t duplicated_items = 0;
  size_t average_heap_size = 0;
  size_t heap_pushes = 0;


  // While the heap is not empty
  while (!heap.empty()) {
    // Get the min element from the heap.
    T j = heap.pop();
    outer_iter++;
    average_heap_size += heap.size();
    auto& j_data = graph.getData(j.src, flag);

    if (j_data.dist < j.dist) {
      duplicated_items++;
      continue;
    }

    // If the min element is not fixed
    if (!j_data.fixed) {
      // Set the element to fixed
      j_data.fixed = true;
      galois::gstl::Vector<WorkItem> initial;
      push_back_edges_of_node(j.src, graph, initial, edgeRange);
      galois::for_each(galois::iterate(initial),
                       r_set_edge_loop,
                       galois::wl<PSchunk>(),
                       galois::no_conflicts(),
                       galois::loopname("sssp_parsp1_inner"));
    }
  }

  galois::runtime::reportStat_Single("SSSP-parSP1", "Outer loop iterations", outer_iter);
  galois::runtime::reportStat_Single("SSSP-parSP1", "Heap pushes", heap_pushes);
  galois::runtime::reportStat_Single("SSSP-parSP1", "Duplicated heap pushes", duplicated_items);
  galois::runtime::reportStat_Single("SSSP-parSP1", "Average heap size", (average_heap_size * 1.0) / outer_iter);
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
  std::atomic<uint32_t> last_heap_pop_dist(0);
  std::atomic<int> r_set_counter(1);

  galois::for_each(galois::iterate({T{source, 0}}),
                   [&](WorkItem z_item, auto& ctx) {

                     const auto& z_data = graph.getData(z_item.src, flag);

                     for (auto& z_edge : edgeRange(z_item.src)) {
                       auto k = graph.getEdgeDst(z_edge);
                       auto& k_data = graph.getData(k, flag);

                       if (k_data.fixed) continue;

                       auto z_k_data = graph.getEdgeData(z_edge, flag);
                       Dist k_dist = k_data.dist;
                       Dist new_k_dist = z_data.dist + z_k_data.dist;
                       bool changed = false;

                       do {
                         if (new_k_dist > k_dist) break;
                         if (k_data.dist.compare_exchange_weak(k_dist, new_k_dist, std::memory_order_relaxed)) {
                           k_dist = new_k_dist;
                           changed = true;
                         }
                         k_dist = k_data.dist;
                       } while(UNLIKELY(!changed));

                       // If the k vertex is now fixed, push the edges to the r set, otherwise push k
                       // to the heap.
                       if (--k_data.pred <= 0 || k_dist <= last_heap_pop_dist.load() + k_data.min_in_weight) {
                         k_data.fixed = true;
                         ++r_set_counter;
                         pushWrap(ctx, k, k_data.dist);
                       } else if (changed){
                         pushWrap(heap, k, k_data.dist);
                       }
                     }

                     uint32_t min_dist;
                     if (--r_set_counter == 0) {

                       uint32_t items_pushed = 0;
                       while(!heap.empty()) {
                         auto j = heap.top();
                         auto& j_data = graph.getData(j.src, flag);

                         if (j_data.dist < j.dist) {
                           heap.pop();
                           continue;
                         }

                         if (j_data.fixed) {
                           heap.pop();
                           continue;
                         }

                         if (items_pushed == 0) {
                           j = heap.pop();
                           last_heap_pop_dist = j.dist;
                           j_data.fixed = true;
                           min_dist = j_data.dist;
                           items_pushed++;
                           ++r_set_counter;
                           pushWrap(ctx, j.src, j.dist);
                         // Push all elements with the same weights
                         } else if (j_data.dist == min_dist) {
                           j = heap.pop();
                           j_data.fixed = true;
                           items_pushed++;
                           ++r_set_counter;
                           pushWrap(ctx, j.src, j.dist);
                         } else break;
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
void parSP2EdgesAlgo(Graph& graph,
                     const GNode& source,
                     const P& pushWrap,
                     const R& edgeRange) {

  using Heap = galois::ThreadSafeMinHeap<T>;
  using PSchunk = galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>;
  using Dist = typename Traits::Dist;

  constexpr galois::MethodFlag flag = galois::MethodFlag::UNPROTECTED;

  struct WorkItem {
    GNode src;
    Dist dist;
    GNode dst;
  };

  Heap heap;
  std::atomic<uint32_t> last_heap_pop_dist(0);
  std::atomic<int> r_set_counter(0);
  std::vector<WorkItem> initial;
  for (const auto& e: edgeRange(source)) {
    ++r_set_counter;
    initial.push_back(WorkItem{source, graph.getEdgeData(e).dist, graph.getEdgeDst(e)});
  }

  galois::for_each(galois::iterate(initial),
                   [&](WorkItem z_item, auto& ctx) {

                     auto& k_data = graph.getData(z_item.dst, flag);

                     if (!k_data.fixed) {
                       const auto& z_data = graph.getData(z_item.src, flag);
                       Dist k_dist = k_data.dist;
                       Dist new_k_dist = z_data.dist + z_item.dist;
                       bool changed = false;

                       while(new_k_dist < k_dist) {
                         if (k_data.dist.compare_exchange_weak(k_dist, new_k_dist, std::memory_order_relaxed)) {
                           k_dist = new_k_dist;
                           changed = true;
                           break;
                         }
                         k_dist = k_data.dist;
                       }

                       // If the k vertex is now fixed, push the edges to the r set, otherwise push k
                       // to the heap.
                       if (--k_data.pred <= 0 || k_dist <= last_heap_pop_dist.load() + k_data.min_in_weight) {
                         k_data.fixed = true;
                         for (auto& e: edgeRange(z_item.dst)) {
                           ++r_set_counter;
                           ctx.push_back(WorkItem{z_item.dst, graph.getEdgeData(e, flag).dist, graph.getEdgeDst(e)});
                         }
                       } else if (changed){
                         pushWrap(heap, z_item.dst, k_data.dist);
                       }
                     }

                     uint32_t min_dist;
                     if (--r_set_counter == 0) {
                       bool pushed_item = false;
                       while(!heap.empty()) {
                         auto j = heap.top();
                         auto& j_data = graph.getData(j.src, flag);

                         if (j_data.fixed || j_data.dist < j.dist) {
                           heap.pop();
                           continue;
                         }

                         if (!pushed_item) {
                           j = heap.pop();
                           last_heap_pop_dist = j.dist;
                           j_data.fixed = true;
                           min_dist = j_data.dist;
                           pushed_item = true;
                           for (auto& e: edgeRange(j.src)) {
                             ++r_set_counter;
                             ctx.push_back(WorkItem{j.src, graph.getEdgeData(e, flag).dist, graph.getEdgeDst(e)});
                           }
                         // Push all elements with the same weights
                         } else if (j_data.dist == min_dist) {
                           j = heap.pop();
                           j_data.fixed = true;
                           for (auto& e: edgeRange(j.src)) {
                             ++r_set_counter;
                             ctx.push_back(WorkItem{j.src, graph.getEdgeData(e, flag).dist, graph.getEdgeDst(e)});
                           }
                         } else break;
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
                         const auto newDist = sdata.dist + graph.getEdgeData(e).dist;
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
                         const auto newDist = sdata.dist + graph.getEdgeData(e).dist;
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
                   data.min_in_weight = SSSP::DIST_INFINITY;
                   data.node = n;
                 });

  graph.getData(source).dist = 0;

  profiling_start = clock_type::now();
  explored_nodes = 0;
  fixed_nodes = 0;
}

struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '/';
    }
};

std::string
basename( std::string const& pathname )
{
    return std::string(
        std::find_if( pathname.rbegin(), pathname.rend(),
                      MatchPathSeparator() ).base(),
        pathname.end() );
}

std::string
removeExtension( std::string const& filename )
{
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
void verify_and_report(Graph& graph, GNode& source, GNode& report, std::string algo_name) {
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

  std::string graph_name = removeExtension(basename(filename));

  std::string prof_filename = "per_node_profile_" + algo_name + "_" + graph_name + ".csv";

  std::cout << "Outputting profile to file: " << prof_filename;

  // For profiling graphs, output the metrics in csv form
  std::ofstream profile;
  profile.open(prof_filename);
  profile << "node_idx;pred;dist;fixed;min_in_weight;num_visits;first_visit_usec;first_visit_idx";
  profile << "last_visit_usec;visits_after_fixed;visits_after_explored;first_explored_usec;first_explored_idx;";
  profile << "num_explores;fixed_usec;fixed_idx;" << std::endl;

  for (auto& node : graph) {
    auto& d = graph.getData(node);
    profile << d.node << ";" << d.pred << ";" << d.dist << ";" << d.fixed << ";" << d.min_in_weight;
    profile << d.num_visits << ";" << d.first_visit_usec << ";" << d.first_visit_idx << ";";
    profile << d.last_visit_usec << ";" << d.visits_after_fixed << ";" << d.visits_after_explored << ";";
    profile << d.first_explored_usec << ";" << d.first_explored_idx << ";";
    profile << d.num_explores << ";" << d.fixed_usec << ";" << d.fixed_idx << ";" << std::endl;
  }

  profile.close();
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

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case deltaStep: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<ParGType>(graph, edgeRange);

      Tmain.start();
      deltaStepAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    // case deltaStepSP1: {
    //   ParGType::Graph graph;
    //   ParGType::GNode source, report;
    //   init_graph<ParGType>(graph, source, report);
    //   auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
    //   calc_graph_predecessors<ParGType>(graph, edgeRange);

    //   Tmain.start();
    //   deltaStepAlgoSP1<ParGType, ParGType::Traits::UpdateRequest>(
    //       graph, source,
    //       ParGType::Traits::ReqPushWrap(),
    //       edgeRange);
    //   Tmain.stop();

    //   verify_and_report<ParGType>(graph, source, report);
    //   break;
    // }
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

      verify_and_report<SerGType>(graph, source, report, ALGO_NAMES[algo]);
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

      verify_and_report<SerGType>(graph, source, report, ALGO_NAMES[algo]);
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

      verify_and_report<SerGType>(graph, source, report, ALGO_NAMES[algo]);
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

      verify_and_report<SerGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case topo: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      topoAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case topoTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);

      Tmain.start();
      topoTileAlgo<ParGType>(graph, source);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case serSP1: {
      SerGType::Graph graph;
      SerGType::GNode source, report;
      init_graph<SerGType>(graph, source, report);
      auto edgeRange = SerGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<SerGType>(graph, edgeRange);

      Tmain.start();
      serSP1Algo<SerGType, SerGType::Traits::UpdateRequest>(
          graph, source,
          SerGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<SerGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case serSP2: {
      using GType = SerGType;
      GType::Graph graph;
      GType::GNode source, report;
      init_graph<GType>(graph, source, report);
      auto edgeRange = GType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<GType>(graph, edgeRange);

      Tmain.start();
      serSP2Algo<GType, GType::Traits::UpdateRequest>(
          graph, source,
          GType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();
      verify_and_report<GType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case parSP1: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<ParGType>(graph, edgeRange);

      Tmain.start();
      parSP1Algo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case parSP2V: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<ParGType>(graph, edgeRange);

      Tmain.start();
      parSP2VerticesAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    case parSP2E: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      init_graph<ParGType>(graph, source, report);
      auto edgeRange = ParGType::Traits::OutEdgeRangeFn{graph};
      calc_graph_predecessors<ParGType>(graph, edgeRange);

      Tmain.start();
      parSP2EdgesAlgo<ParGType, ParGType::Traits::UpdateRequest>(
          graph, source,
          ParGType::Traits::ReqPushWrap(),
          edgeRange);
      Tmain.stop();

      verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      break;
    }
    default:
      std::abort();
  }

  return 0;
}
