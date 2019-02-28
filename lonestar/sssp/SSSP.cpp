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
#include "time.h"

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
static cll::opt<unsigned int>
    APSP("APSP",
              cll::desc("Run All Pairs Shortest Path version of the Algorithm if supported"),
              cll::init(0));
static cll::opt<unsigned int>
    VERIFY("verify",
              cll::desc("Number of sources to verify in All Pairs Shortest Path"),
              cll::init(2));


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
static std::atomic<uint64_t> heap_operations_global;
static std::atomic<uint64_t> heap_pushes_global;

enum VertexSource {
  R_SET,
  HEAP
};

struct VertexContext {
  VertexSource vertex_source;
  size_t set_size;
};

template<typename GraphTypes>
struct SerNodeData {
  typename GraphTypes::SerGraphType::GraphNode node;
  int pred;
  uint32_t dist;
  bool fixed;
  uint32_t min_in_weight;
  SerNodeData(): node(-1), pred(0), dist(0), fixed(false) {}

  inline void visit() {}
  inline void explore(VertexContext context = VertexContext()) {}
  inline void fix(VertexContext context = VertexContext()) {}
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
  inline void explore(VertexContext context = VertexContext()) {}
  inline void fix(VertexContext context = VertexContext()) {}
};

template<typename GraphTypes>
struct APSPNodeData {
  typename GraphTypes::APSPGraphType::GraphNode node;
  std::vector<int> pred;
  std::vector<uint32_t> dist;
  std::vector<bool> fixed;
  uint32_t min_in_weight;
  uint32_t min_out_weight;
};

constexpr static const uint32_t DIST_INF = std::numeric_limits<uint32_t>::max() / 2 - 1;

template<typename GraphTypes>
struct ProfilingNodeData {
  typename GraphTypes::ParGraphType::GraphNode node;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;
  uint32_t min_in_weight;
  uint32_t min_out_weight;

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
  VertexContext ctx_explored;
  VertexContext ctx_fixed;


  ProfilingNodeData(): node(-1) {
    reset();
  }

  inline void explore(VertexContext context = VertexContext()) {
    if (first_explored_usec == 0) {
      num_explores++;
      ctx_explored = context;
      timestamp now = clock_type::now();
      first_explored_usec = (now - profiling_start).count();
      first_explored_idx = explored_nodes++;
    }
  }

  inline void fix(VertexContext context = VertexContext()) {
    if (fixed) return;
    timestamp now = clock_type::now();
    fixed = true;
    fixed_usec = (now - profiling_start).count();
    fixed_idx = fixed_nodes++;
    ctx_fixed = context;
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

  void reset(bool reset_dist = true) {
    pred = 0;
    if (reset_dist) dist = DIST_INF;
    fixed = false;
    min_in_weight = DIST_INF;
    min_out_weight = DIST_INF;
    num_visits = 0;
    first_visit_usec = 0;
    first_visit_idx = 0;
    last_visit_usec = 0;
    visits_after_fixed = 0;
    visits_after_explored = 0;
    first_explored_usec = 0;
    first_explored_idx = 0;
    num_explores = 0;
    fixed_usec = 0;
    fixed_idx = 0;
    ctx_explored = VertexContext();
    ctx_fixed = VertexContext();
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
//   typedef APSPNodeData<GraphTypes> APSPNodeDataType;
//   typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::APSPNodeDataType, uint32_t>::
//   with_no_lockable<true>::type ::with_numa_alloc<true>::type APSPGraphType;
// };


// PROFILING GRAPH TYPES
struct GraphTypes {
  typedef ProfilingNodeData<GraphTypes> SerNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::SerNodeDataType, ProfilingEdgeData, false, false, false, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type SerGraphType;
  typedef ProfilingNodeData<GraphTypes> ParNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::ParNodeDataType, ProfilingEdgeData, false, false, false, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type ParGraphType;
  typedef APSPNodeData<GraphTypes> APSPNodeDataType;
  typedef galois::graphs::LC_CSR_Graph<typename GraphTypes::APSPNodeDataType, ProfilingEdgeData, false, false, false, uint32_t>::
  with_no_lockable<true>::type ::with_numa_alloc<true>::type APSPGraphType;
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
                  const P& pushWrap, const R& edgeRange, uint32_t index = 0) {

  using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

  SerialBucketWL<T, UpdateRequestIndexer> wl(UpdateRequestIndexer{stepShift});
  graph.getData(source).dist.at(index) = 0;

  pushWrap(wl, source, 0);

  while (!wl.empty()) {

    auto& curr = wl.minBucket();

    while (!curr.empty()) {
      auto item = curr.front();
      curr.pop_front();

      if (graph.getData(item.src).dist.at(index) < item.dist) {
       continue;
      }

      for (auto e : edgeRange(item)) {
        GNode dst   = graph.getEdgeDst(e);
        auto& ddata = graph.getData(dst);

        const auto newDist = item.dist + graph.getEdgeData(e).dist;

        if (newDist < ddata.dist.at(index)) {
          ddata.dist.at(index) = newDist;
          pushWrap(wl, dst, newDist);
        }
      }
    }

    wl.goToNextBucket();
  }

  if (!wl.allEmpty()) {
    std::abort();
  }
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
                  const R& edgeRange, uint32_t index = 0) {

  using WL = galois::MinHeap<T>;
  graph.getData(source).dist.at(index) = 0;
  uint32_t heap_operations = 0;
  uint32_t heap_pushes = 0;

  WL wl;
  pushWrap(wl, source, 0);
  heap_operations++;
  heap_pushes++;

  while (!wl.empty()) {

    T item = wl.top();
    wl.pop();
    heap_operations++;

    auto& item_data = graph.getData(item.src);
    if (item_data.dist.at(index) < item.dist) {
      continue;
    }

    for (auto e : edgeRange(item)) {

      GNode dst   = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst);

      const auto newDist = item.dist + graph.getEdgeData(e).dist;

      if (newDist < ddata.dist.at(index)) {
        ddata.dist.at(index) = newDist;
        pushWrap(wl, dst, newDist);
	heap_operations++;
	heap_pushes++;
      }
    }
  }
  heap_operations_global += heap_operations;
  heap_pushes_global += heap_pushes;
}


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
    auto& v_data = graph.getData(vertex);
    for (auto edge : edgeRange(vertex)) {
      auto edge_dist = edgeDist(graph, edge);
      auto& k_data = graph.getData(graph.getEdgeDst(edge));
      k_data.pred++;
      k_data.min_in_weight = std::min(k_data.min_in_weight, edge_dist);
      if (v_data.min_out_weight > edge_dist) {
        v_data.min_out_weight = edge_dist;
      }
    }
  }
}

template<typename GType, typename R>
void calc_graph_predecessorsAPSP(typename GType::Graph& graph, R& edgeRange) {
  detail::EdgeDist<GType, typename GType::Graph::edge_data_type> edgeDist;
  // Fill the pred array
  for (auto vertex : graph) {
    auto& v_data = graph.getData(vertex);
    for (auto edge : edgeRange(vertex)) {
      auto edge_dist = edgeDist(graph, edge);
      auto& k_data = graph.getData(graph.getEdgeDst(edge));
      k_data.pred.at(0)++;
      k_data.min_in_weight = std::min(k_data.min_in_weight, edge_dist);
      if (v_data.min_out_weight > edge_dist) {
        v_data.min_out_weight = edge_dist;
      }
    }
  }
  if(APSP){
    for (auto vertex : graph) {
      auto& v_data = graph.getData(vertex);
      for(auto index : graph) {
        v_data.pred.at(index) = v_data.pred.at(0);
      }
    }
  }
}


template <typename GType,
          typename T,
          typename P,
          typename R,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename GNodeData = typename GType::Graph::node_data_type,
          typename Traits = typename GType::Traits>
void serSP1Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange, uint32_t index = 0) {
  using Heap = galois::MinHeap<T>;

  Heap heap;
  graph.getData(source).dist.at(index) = 0;
  pushWrap(heap, source, 0);

  // The set of nodes which have been fixed but not explored
  galois::gstl::Vector<GNode> r_set;

  // While the heap is not empty
  while (!heap.empty()) {
    // Get the min element from the heap.
    T j = heap.pop();

    auto& j_data = graph.getData(j.src);

    if (j_data.dist.at(index) < j.dist) {
      continue;
    }

    // If the min element is not fixed
    if (!j_data.fixed.at(index)) {
      // Set the element to fixed
      j_data.fixed.at(index) = true;
      GNode z = j.src;

      // Inner loop, go through the all the elements in R
      while(true) {
        auto& z_data = graph.getData(z);
        // Get all the vertices that have edges from z
        for (auto e : edgeRange(z)) {
          GNode k = graph.getEdgeDst(e);
          auto& k_data = graph.getData(k);
          // If k vertex is not fixed, process the edge between z and k.
          if (!k_data.fixed.at(index)){
            auto z_k_dist = graph.getEdgeData(e).dist;
            k_data.pred.at(index)--;
            if (k_data.pred.at(index) <= 0) {
	      k_data.fixed.at(index) = true;
              r_set.push_back(k);
            }
            if (k_data.dist.at(index) > z_data.dist.at(index) + z_k_dist) {
              k_data.dist.at(index) = z_data.dist.at(index) + z_k_dist;
              if (!k_data.fixed.at(index)) {
                pushWrap(heap, k_data.node, k_data.dist.at(index));
              }
            }
          }
        }

        if (r_set.empty()) break;
        z = r_set.back();
        r_set.pop_back();
      }
    }
  }
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
                const P& pushWrap, const R& edgeRange, uint32_t index=0) {
  using Heap = galois::MinHeap<T>;
  using Dist = typename Traits::Dist;
  uint32_t heap_operations = 0;
  uint32_t heap_pushes = 0;

  Heap heap;
  graph.getData(source).dist.at(index) = 0;
  pushWrap(heap, source, 0);
  heap_operations++;
  heap_pushes++;
  auto d = 0;

  galois::gstl::Deque<GNode> r_set;
  galois::gstl::Vector<GNodeData*> q_set;

  while (!heap.empty()) {
    T j = heap.pop();
    heap_operations++;

    auto& j_data = graph.getData(j.src);

    if (j_data.dist.at(index) < j.dist) {
      continue;
    }

    if (!j_data.fixed.at(index)) {
      j_data.fixed.at(index) = true;
      GNode z = j.src;
      d = j_data.dist.at(index);

      while(true) {
        auto& z_data = graph.getData(z);
        for (auto e : edgeRange(z)) {
          GNode k = graph.getEdgeDst(e);
          auto& k_data = graph.getData(k);
          if (!k_data.fixed.at(index)) {

            Dist z_k_dist = graph.getEdgeData(e).dist;
	    k_data.pred.at(index)--;
	    auto k_dist = z_data.dist.at(index) + z_k_dist;

            if ((k_data.pred.at(index) <=0) || k_dist <= (d + k_data.min_in_weight)){
                k_data.fixed.at(index) = true;
                r_set.push_back(k);
            }

            if (k_data.dist.at(index) > z_data.dist.at(index) + z_k_dist){
              k_data.dist.at(index) = z_data.dist.at(index) + z_k_dist;

	      if (!k_data.fixed.at(index)) {
		q_set.push_back(&k_data);
              }
	    }
          }
        }

        if (r_set.empty()) break;
        z = r_set.front();
        r_set.pop_front();
      }

      for (auto& z : q_set) {
        auto& z_data = *z;
        pushWrap(heap, z_data.node, z_data.dist.at(index));
        heap_operations++;
        heap_pushes++;
      }
      q_set.clear();
    }
  }
  heap_operations_global += heap_operations;
  heap_pushes_global += heap_pushes;
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
void init_graph(Graph& graph, GNode& source, GNode& report, bool reset_dist = true) {

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
                 [&](GNode n) {
                   auto& data = graph.getData(n);
                   data.reset(reset_dist);
                   data.node = n;
                 });

  graph.getData(source).dist = 0;

  profiling_start = clock_type::now();
  explored_nodes = 0;
  fixed_nodes = 0;
}

template <typename GType,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void init_graphAPSP(Graph& graph) {

  using SSSP = typename Traits::SSSP;

  galois::do_all(galois::iterate(graph),
                 [&graph](GNode n) {
                   auto& data = graph.getData(n);
                   data.min_in_weight = SSSP::DIST_INFINITY;
                   data.node = n;
		   if(APSP){
                     data.dist.resize(graph.size());
                     data.pred.resize(graph.size());
                     data.fixed.resize(graph.size());
		     for(GNode index : graph){
                       data.dist.at(index) = SSSP::DIST_INFINITY;
		       data.pred.at(index) = 0;
		       data.fixed.at(index) = false;
		     }
		   }
		   else{
                     data.dist.resize(1);
                     data.pred.resize(1);
                     data.fixed.resize(1);
                     data.dist.at(0) = SSSP::DIST_INFINITY;
                     data.pred.at(0) = 0;
                     data.fixed.at(0) = false;
		   }
                 });
}

template <typename GType,
	  typename T,
	  typename P,
          typename Graph = typename GType::Graph,
          typename GNode = typename GType::GNode,
          typename Traits = typename GType::Traits>
void all_pairs_shortest_path(Graph& graph_ref, uint32_t algo, std::string file, T& edgeRange, P& pushWrap) {
  using SSSP = typename Traits::SSSP;

  heap_operations_global = 0;
  heap_pushes_global = 0;

  init_graphAPSP<GType>(graph_ref);
  galois::gstl::Vector<GNode> ctx;
  calc_graph_predecessorsAPSP<GType>(graph_ref, edgeRange);

  for(GNode vertex : graph_ref){
    ctx.push_back(vertex);
  }

  auto parallel_loop = [&](GNode source, auto& ctx){
    switch(algo){
      case dijkstra : {
	  dijkstraAlgo<GType, GType::Traits::UpdateRequest>(graph_ref, source, pushWrap, edgeRange, source);
          break;
        }
      case serDelta : {
          serDeltaAlgo<GType, GType::Traits::UpdateRequest>(graph_ref, source, pushWrap, edgeRange, source);
          break;
        }
      case serSP1 : {
          serSP1Algo<GType, GType::Traits::UpdateRequest>(graph_ref, source, pushWrap, edgeRange, source);
          break;
        }
      case serSP2 : {
          serSP2Algo<GType, GType::Traits::UpdateRequest>(graph_ref, source, pushWrap, edgeRange, source);
          break;
        }
      default : break;
    }
  };

  galois::for_each(galois::iterate(ctx),
                  parallel_loop,
                  galois::wl<galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>>(),
                  galois::no_conflicts(), 
                  galois::loopname("All Pairs Shortest Path run 1"));

  std::cout<< "Total heap operations for run 1 " << heap_operations_global << std::endl;
  std::cout<< "Total heap pushes for run 1 " << heap_pushes_global << std::endl;
 

  heap_operations_global = 0;
  heap_pushes_global = 0;

  galois::do_all(galois::iterate(graph_ref),
                 [&graph_ref](GNode n) {
                   auto& data = graph_ref.getData(n);
                   if(APSP){
                     for(GNode index : graph_ref){
                       data.dist.at(index) = SSSP::DIST_INFINITY;
                       data.pred.at(index) = 0;
                       data.fixed.at(index) = false;
                     }
                   }
                 });

  calc_graph_predecessorsAPSP<GType>(graph_ref, edgeRange);

  galois::for_each(galois::iterate(ctx),
                  parallel_loop,
                  galois::wl<galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>>(),
                  galois::no_conflicts(), 
                  galois::loopname("All Pairs Shortest Path run 2"));

  std::cout<< "Total heap operations for run 2 " << heap_operations_global << std::endl;
  std::cout<< "Total heap pushes for run 2 " << heap_pushes_global << std::endl;
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
void verify_and_reportAPSP(Graph& graph, GNode& source, GNode& report, std::string algo_name, uint32_t index = 0) {
  galois::reportPageAlloc("MeminfoPost");

  std::cout << "Node " << reportNode << " has distance " << graph.getData(report).dist.at(index) << "\n";
  bool verify = true;

  if(algo == serDelta || algo == dijkstra || algo == serSP1 || algo == serSP2){
    if (!skipVerify && !APSP) {
      using SerGType = GraphType<GraphTypes::SerGraphType>;
      SerGType::Graph graph_verify;
      galois::graphs::readGraph(graph_verify, filename);
      for(GNode vertex:graph_verify){ 
        graph_verify.getData(vertex).dist = graph.getData(vertex).dist.at(index);
      }
      if (SerGType::Traits::SSSP::verify(graph_verify, source)) {
        std::cout << "Verification successful.\n";
      } else {
        GALOIS_DIE("Verification failed");
      }
    }
    else if(!skipVerify && APSP){
      srand(time(NULL));
      using SerGType = GraphType<GraphTypes::SerGraphType>;
      SerGType::Graph graph_verify;
      galois::graphs::readGraph(graph_verify, filename);
      for(int i=0; i<VERIFY ; i++){
        index = rand()%graph_verify.size();

        for(GNode vertex:graph_verify)
          graph_verify.getData(vertex).dist = graph.getData(vertex).dist.at(index);
	
        if(SerGType::Traits::SSSP::verify(graph_verify, index)){
          verify &= true;
	  std::cout << "Verified graph successfully with source as " << index << std::endl;
        }
        else{ 
          verify &= false;
	  std::cout << "Verified graph failed with source as " << index << std::endl;
	}
      }
      if(verify == true)
	std::cout << "Verification successfull" << std::endl;
      else
	std::cout << "Verification failed, some graphs in All Pairs Shortest Path have incorrect distances" << std::endl;
    }
  }
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
  profile << "num_explores;fixed_usec;fixed_idx;fixed_source;fixed_source_size;explored_source;explored_source_size;" << std::endl;

  for (auto& node : graph) {
    auto& d = graph.getData(node);
    if (d.dist == 2147483646) continue;
    profile << d.node << ";" << d.pred << ";" << d.dist << ";" << d.fixed << ";" << d.min_in_weight;
    profile << d.num_visits << ";" << d.first_visit_usec << ";" << d.first_visit_idx << ";";
    profile << d.last_visit_usec << ";" << d.visits_after_fixed << ";" << d.visits_after_explored << ";";
    profile << d.first_explored_usec << ";" << d.first_explored_idx << ";";
    profile << d.num_explores << ";" << d.fixed_usec << ";" << d.fixed_idx << ";";
    profile << d.ctx_fixed.vertex_source << ";" << d.ctx_fixed.set_size << ";";
    profile << d.ctx_explored.vertex_source << ";" << d.ctx_explored.set_size << std::endl;
  }

  profile.close();
}


int main(int argc, char** argv) {
  using ParGType = GraphType<GraphTypes::ParGraphType>;
  using SerGType = GraphType<GraphTypes::SerGraphType>;
  using APSPGType = GraphType<GraphTypes::APSPGraphType>;

  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;
  galois::StatTimer Tmain;

  switch (algo) {
    case deltaTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
        init_graph<ParGType>(graph, source, report);

        Tmain.start();
        deltaStepAlgo<ParGType, ParGType::Traits::SrcEdgeTile>(
            graph, source,
            ParGType::Traits::SrcEdgeTilePushWrap{graph},
            ParGType::Traits::TileRangeFn());
        Tmain.stop();

        verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case deltaStep: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
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
      }
      break;
    }
    case serDeltaTile: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      Graph graph;
      GType::GNode source, report;
      auto pushWrap = GType::Traits::SrcEdgeTilePushWrap{graph};
      auto edgeRange =  GType::Traits::TileRangeFn();

      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);
        auto it = graph.begin();
        std::advance(it, startNode);
        source = *it;
        it = graph.begin();
        std::advance(it, reportNode);
        report = *it;

        Tmain.start();
        serDeltaAlgo<GType, GType::Traits::SrcEdgeTile>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();

        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case serDelta: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      Graph graph;
      GType::GNode source, report;
      auto edgeRange = GType::Traits::OutEdgeRangeFn{graph};
      auto pushWrap = GType::Traits::ReqPushWrap();
      auto it = graph.begin();
      std::advance(it, startNode);
      source = *it;
      it = graph.begin();
      std::advance(it, reportNode);
      report = *it;

      if(APSP){
        galois::graphs::readGraph(graph, filename);
        all_pairs_shortest_path<GType>(graph, algo, filename, edgeRange, pushWrap);
        if(VERIFY > graph.size() && APSP){
          std::cout << "Number of nodes to verify is greater than Graph size" << std::endl;
          std::cout << "Graph size: "<< graph.size() << std::endl;
          std::cout << "Number of nodes to verify specified:  "<< VERIFY << std::endl;
          GALOIS_DIE();
        }
        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);
        Tmain.start();
        serDeltaAlgo<GType, GType::Traits::UpdateRequest>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();
	
	verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case dijkstraTile: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      GType::Graph graph;
      GType::GNode source, report;
      auto pushWrap = GType::Traits::SrcEdgeTilePushWrap{graph};
      auto edgeRange = GType::Traits::TileRangeFn();

      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);
        auto it = graph.begin();
        std::advance(it, startNode);
        source = *it;
        it = graph.begin();
        std::advance(it, reportNode);
        report = *it;

        Tmain.start();
        dijkstraAlgo<GType, GType::Traits::SrcEdgeTile>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();

        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case dijkstra: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      Graph graph;
      GType::GNode source, report;
      auto edgeRange = GType::Traits::OutEdgeRangeFn{graph};
      auto pushWrap = GType::Traits::ReqPushWrap();
      auto it = graph.begin();
      std::advance(it, startNode);
      source = *it;
      it = graph.begin();
      std::advance(it, reportNode);
      report = *it;

      if(APSP){
        galois::graphs::readGraph(graph, filename);
        all_pairs_shortest_path<GType>(graph, algo, filename, edgeRange, pushWrap);
        if(VERIFY > graph.size() && APSP){
          std::cout << "Number of nodes to verify is greater than Graph size" << std::endl;
          std::cout << "Graph size: "<< graph.size() << std::endl;
          std::cout << "Number of nodes to verify specified:  "<< VERIFY << std::endl;
          GALOIS_DIE();
        }
        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);

        Tmain.start();
        dijkstraAlgo<GType, GType::Traits::UpdateRequest>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();
        std::cout<< "Total heap operations " << heap_operations_global << std::endl;
        std::cout<< "Total heap pushes " << heap_pushes_global << std::endl;

        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case topo: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
        init_graph<ParGType>(graph, source, report);

        Tmain.start();
        topoAlgo<ParGType>(graph, source);
        Tmain.stop();

        verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case topoTile: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
        init_graph<ParGType>(graph, source, report);

        Tmain.start();
        topoTileAlgo<ParGType>(graph, source);
        Tmain.stop();

        verify_and_report<ParGType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case serSP1: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      Graph graph;
      GType::GNode source, report;
      auto edgeRange = GType::Traits::OutEdgeRangeFn{graph};
      auto pushWrap = GType::Traits::ReqPushWrap();
      if(APSP){
        galois::graphs::readGraph(graph, filename);
        all_pairs_shortest_path<GType>(graph, algo, filename, edgeRange, pushWrap);
        if(VERIFY > graph.size() && APSP){
          std::cout << "Number of nodes to verify is greater than Graph size" << std::endl;
          std::cout << "Graph size: "<< graph.size() << std::endl;
          std::cout << "Number of nodes to verify specified:  "<< VERIFY << std::endl;
          GALOIS_DIE();
        }
        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);
        auto it = graph.begin();
        std::advance(it, startNode);
        source = *it;
        it = graph.begin();
        std::advance(it, reportNode);
        report = *it;
	
        calc_graph_predecessorsAPSP<GType>(graph, edgeRange);

        Tmain.start();
        serSP1Algo<GType, GType::Traits::UpdateRequest>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();

        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case serSP2: {
      using GType = APSPGType;
      using Graph = GType::Graph;
      Graph graph;
      GType::GNode source, report;
      auto edgeRange = GType::Traits::OutEdgeRangeFn{graph};
      auto pushWrap = GType::Traits::ReqPushWrap();
      auto it = graph.begin();
      std::advance(it, startNode);
      source = *it;
      it = graph.begin();
      std::advance(it, reportNode);
      report = *it;
      if(APSP){
        galois::graphs::readGraph(graph, filename);
        all_pairs_shortest_path<GType>(graph, algo, filename, edgeRange, pushWrap);
        if(VERIFY > graph.size() && APSP){
          std::cout << "Number of nodes to verify is greater than Graph size" << std::endl;
          std::cout << "Graph size: "<< graph.size() << std::endl;
          std::cout << "Number of nodes to verify specified:  "<< VERIFY << std::endl;
          GALOIS_DIE();
        }
        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      else{
        galois::graphs::readGraph(graph, filename);
        init_graphAPSP<GType>(graph);
        calc_graph_predecessorsAPSP<GType>(graph, edgeRange);

        Tmain.start();
        serSP2Algo<GType, GType::Traits::UpdateRequest>(
            graph, source,
            pushWrap,
            edgeRange);
        Tmain.stop();
        std::cout<< "Total heap operations " << heap_operations_global << std::endl;
        std::cout<< "Total heap pushes " << heap_pushes_global << std::endl;
        verify_and_reportAPSP<GType>(graph, source, report, ALGO_NAMES[algo]);
      }
      break;
    }
    case parSP1: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
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
      }
      break;
    }
    case parSP2V: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
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
      }
      break;
    }
    case parSP2E: {
      ParGType::Graph graph;
      ParGType::GNode source, report;
      if(APSP){
        std::cout << "All Pairs Shortest Path implementation is not supported by " << ALGO_NAMES[algo] << std::endl;
	std::abort();
      }
      else{
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
      }
      break;
    }
    default:
      std::abort();
  }

  return 0;
}
