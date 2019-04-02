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

#include "stack_vector.h"
#include "bf/all.hpp"

#include <fstream>
#include <iostream>
#include <unordered_set>
#include <cstdlib>

#ifndef NDEBUG
#define D(x) x
#else
#define D(x)
#endif

#ifndef STACK_BITVECTOR_SIZE
 #define STACK_BITVECTOR_SIZE 16384
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
static cll::opt<bool>
    apsp("apsp",
              cll::desc("Run All Pairs Shortest Path version of the algorithm, if supported"),
              cll::init(false));
static cll::opt<bool>
    prof("prof",
         cll::desc("Run a profiling run of the chosen algorithm"),
         cll::init(false));

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

const char* const ALGO_NAMES[] = {"deltaTile", "deltaStep", "serDeltaTile",
                                  "serDelta",  "dijkstraTile", "dijkstra",
                                  "topo", "topoTile", "serSP1", "serSP2"};

static cll::opt<Algo>
    algo("algo", cll::desc("Choose an algorithm:"),
         cll::values(clEnumVal(deltaTile, "deltaTile"),
                     clEnumVal(deltaStep, "deltaStep"),
                     clEnumVal(serDeltaTile, "serDeltaTile"),
                     clEnumVal(serDelta, "serDelta"),
                     clEnumVal(dijkstraTile, "dijkstraTile"),
                     clEnumVal(dijkstra, "dijkstra"),
                     clEnumVal(topo, "topo"),
                     clEnumVal(topoTile, "topoTile"),
                     clEnumVal(serSP1, "serSP1"),
                     clEnumVal(serSP2, "serSP2"),
                     clEnumValEnd),
         cll::init(deltaTile));

constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

typedef std::chrono::steady_clock clock_type;
typedef std::chrono::time_point<clock_type> timestamp;
typedef uint32_t GNode;

static timestamp profiling_start;

static std::atomic<uint64_t> fixed_nodes;
static std::atomic<uint64_t> visited_nodes;
static std::atomic<uint64_t> explored_nodes;

constexpr static const uint32_t DIST_INF = std::numeric_limits<uint32_t>::max() / 2 - 1;

enum VertexSource {
  R_SET,
  HEAP
};

struct VertexContext {
  VertexSource vertex_source;
  size_t set_size;
};

struct NodeConstants {
  GNode node;
  uint32_t total_pred;
  uint32_t min_in_weight;
  uint32_t min_out_weight;

  NodeConstants() : node(-1), total_pred(0), min_in_weight(DIST_INF), min_out_weight(DIST_INF) {}
};

struct SerSourceData {
  typedef uint32_t Dist;
  NodeConstants* node_constants;
  int pred;
  Dist dist;
  bool fixed;

  SerSourceData()
      : pred(0),
        dist(DIST_INF),
        fixed(false) {}

  inline void reset(bool reset_dist) {
    if (reset_dist) dist = DIST_INF;
    pred = node_constants->total_pred;
    fixed = false;
  }
};


struct ParSourceData {
  typedef std::atomic<uint32_t> Dist;
  NodeConstants* node_constants;
  std::atomic<int> pred;
  Dist dist;
  std::atomic<bool> fixed;

  ParSourceData()
      : pred(0),
        dist(DIST_INF),
        fixed(false) {}

  inline void reset(bool reset_dist) {
    if (reset_dist) dist = DIST_INF;
    pred = node_constants->total_pred;
    fixed = false;
  }
};

template<typename SourceData,
         bool single_source = true>
struct NodeData {
  typedef SourceData SourceDataType;
  int size;
  NodeConstants node_constants;

  inline SourceData& source_data_at(int index);
  inline void reset(bool reset_dist = true);

  inline uint32_t dist_at(int index = 0) {
    return source_data_at(index).dist;
  }

  inline int pred(int index = 0) {
    return source_data_at(index).pred;
  }

  inline bool fixed(int index = 0) {
    return source_data_at(index).fixed;
  }

  inline void fix(VertexContext context = VertexContext(), int index = 0) {
    source_data_at(index).fixed = true;
  }

  inline void visit(int index = 0) {}
  inline void explore(VertexContext context = VertexContext(), int index = 0) {}
};

template<typename SourceData>
struct NodeData<SourceData, true> {
  typedef SourceData SourceDataType;
  int size;
  NodeConstants node_constants;
  SourceData source_data;
  inline SourceData& source_data_at(int index);
  inline void reset(bool reset_dist = true);
  NodeData(int size_ = 1);
};

template<typename SourceData>
SourceData& NodeData<SourceData, true>::source_data_at(int index) {
  return source_data;
}

template<typename SourceData>
inline void NodeData<SourceData, true>::reset(bool reset_dist) {
  source_data.reset(reset_dist);
}

template<typename SourceData>
NodeData<SourceData, true>::NodeData(int size_)
    : size(size_) {
  source_data.node_constants = &node_constants;
  reset();
}

template<typename SourceData>
struct NodeData<SourceData, false> {
  typedef SourceData SourceDataType;
  int size;
  NodeConstants node_constants;
  std::unique_ptr<SourceData[]> source_data;
  inline SourceData& source_data_at(int index);
  inline void reset(bool reset_dist = true);
};

template<typename SourceData>
SourceData& NodeData<SourceData, false>::source_data_at(int index) {
  return source_data[index];
}

template<typename SourceData>
inline void NodeData<SourceData, false>::reset(bool reset_dist) {
  if (!source_data) {
    source_data.reset(new SourceData[size]);
  }
  for (int i = 0; i < size; i++) {
    source_data_at(i).node_constants = &node_constants;
    source_data_at(i).reset(reset_dist);
  }
}

struct ProfilingSourceData {
  typedef std::atomic<uint32_t> Dist;
  NodeConstants* node_constants;
  std::atomic<int> pred;
  std::atomic<uint32_t> dist;
  std::atomic<bool> fixed;

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

  void reset(bool reset_dist = true) {
    pred = node_constants->total_pred;
    if (reset_dist) dist = DIST_INF;
    fixed = false;
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

template<>
struct NodeData<ProfilingSourceData, true> {
  int size;
  NodeConstants node_constants;
  std::unique_ptr<ProfilingSourceData[]> source_data;
  typedef ProfilingSourceData SourceDataType;

  NodeData(int size_ = 1) : size(size_),
                            source_data(new ProfilingSourceData[size]) {
    reset();
  }

  ProfilingSourceData& source_data_at(int index) {
    return source_data[index];
  }

  inline void explore(VertexContext context = VertexContext(), int index = 0) {
    auto& source_data = source_data_at(index);
    if (source_data.first_explored_usec == 0) {
      source_data.num_explores++;
      source_data.ctx_explored = context;
      timestamp now = clock_type::now();
      source_data.first_explored_usec = (now - profiling_start).count();
      source_data.first_explored_idx = explored_nodes++;
    }
  }

  inline void fix(VertexContext context = VertexContext(), int index = 0) {
    auto& source_data = source_data_at(index);
    if (source_data.fixed) return;
    timestamp now = clock_type::now();
    source_data.fixed = true;
    source_data.fixed_usec = (now - profiling_start).count();
    source_data.fixed_idx = fixed_nodes++;
    source_data.ctx_fixed = context;
  }

  inline void visit(int index = 0) {
    auto& source_data = source_data_at(index);
    if (source_data.num_visits == 0) {
      timestamp now = clock_type::now();
      source_data.first_visit_usec = (now - profiling_start).count();
      source_data.first_visit_idx = visited_nodes++;
    }
    source_data.num_visits++;
    if (source_data.num_explores > 0) {
      source_data.visits_after_explored++;
    }
    if (source_data.fixed) {
      timestamp now = clock_type::now();
      source_data.visits_after_fixed++;
      source_data.last_visit_usec = (now - profiling_start).count();
    }
  }

  inline void reset(bool reset_dist = true) {
    for (int i = 0; i < size; i++) {
      source_data[i].reset(reset_dist);
    }
  }
};

typedef typename galois::graphs::LC_CSR_Graph<NodeData<SerSourceData, true>,
                                              uint32_t>
::with_no_lockable<true>::type
::with_numa_alloc<true>::type SerialSSGraph;

typedef typename galois::graphs::LC_CSR_Graph<NodeData<ParSourceData, true>,
                                              uint32_t>
::with_no_lockable<true>
::type::with_numa_alloc<true>::type ParallelSSGraph;

typedef typename galois::graphs::LC_CSR_Graph<NodeData<SerSourceData, false>,
                                              uint32_t>
::with_no_lockable<true>::type
::with_numa_alloc<true>::type SerialASGraph;

typedef typename galois::graphs::LC_CSR_Graph<NodeData<ParSourceData, false>,
                                              uint32_t>
::with_no_lockable<true>
::type::with_numa_alloc<true>::type ParallelASGraph;

typedef typename galois::graphs::LC_CSR_Graph<NodeData<ProfilingSourceData, true>,
                                              uint32_t>
::with_no_lockable<true>
::type::with_numa_alloc<true>::type ProfilingGraph;

template<typename Graph>
struct GraphTraits {
  typedef BFS_SSSP<Graph, uint32_t, true, EDGE_TILE_SIZE> SSSP;
  typedef typename SSSP::Dist Dist;
  typedef typename SSSP::UpdateRequest UpdateRequest;
  typedef typename SSSP::UpdateRequestIndexer UpdateRequestIndexer;
  typedef typename SSSP::SrcEdgeTile SrcEdgeTile;
  typedef typename SSSP::SrcEdgeTileMaker SrcEdgeTileMaker;
  typedef typename SSSP::SrcEdgeTilePushWrap SrcEdgeTilePushWrap;
  typedef typename SSSP::ReqPushWrap ReqPushWrap;
  typedef typename SSSP::OutEdgeRangeFn OutEdgeRangeFn;
  typedef typename SSSP::TileRangeFn TileRangeFn;
  typedef typename Graph::node_data_type NodeData;
  typedef typename NodeData::SourceDataType SourceData;
};


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


template <typename Graph>
void init_graph(Graph& graph, GNode& source, GNode& report) {

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
  galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());
  galois::reportPageAlloc("MeminfoPre");
}

template<typename Graph>
void calc_graph_predecessors(Graph& graph) {
  // Fill the pred array
  for (auto vertex : graph) {
    auto& v_data = graph.getData(vertex);
    for (auto edge : graph.edges(vertex, galois::MethodFlag::UNPROTECTED)) {
      auto edge_dist = graph.getEdgeData(edge);
      auto& k_data = graph.getData(graph.getEdgeDst(edge));
      k_data.node_constants.total_pred++;
      k_data.node_constants.min_in_weight = std::min(k_data.node_constants.min_in_weight, edge_dist);
      if (v_data.node_constants.min_out_weight > edge_dist) {
        v_data.node_constants.min_out_weight = edge_dist;
      }
    }
  }
}

template<typename Graph>
void dump_profile(Graph& graph, std::ofstream& profile) {}

template<>
void dump_profile<ProfilingGraph>(ProfilingGraph& graph, std::ofstream& profile) {
  for (auto& node : graph) {
    // TODO support profiling beyond the first source
    auto& d =  graph.getData(node).source_data_at(0);
    if (d.dist == 2147483646) continue;
    profile << d.node_constants->node << ";" << d.pred << ";" << d.dist << ";" << d.fixed << ";" << d.node_constants->min_in_weight;
    profile << d.num_visits << ";" << d.first_visit_usec << ";" << d.first_visit_idx << ";";
    profile << d.last_visit_usec << ";" << d.visits_after_fixed << ";" << d.visits_after_explored << ";";
    profile << d.first_explored_usec << ";" << d.first_explored_idx << ";";
    profile << d.num_explores << ";" << d.fixed_usec << ";" << d.fixed_idx << ";";
    profile << d.ctx_fixed.vertex_source << ";" << d.ctx_fixed.set_size << ";";
    profile << d.ctx_explored.vertex_source << ";" << d.ctx_explored.set_size << std::endl;
  }
}

template<typename T>
class AlgoRunner {
 public:
  virtual void init() = 0;
  virtual std::string name() = 0;
  virtual uint64_t node_dist(GNode node, uint32_t index = 0) = 0;
  virtual void reset(bool reset_distances = true) = 0;
  // TODO make children return an iterable.
  virtual void write_profile(std::ofstream& profile, uint32_t index) = 0;


  virtual bool verify(uint32_t index = 0) {
    galois::reportPageAlloc("MeminfoPost");

    std::cout << "Node " << reportNode << " has distance "
              << node_dist(reportNode) << "\n";

    if (skipVerify) {
      std::cout << "Skipping verification." << std::endl;
      return true;
    }

    if (do_verify()) {
      std::cout << "Verification successful.\n";
    } else {
      GALOIS_DIE("Verification failed");
    }

    if (!prof) return true;

    std::string graph_name = removeExtension(basename(filename));
    std::string prof_filename = "per_node_profile_" + name() + "_" + graph_name + ".csv";
    std::cout << "Outputting profile to file: " << prof_filename << std::endl;

    // For profiling graphs, output the metrics in csv form
    std::ofstream profile;
    profile.open(prof_filename);
    profile << "node_idx;pred;dist;fixed;min_in_weight;num_visits;first_visit_usec;first_visit_idx";
    profile << "last_visit_usec;visits_after_fixed;visits_after_explored;first_explored_usec;first_explored_idx;";
    profile << "num_explores;fixed_usec;fixed_idx;fixed_source;fixed_source_size;explored_source;explored_source_size;" << std::endl;
    write_profile(profile, index);

    profile.close();

    return true;
  }

  void run(std::string run_name) {
    std::cout << "Running " << name() << " algorithm for run: " << run_name << std::endl;
    galois::StatTimer tmain;
    tmain.start();
    profiling_start = clock_type::now();
    explored_nodes = 0;
    fixed_nodes = 0;
    if (!apsp) {
      do_run(source);
    } else {
      do_run_apsp();
    }
    tmain.stop();
    galois::runtime::reportStat_Single("AlgoRunner[" + run_name + "]", "Total time", tmain.get());
  }


 protected:
  virtual bool do_verify() = 0;
  virtual void do_run(GNode source) = 0;
  virtual void do_run_apsp() = 0;

  GNode source;
  GNode report;
};

template<typename Graph,
         typename T = typename GraphTraits<Graph>::SSSP::UpdateRequest>
class GraphAlgoBase : public AlgoRunner<typename GraphTraits<Graph>::SSSP::UpdateRequest> {
 public:
  virtual uint64_t node_dist(GNode node, uint32_t index = 0) {
    return graph.getData(node).source_data_at(index).dist;
  }

  virtual void init() {
    init_graph(graph, this->source, this->report);
    calc_graph_predecessors(graph);
    reset(true);
  }

  virtual void reset(bool reset_distances = true) {
    galois::do_all(galois::iterate(graph),
                   [&](GNode n) {
                     auto& data = graph.getData(n);
                     data.node_constants.node = n;
                     if (!apsp && data.size != 1) GALOIS_DIE("Not APSP, only need one source data per node");
                     if (apsp) data.size = graph.size();
                     data.reset(reset_distances);
                   });
  }

  virtual void do_run_apsp() {
    auto parallel_loop = [this](GNode source, auto& ctx) {
                           this->do_run(source);
                         };
    galois::for_each(galois::iterate(graph),
                     parallel_loop,
                     galois::wl<galois::worklists::PerSocketChunkFIFO<CHUNK_SIZE>>(),
                     galois::no_conflicts(),
                     galois::loopname("All Pairs Shortest Path"));
  }

  virtual bool do_verify() {
    if (!apsp) {
      return GraphTraits<Graph>::SSSP::verify(graph, this->source);
    } else {
      GraphTraits<Graph>::SSSP::verify(graph, this->source);
      int random_node  = rand() % graph.size();
      std::cout << "Randomly veryfying node: " << random_node << std::endl;
      return GraphTraits<Graph>::SSSP::verify(graph, random_node);
    }
  }

  virtual void write_profile(std::ofstream& profile, uint32_t index) {
    dump_profile(graph, profile);
  }

  Graph graph;
};

template <typename Graph,
          typename T,
          typename P,
          typename R,
          typename Traits = GraphTraits<Graph>>
void deltaStepAlgo(Graph& graph,
                   GNode source,
                   const P& pushWrap,
                   const R& edgeRange,
                   int index = 0) {

  using Dist                 = typename Traits::Dist;
  using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> BadWork;
  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> WLEmptyWork;

  namespace gwl = galois::worklists;

  using PSchunk = gwl::PerSocketChunkFIFO<CHUNK_SIZE>;
  using OBIM = galois::worklists::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

  graph.getData(source).source_data_at(index).dist = 0;

  galois::InsertBag<T> initBag;
  pushWrap(initBag, source, 0, "parallel");

  galois::for_each(galois::iterate(initBag),
                   [&](const T& item, auto& ctx) {
                     constexpr galois::MethodFlag flag =
                         galois::MethodFlag::UNPROTECTED;
                     const auto& sdata = graph.getData(item.src, flag).source_data_at(index);

                     if (sdata.dist < item.dist) {
                       if (TRACK_WORK)
                         WLEmptyWork += 1;
                       return;
                     }

                     for (auto ii : edgeRange(item)) {

                       GNode dst          = graph.getEdgeDst(ii);
                       auto& ddata        = graph.getData(dst, flag).source_data_at(index);
                       Dist ew            = graph.getEdgeData(ii, flag);
                       const Dist newDist = sdata.dist + ew;


                       while (true) {
                         Dist oldDist = ddata.dist;

                         if (oldDist <= newDist) {
                           break;
                         }

                         if (ddata.dist.compare_exchange_weak(
                                 oldDist, newDist, std::memory_order_relaxed)) {

                           if (TRACK_WORK) {
                             //! [per-thread contribution of self-defined stats]
                             if (oldDist != Traits::SSSP::DIST_INFINITY) {
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

template<typename Graph,
         typename T = typename GraphTraits<Graph>::SSSP::UpdateRequest>
class DeltaStepAlgoRunner : public GraphAlgoBase<Graph, T> {
 public:
  using PushWrap = typename GraphTraits<Graph>::ReqPushWrap;
  using EdgeRange = typename GraphTraits<Graph>::OutEdgeRangeFn;

  DeltaStepAlgoRunner() : graph(GraphAlgoBase<Graph, T>::graph), er({graph}) {}

  virtual std::string name() {
    std::string ret = "Deltastep";
    return ret + (apsp ? "[APSP]" : "[SSSP]") + (prof ? "-prof" : "");
  }

  virtual void do_run(GNode source) {
    deltaStepAlgo<Graph, T>(graph, source, pw, er);
  }

  Graph& graph;
  PushWrap pw;
  EdgeRange er;
};

template <typename Graph,
          typename T,
          typename P,
          typename R,
          typename Traits = GraphTraits<Graph>>
void serDeltaAlgo(Graph& graph,
                  const GNode& source,
                  const P& pushWrap,
                  const R& edgeRange,
                  uint32_t index = 0) {

  using UpdateRequestIndexer = typename Traits::UpdateRequestIndexer;

  SerialBucketWL<T, UpdateRequestIndexer> wl(UpdateRequestIndexer{stepShift});
  graph.getData(source)->source_data_at(index).dist = 0;

  pushWrap(wl, source, 0);

  while (!wl.empty()) {

    auto& curr = wl.minBucket();

    while (!curr.empty()) {
      auto item = curr.front();
      curr.pop_front();

      auto& item_data = graph.getData(item.src).source_data_at(index);
      if (item_data.dist < item.dist) {
       continue;
      }

      for (auto e : edgeRange(item)) {
        GNode dst   = graph.getEdgeDst(e);
        auto& ddata = graph.getData(dst).source_data_at(index);

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
}


template <typename Graph,
          typename T,
          typename P,
          typename R>
void dijkstraAlgo(Graph& graph,
                  const GNode& source,
                  const P& pushWrap,
                  const R& edgeRange) {

  using WL = galois::MinHeap<T>;
  graph.getData(source).source_data_at(source).dist = 0;

  WL wl;
  pushWrap(wl, source, 0);

  while (!wl.empty()) {

    T item = wl.top();
    wl.pop();

    auto& item_data = graph.getData(item.src).source_data_at(source);
    if (item_data.dist < item.dist) {
      continue;
    }

    for (auto e : edgeRange(item)) {

      GNode dst   = graph.getEdgeDst(e);
      auto& ddata = graph.getData(dst).source_data_at(source);
      auto& edata =  graph.getEdgeData(e);

      const auto newDist = item.dist + edata;

      if (newDist < ddata.dist) {
        ddata.dist = newDist;
        pushWrap(wl, dst, newDist);
      }
    }
  }
}


template<typename Graph,
         typename T = typename GraphTraits<Graph>::SSSP::UpdateRequest>
class DijkstraAlgoRunner : public GraphAlgoBase<Graph, T> {
 public:
  using PushWrap = typename GraphTraits<Graph>::ReqPushWrap;
  using EdgeRange = typename GraphTraits<Graph>::OutEdgeRangeFn;

  DijkstraAlgoRunner() : graph(GraphAlgoBase<Graph, T>::graph), er({graph}) {}

  virtual std::string name() {
    std::string ret = "Dijkstra";
    return ret + (apsp ? "[APSP]" : "[SSSP]") + (prof ? "-prof" : "");
  }

  virtual void do_run(GNode source) {
    dijkstraAlgo<Graph, T>(graph, source, pw, er);
  }

  Graph& graph;
  PushWrap pw;
  EdgeRange er;
};


template<typename Graph, typename HeapEntryType>
struct GarbageCollectFixedNodes {
  GarbageCollectFixedNodes(Graph& graph_, my_bitvector<STACK_BITVECTOR_SIZE>& fixed) : graph(graph_), fixed(fixed) {}
  Graph& graph;
  my_bitvector<STACK_BITVECTOR_SIZE>& fixed;
  inline bool operator()(HeapEntryType& item) const {
    return fixed[item.node];
  }
};

template<typename SourceData>
struct HeapNode {
  typedef typename SourceData::Dist Dist;
  GNode node;
  Dist dist;
  SourceData* sdata;

  friend bool operator<(const HeapNode& left,
                        const HeapNode& right) {
    return left.dist == right.dist ? left.node < right.node
                                                : left.dist < right.dist;
  }

  friend bool operator==(const HeapNode& left,
                         const HeapNode& right) {
    return left.node == right.node;
  }

  HeapNode& operator=(const HeapNode& other) {
    node = other.node;
    uint32_t dist_as_int = other.dist;
    dist = dist_as_int;
    sdata = other.sdata;
    return *this;
  }

  HeapNode(GNode node, uint32_t dist_, SourceData* sdata)
      : node(node), sdata(sdata){
    dist = dist_;
  }

  HeapNode(const HeapNode& hnode) {
    *this = hnode;
  }
};

template <typename Graph,
          typename T,
          typename P,
          typename R,
          typename Traits = GraphTraits<Graph>,
          typename Cmp = std::less<typename Traits::Dist>>
void serSP2Algo(Graph& graph, const GNode& source,
                const P& pushWrap, const R& edgeRange) {

  using SourceData = typename Traits::NodeData::SourceDataType;
  using HNode = HeapNode<SourceData>;

  using Heap = galois::GarbageCollectingMinHeap<HNode,
                                                GarbageCollectFixedNodes<Graph, HNode>>;

  // The set of nodes which have been fixed but not explored
  size_t nodes_fixed = 0;

  StackVector<HNode, 128> r_set;

  my_bitvector<STACK_BITVECTOR_SIZE> fixed;
  GarbageCollectFixedNodes<Graph, HNode> gc(graph, fixed);

  SourceData* sdata = &graph.getData(source).source_data_at(source);
  sdata->dist = 0;

  Heap heap(gc);
  heap.push(HNode(source, 0, sdata));

  HNode min(0, 0, nullptr);
  Cmp cmp;

  // While the heap is not empty
  while (LIKELY(nodes_fixed < graph.size()) && (!heap.empty() || !r_set->empty())) {

    if (LIKELY(!heap.empty())) {
      HNode node = heap.top();

      if (fixed[node.node] || node.sdata->dist < node.dist) {
        heap.pop();
        // If we got a min, go do some work.
        if (!r_set->empty()) goto mainloop;
        continue;
      }

      heap.pop();
      min = node;
      fixed[min.node] = true;
      r_set->push_back(min);
    }

    if (LIKELY(!heap.empty()) && heap.top().dist == min.dist) {
      continue;
    }

    mainloop:
    // Inner loop, go through the all the elements in R
    while(r_set->size() > 0) {
      HNode z = r_set->back();
      r_set->pop_back();

      // Get all the vertices that have edges from z
      for (auto e : edgeRange(z.node)) {
        GNode knode = graph.getEdgeDst(e);

        if (fixed[knode]) continue;

        auto& kndata = graph.getData(knode);
        SourceData& k = kndata.source_data_at(source);
        HNode khnode(knode, k.dist, &k);

        // If k vertex is not fixed, process the edge between z and k.
        auto z_k_dist = graph.getEdgeData(e);

        bool changed = false;
        if (cmp(z.dist + z_k_dist, k.dist)) {
          k.dist = z.dist + z_k_dist;
          uint32_t dist_as_int = k.dist;
          khnode.dist = dist_as_int;
          changed = true;
          if (k.dist < min.dist) {
            min = khnode;
          }
        }

        if (--k.pred <= 0 || cmp(k.dist, (min.dist + kndata.node_constants.min_in_weight))) {
          fixed[knode] = true;
          r_set->push_back(khnode);
        } else if (changed) {
          heap.push(HNode(knode, k.dist, &k));
        }
      }

      if (r_set->empty() && !fixed[min.node] && cmp(min.dist, heap.top().dist)) {
        // We're done, but before we break, let's just check whether we have the new min in the q set
        // That is, if the heap is not empty and the current min is higher than the min in the q
        // set no point in pushing back to the heap, where it would have to bubble up.
        fixed[min.node] = true;
        r_set->push_back(min);
      }
    }
  }
}


template<typename Graph,
         typename T = typename GraphTraits<Graph>::SSSP::UpdateRequest>
class SerSP2AlgoRunner : public GraphAlgoBase<Graph, T> {
 public:
  using PushWrap = typename GraphTraits<Graph>::ReqPushWrap;
  using EdgeRange = typename GraphTraits<Graph>::OutEdgeRangeFn;

  SerSP2AlgoRunner() : graph(GraphAlgoBase<Graph, T>::graph), er({graph}) {}

  virtual std::string name() {
    std::string ret = "SerSP2";
    return ret + (apsp ? "[APSP]" : "[SSSP]") + (prof ? "-prof" : "");
  }

  virtual void do_run(GNode source) {
    serSP2Algo<Graph, T>(graph, source, pw, er);
  }

  Graph& graph;
  PushWrap pw;
  EdgeRange er;
};

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


template<typename AlgoRunner>
void run_benchmark(AlgoRunner& runner) {

  runner.init();

  // Do a dry run without verification to warm things up.
  runner.run("init");
  // Reset everything before the main run.
  runner.reset();

  // Main bench run.
  runner.run("main");
  runner.verify();

  // Reset everything but the distances for a final baseline run.
  //runner.reset(false);
  //runner.run("baseline");
}


int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;
  galois::StatTimer Tmain;


  switch (algo) {
    case deltaStep: {
      if (!apsp && !prof) {
        DeltaStepAlgoRunner<ParallelSSGraph> runner;
        run_benchmark(runner);
      } else if (!apsp && prof) {
        DeltaStepAlgoRunner<ProfilingGraph> runner;
        run_benchmark(runner);
      } else if (apsp && !prof) {
        DeltaStepAlgoRunner<ParallelASGraph> runner;
        run_benchmark(runner);
      } else {
        std::abort();
      }
      break;
    }
    case dijkstra: {
      if (!apsp && !prof) {
        DijkstraAlgoRunner<SerialSSGraph> runner;
        run_benchmark(runner);
      } else if (!apsp && prof) {
        DijkstraAlgoRunner<ProfilingGraph> runner;
        run_benchmark(runner);
      } else if (apsp && !prof) {
        DijkstraAlgoRunner<SerialASGraph> runner;
        run_benchmark(runner);
      } else {
        std::abort();
      }
      break;
    }
    case serSP2: {
      if (!apsp && !prof) {
        SerSP2AlgoRunner<SerialSSGraph> runner;
        run_benchmark(runner);
      } else if (!apsp && prof) {
        SerSP2AlgoRunner<ProfilingGraph> runner;
        run_benchmark(runner);
      } else if (apsp && !prof) {
        SerSP2AlgoRunner<SerialASGraph> runner;
        run_benchmark(runner);
      } else {
        std::abort();
      }
      break;
    }
    default:
      std::cout << "Unsupported algorithm. Aborting."<< std::endl;
      return 1;
  }

  return 0;
}
