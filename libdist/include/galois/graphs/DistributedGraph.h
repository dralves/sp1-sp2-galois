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

/**
 * @file DistributedGraph.h
 *
 * Contains the implementation for DistGraph. Command line argument definitions
 * are found in DistributedGraph.cpp.
 */

#ifndef _GALOIS_DIST_HGRAPH_H_
#define _GALOIS_DIST_HGRAPH_H_

#include <unordered_map>
#include <fstream>

#include "galois/runtime/GlobalObj.h"
#include "galois/graphs/BufferedGraph.h"
#include "galois/graphs/B_LC_CSR_Graph.h"
#include "galois/runtime/DistStats.h"
#include "galois/graphs/OfflineGraph.h"
#include "galois/runtime/SyncStructures.h"
#include "galois/runtime/DataCommMode.h"
#include "galois/DynamicBitset.h"

#ifdef __GALOIS_HET_CUDA__
#include "galois/cuda/HostDecls.h"
#endif
#ifdef __GALOIS_HET_OPENCL__
#include "galois/opencl/CL_Header.h"
#endif

#include "galois/runtime/BareMPI.h"

#include "llvm/Support/CommandLine.h"

/*
 * Headers for boost serialization
 */
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>

namespace cll = llvm::cl;

/**
 * Enums specifying how masters are to be distributed among hosts.
 */
enum MASTERS_DISTRIBUTION {
  //! balance nodes
  BALANCED_MASTERS,
  //! balance edges
  BALANCED_EDGES_OF_MASTERS,
  //! balance nodes and edges
  BALANCED_MASTERS_AND_EDGES
};

//! Specifies if synchronization should be partition agnostic
extern cll::opt<bool> partitionAgnostic;
//! Specifies what format to send metadata in
extern cll::opt<DataCommMode> enforce_metadata;
//! Specifies how to distribute masters among hosts
extern cll::opt<MASTERS_DISTRIBUTION> masters_distribution;
//! Specifies how much weight to give to a node when
extern cll::opt<uint32_t> nodeWeightOfMaster;
//! Specifies how much weight to give to a node when
extern cll::opt<uint32_t> edgeWeightOfMaster;
//! Specifies how much weight to give to a node when
extern cll::opt<uint32_t> nodeAlphaRanges;
//! Specifies number of threads doing I/O
extern cll::opt<unsigned> numFileThreads;
//! Specifies the size of the buffer used for
extern cll::opt<unsigned> edgePartitionSendBufSize;

//! Enumeration for specifiying write location for sync calls
enum WriteLocation {
  //! write at source
  writeSource,
  //! write at destination
  writeDestination,
  //! write at source and/or destination
  writeAny
};
//! Enumeration for specifiying read location for sync calls
enum ReadLocation {
  //! read at source
  readSource,
  //! read at destination
  readDestination,
  //! read at source and/or destination
  readAny
};

namespace galois {
namespace graphs {

/**
 * Base DistGraph class that all distributed graphs extend from.
 *
 * @tparam NodeTy type of node data for the graph
 * @tparam EdgeTy type of edge data for the graph
 * @tparam WithInEdges controls whether or not it is possible to store in-edges
 * in addition to outgoing edges in this graph
 */
template <typename NodeTy, typename EdgeTy, bool WithInEdges = false>
class DistGraph : public galois::runtime::GlobalObject {
private:
  //! Graph name used for printing things
  constexpr static const char* const GRNAME = "dGraph";

  //! testing
  using GraphTy = typename std::conditional<
      WithInEdges, galois::graphs::B_LC_CSR_Graph<NodeTy, EdgeTy, true>,
      galois::graphs::LC_CSR_Graph<NodeTy, EdgeTy, true>>::type;

  //! Tracks current round number; has to be set manually by the user
  bool round;

protected:
  //! The internal graph used by DistGraph to represent the graph
  GraphTy graph;
  //! Synchronization type
  enum SyncType {
    syncReduce,   //!< Reduction sync
    syncBroadcast //!< Broadcast sync
  };

  //! Marks if the graph is transposed or not. Important when determining
  //! how to do synchronization
  bool transposed;

  // global graph variables
  uint64_t numGlobalNodes; //!< Total nodes in the global unpartitioned graph.
  uint64_t numGlobalEdges; //!< Total edges in the global unpartitioned graph.

  const unsigned id; //!< Copy of net.ID, which is the ID of the machine.
  const uint32_t
      numHosts; //!< Copy of net.Num, which is the total number of machines

  // local graph
  // size() = Number of nodes created on this host (masters + mirrors)
  uint32_t numOwned;    //!< Number of nodes owned (masters) by this host.
                        //!< size() - numOwned = mirrors on this host
  uint32_t beginMaster; //!< Local id of the beginning of master nodes.
                        //!< beginMaster + numOwned = local id of the end of
                        //!< master nodes
  uint32_t numNodesWithEdges; //!< Number of nodes (masters + mirrors) that have
                              //!< outgoing edges

  //! Information that converts GID to host that has master proxy of that node
  std::vector<std::pair<uint64_t, uint64_t>> gid2host;

  uint64_t last_nodeID_withEdges_bipartite; //!< used only for bipartite graphs

  // memoization optimization
  //! Mirror nodes from different hosts. For reduce
  std::vector<std::vector<size_t>> mirrorNodes;
  //! Master nodes on different hosts. For broadcast
  std::vector<std::vector<size_t>> masterNodes;

  //! Typedef used so galois::runtime::BITVECTOR_STATUS doesn't have to be
  //! written
  using BITVECTOR_STATUS = galois::runtime::BITVECTOR_STATUS;

  //! A pointer set during sync_on_demand calls that points to the status
  //! of a bitvector with regard to where data has been synchronized
  //! @todo pass the flag as function paramater instead
  BITVECTOR_STATUS* currentBVFlag;

private:
  // vector for determining range objects for master nodes + nodes
  // with edges (which includes masters)
  //! represents split of all nodes among threads to balance edges
  std::vector<uint32_t> allNodesRanges;
  //! represents split of master nodes among threads to balance edges
  std::vector<uint32_t> masterRanges;
  //! represents split of nodes with edges (includes masters) among threads to
  //! balance edges
  std::vector<uint32_t> withEdgeRanges;

  //! represents split of all nodes among threads to balance in-edges
  std::vector<uint32_t> allNodesRangesIn;
  //! represents split of master nodes among threads to balance in-edges
  std::vector<uint32_t> masterRangesIn;

  using NodeRangeType =
      galois::runtime::SpecificRange<boost::counting_iterator<size_t>>;

  //! Vector of ranges that stores the 3 different range objects that a user is
  //! able to access
  std::vector<NodeRangeType> specificRanges;
  //! Like specificRanges, but for in edges
  std::vector<NodeRangeType> specificRangesIn;

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
  std::vector<MPI_Group> mpi_identity_groups;
#endif

protected:
  //! Prints graph statistics.
  void printStatistics() {
    if (id == 0) {
      galois::gPrint("Total nodes: ", numGlobalNodes, "\n");
      galois::gPrint("Total edges: ", numGlobalEdges, "\n");
    }
    galois::gPrint("[", id, "] Master nodes: ", numOwned, "\n");
    galois::gPrint("[", id, "] Mirror nodes: ", size() - numOwned, "\n");
    galois::gPrint("[", id, "] Nodes with edges: ", numNodesWithEdges, "\n");
    galois::gPrint("[", id, "] Edges: ", sizeEdges(), "\n");

    // reports the number of nodes + edge as well
    galois::runtime::reportStatCond_Tsum<MORE_DIST_STATS>(
        "dGraph", "TotalNodes", numOwned);
    galois::runtime::reportStatCond_Tsum<MORE_DIST_STATS>(
        "dGraph", "TotalEdges", sizeEdges());
  }

  //! Increments evilPhase, a phase counter used by communication.
  void inline increment_evilPhase() {
    ++galois::runtime::evilPhase;
    if (galois::runtime::evilPhase >=
        std::numeric_limits<int16_t>::max()) { // limit defined by MPI or LCI
      galois::runtime::evilPhase = 1;
    }
  }

public:
  /****** VIRTUAL FUNCTIONS *********/
  //! Converts a global node id to a local node id
  //! @returns local id from passed in global id
  virtual uint32_t G2L(uint64_t) const = 0;
  //! Converts a local node id to a local node id
  //! @returns global id from passed in local id
  virtual uint64_t L2G(uint32_t) const = 0;
  //! Returns vertex cut status
  //! @returns True if the cut being used is a vertex cut
  virtual bool is_vertex_cut() const = 0;
  //! Determines which host has the master for a particular node
  //! @returns Host id of node in question
  virtual unsigned getHostID(uint64_t) const = 0;
  //! Determine if a node has a master on this host.
  //! @returns True if passed in global id has a master on this host
  virtual bool isOwned(uint64_t) const = 0;
  //! Determine if a node has a proxy on this host
  //! @returns True if passed in global id has a proxy on this host
  virtual bool isLocal(uint64_t) const = 0;
  //! Serialize this graph
  virtual void boostSerializeLocalGraph(boost::archive::binary_oarchive& ar,
                                        const unsigned int version = 0) const {}
  //! Deserialize a graph
  virtual void boostDeSerializeLocalGraph(boost::archive::binary_iarchive& ar,
                                          const unsigned int version = 0){};
  //! Gets the mirror node ranges on this graph
  //! @returns Range of mirror nodes on this graph
  virtual std::vector<std::pair<uint32_t, uint32_t>>
  getMirrorRanges() const = 0;

  // Requirement: For all X and Y,
  // On X, nothingToSend(Y) <=> On Y, nothingToRecv(X)
  // Note: templates may not be virtual, so passing types as arguments
  /**
   * Determine if we have anything that we need to send to a particular host
   *
   * @param host Host number that we may or may not send to
   * @param syncType Synchronization type to determine which nodes on a
   * host need to be considered
   * @param writeLocation If data is being written to on source or
   * destination (or both)
   * @param readLocation If data is being read from on source or
   * destination (or both)
   * @returns true if there is nothing to send to a host, false otherwise
   */
  virtual bool nothingToSend(unsigned host, SyncType syncType,
                             WriteLocation writeLocation,
                             ReadLocation readLocation) {
    auto& sharedNodes = (syncType == syncReduce) ? mirrorNodes : masterNodes;
    return (sharedNodes[host].size() == 0);
  }

  /**
   * Determine if we have anything that we need to receive from a particular
   * host
   *
   * @param host Host number that we may or may not receive from
   * @param syncType Synchronization type to determine which nodes on a
   * host need to be considered
   * @param writeLocation If data is being written to on source or
   * destination (or both)
   * @param readLocation If data is being read from on source or
   * destination (or both)
   * @returns true if there is nothing to receive from a host, false otherwise
   */
  virtual bool nothingToRecv(unsigned host, SyncType syncType,
                             WriteLocation writeLocation,
                             ReadLocation readLocation) {
    auto& sharedNodes = (syncType == syncReduce) ? masterNodes : mirrorNodes;
    return (sharedNodes[host].size() == 0);
  }

  /**
   * Reset a provided bitset given the type of synchronization performed
   *
   * @param syncType Type of synchronization to consider when doing reset
   * @param bitset_reset_range Function to reset range with
   */
  virtual void reset_bitset(SyncType syncType,
                            void (*bitset_reset_range)(size_t,
                                                       size_t)) const = 0;

private:
  uint32_t num_run;   //!< Keep track of number of runs.
  uint32_t num_round; //!< Keep track of number of rounds.

  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host. Considers
   * ONLY nodes and not edges.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param numNodes_to_divide The total number of nodes you are
   * assigning to different hosts
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  void computeMastersBlockedNodes(galois::graphs::OfflineGraph& g,
                                  uint64_t numNodes_to_divide,
                                  const std::vector<unsigned>& scalefactor,
                                  unsigned DecomposeFactor = 1) {
    if (scalefactor.empty() || (numHosts * DecomposeFactor == 1)) {
      for (unsigned i = 0; i < numHosts * DecomposeFactor; ++i)
        gid2host.push_back(galois::block_range(0U, (unsigned)numNodes_to_divide,
                                               i, numHosts * DecomposeFactor));
    } else { // TODO: not compatible with DecomposeFactor.
      assert(scalefactor.size() == numHosts);

      unsigned numBlocks = 0;

      for (unsigned i = 0; i < numHosts; ++i) {
        numBlocks += scalefactor[i];
      }

      std::vector<std::pair<uint64_t, uint64_t>> blocks;
      for (unsigned i = 0; i < numBlocks; ++i) {
        blocks.push_back(galois::block_range(0U, (unsigned)numNodes_to_divide,
                                             i, numBlocks));
      }

      std::vector<unsigned> prefixSums;
      prefixSums.push_back(0);

      for (unsigned i = 1; i < numHosts; ++i) {
        prefixSums.push_back(prefixSums[i - 1] + scalefactor[i - 1]);
      }

      for (unsigned i = 0; i < numHosts; ++i) {
        unsigned firstBlock = prefixSums[i];
        unsigned lastBlock  = prefixSums[i] + scalefactor[i] - 1;
        gid2host.push_back(
            std::make_pair(blocks[firstBlock].first, blocks[lastBlock].second));
      }
    }
  }

  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host while taking
   * into consideration the only edges of the node to get
   * even blocks.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param numNodes_to_divide The total number of nodes you are
   * assigning to different hosts
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  void computeMastersBalancedEdges(galois::graphs::OfflineGraph& g,
                                   uint64_t numNodes_to_divide,
                                   const std::vector<unsigned>& scalefactor,
                                   unsigned DecomposeFactor = 1) {
    if (edgeWeightOfMaster == 0) {
      edgeWeightOfMaster = 1;
    }

    auto& net = galois::runtime::getSystemNetworkInterface();

    gid2host.resize(numHosts * DecomposeFactor);
    for (unsigned d = 0; d < DecomposeFactor; ++d) {
      auto r = g.divideByNode(0, edgeWeightOfMaster, (id + d * numHosts),
                              numHosts * DecomposeFactor, scalefactor);
      gid2host[id + d * numHosts].first  = *(r.first.first);
      gid2host[id + d * numHosts].second = *(r.first.second);
    }

    for (unsigned h = 0; h < numHosts; ++h) {
      if (h == id)
        continue;
      galois::runtime::SendBuffer b;
      for (unsigned d = 0; d < DecomposeFactor; ++d) {
        galois::runtime::gSerialize(b, gid2host[id + d * numHosts]);
      }
      net.sendTagged(h, galois::runtime::evilPhase, b);
    }
    net.flush();
    unsigned received = 1;
    while (received < numHosts) {
      decltype(net.recieveTagged(galois::runtime::evilPhase, nullptr)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase, nullptr);
      } while (!p);
      assert(p->first != id);
      auto& b = p->second;
      for (unsigned d = 0; d < DecomposeFactor; ++d) {
        galois::runtime::gDeserialize(b, gid2host[p->first + d * numHosts]);
      }
      ++received;
    }
    increment_evilPhase();
  }

  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host while taking
   * into consideration the edges of the node AND the node itself.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param numNodes_to_divide The total number of nodes you are
   * assigning to different hosts
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1. Ignored
   * in this function currently.
   *
   * @todo make this function work with decompose factor
   */
  void computeMastersBalancedNodesAndEdges(
      galois::graphs::OfflineGraph& g, uint64_t numNodes_to_divide,
      const std::vector<unsigned>& scalefactor, unsigned DecomposeFactor = 1) {
    if (nodeWeightOfMaster == 0) {
      nodeWeightOfMaster = g.sizeEdges() / g.size(); // average degree
    }

    if (edgeWeightOfMaster == 0) {
      edgeWeightOfMaster = 1;
    }
    auto& net = galois::runtime::getSystemNetworkInterface();
    gid2host.resize(numHosts);
    auto r = g.divideByNode(nodeWeightOfMaster, edgeWeightOfMaster, id,
                            numHosts, scalefactor);
    gid2host[id].first  = *r.first.first;
    gid2host[id].second = *r.first.second;
    for (unsigned h = 0; h < numHosts; ++h) {
      if (h == id)
        continue;
      galois::runtime::SendBuffer b;
      galois::runtime::gSerialize(b, gid2host[id]);
      net.sendTagged(h, galois::runtime::evilPhase, b);
    }
    net.flush();
    unsigned received = 1;
    while (received < numHosts) {
      decltype(net.recieveTagged(galois::runtime::evilPhase, nullptr)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase, nullptr);
      } while (!p);
      assert(p->first != id);
      auto& b = p->second;
      galois::runtime::gDeserialize(b, gid2host[p->first]);
      ++received;
    }
    increment_evilPhase();
  }

protected:
  /**
   * Wrapper call that will call into more specific compute masters
   * functions that compute masters based on nodes, edges, or both.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param isBipartite Specifies if the graph is a bipartite graph
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  uint64_t computeMasters(galois::graphs::OfflineGraph& g,
                          const std::vector<unsigned>& scalefactor,
                          bool isBipartite         = false,
                          unsigned DecomposeFactor = 1) {
    galois::Timer timer;
    timer.start();
    g.reset_seek_counters();

    uint64_t numNodes_to_divide = 0;

    if (isBipartite) {
      for (uint64_t n = 0; n < g.size(); ++n) {
        if (std::distance(g.edge_begin(n), g.edge_end(n))) {
          ++numNodes_to_divide;
          last_nodeID_withEdges_bipartite = n;
        }
      }
    } else {
      numNodes_to_divide = g.size();
    }

    // compute masters for all nodes
    switch (masters_distribution) {
    case BALANCED_MASTERS:
      computeMastersBlockedNodes(g, numNodes_to_divide, scalefactor,
                                 DecomposeFactor);
      break;
    case BALANCED_MASTERS_AND_EDGES:
      computeMastersBalancedNodesAndEdges(g, numNodes_to_divide, scalefactor,
                                          DecomposeFactor);
      break;
    case BALANCED_EDGES_OF_MASTERS:
    default:
      computeMastersBalancedEdges(g, numNodes_to_divide, scalefactor,
                                  DecomposeFactor);
      break;
    }

    timer.stop();

    galois::runtime::reportStatCond_Tmax<MORE_DIST_STATS>(
        GRNAME, "MasterDistTime", timer.get());

    galois::gPrint(
        "[", id, "] Master distribution time : ", timer.get_usec() / 1000000.0f,
        " seconds to read ", g.num_bytes_read(), " bytes in ", g.num_seeks(),
        " seeks (", g.num_bytes_read() / (float)timer.get_usec(), " MBPS)\n");
    return numNodes_to_divide;
  }

private:
  /**
   * Initalize MPI related things. The MPI layer itself should have been
   * initialized when the network interface was initiailized.
   */
  void initBareMPI() {
#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
    if (bare_mpi == noBareMPI)
      return;

#ifdef GALOIS_USE_LWCI
    // sanity check of ranks
    int taskRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &taskRank);
    if ((unsigned)taskRank != id)
      GALOIS_DIE("Mismatch in MPI rank");
    int numTasks;
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    if ((unsigned)numTasks != numHosts)
      GALOIS_DIE("Mismatch in MPI rank");
#endif

    // group setup
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    mpi_identity_groups.resize(numHosts);

    for (unsigned x = 0; x < numHosts; ++x) {
      const int g[1] = {(int)x};
      MPI_Group_incl(world_group, 1, g, &mpi_identity_groups[x]);
    }

    if (id == 0) {
      switch (bare_mpi) {
      case nonBlockingBareMPI:
        galois::gPrint("Using non-blocking bare MPI\n");
        break;
      case oneSidedBareMPI:
        galois::gPrint("Using one-sided bare MPI\n");
        break;
      case noBareMPI:
      default:
        GALOIS_DIE("Unsupported bare MPI");
      }
    }
#endif
  }

public:
  //! Type representing a node in this graph
  using GraphNode = typename GraphTy::GraphNode;
  //! iterator type over nodes
  using iterator = typename GraphTy::iterator;
  //! constant iterator type over nodes
  using const_iterator = typename GraphTy::const_iterator;
  //! iterator type over edges
  using edge_iterator = typename GraphTy::edge_iterator;

  /**
   * Constructor for DistGraph. Initializes metadata fields.
   *
   * @param host host number that this graph resides on
   * @param numHosts total number of hosts in the currently executing program
   */
  DistGraph(unsigned host, unsigned numHosts)
      : galois::runtime::GlobalObject(this), round(false), transposed(false),
        id(host), numHosts(numHosts) {
    enforce_data_mode = enforce_metadata;

    masterNodes.resize(numHosts);
    mirrorNodes.resize(numHosts);

    num_run        = 0;
    num_round      = 0;
    numGlobalEdges = 0;
    currentBVFlag  = nullptr;

    initBareMPI();
  }

protected:
  /**
   * Sets up the communication between the different hosts that contain
   * different parts of the graph by exchanging master/mirror information.
   */
  void setup_communication() {
    galois::CondStatTimer<MORE_DIST_STATS> Tcomm_setup("CommunicationSetupTime",
                                                       GRNAME);

    // barrier so that all hosts start the timer together
    galois::runtime::getHostBarrier().wait();

    Tcomm_setup.start();

    // Exchange information for memoization optimization.
    exchange_info_init();

    // convert the global ids stored in the master/mirror nodes arrays to local
    // ids
    // TODO: use 32-bit distinct vectors for masters and mirrors from here on
    for (uint32_t h = 0; h < masterNodes.size(); ++h) {
      galois::do_all(
          galois::iterate(0ul, masterNodes[h].size()),
          [&](uint32_t n) { masterNodes[h][n] = G2L(masterNodes[h][n]); },
#if MORE_COMM_STATS
          galois::loopname(get_run_identifier("MasterNodes").c_str()),
#endif
          galois::no_stats());
    }

    for (uint32_t h = 0; h < mirrorNodes.size(); ++h) {
      galois::do_all(
          galois::iterate(0ul, mirrorNodes[h].size()),
          [&](uint32_t n) { mirrorNodes[h][n] = G2L(mirrorNodes[h][n]); },
#if MORE_COMM_STATS
          galois::loopname(get_run_identifier("MirrorNodes").c_str()),
#endif
          galois::no_stats());
    }

    Tcomm_setup.stop();

    // report masters/mirrors to/from other hosts as statistics
    for (auto x = 0U; x < masterNodes.size(); ++x) {
      if (x == id)
        continue;
      std::string master_nodes_str =
          "MasterNodesFrom_" + std::to_string(id) + "_To_" + std::to_string(x);
      galois::runtime::reportStatCond_Tsum<MORE_DIST_STATS>(
          GRNAME, master_nodes_str, masterNodes[x].size());
    }

    for (auto x = 0U; x < mirrorNodes.size(); ++x) {
      if (x == id)
        continue;
      std::string mirror_nodes_str =
          "MirrorNodesFrom_" + std::to_string(x) + "_To_" + std::to_string(id);
      galois::runtime::reportStatCond_Tsum<MORE_DIST_STATS>(
          GRNAME, mirror_nodes_str, mirrorNodes[x].size());
    }

    send_info_to_host();

    // do not track memory usage of partitioning
    auto& net = galois::runtime::getSystemNetworkInterface();
    net.resetMemUsage();
  }

public:
  /**
   * Do an in-memory transpose while keeping the original graph intact.
   * Find thread ranges as well.
   *
   * Only does something if WithInEdges is enabled.
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  void constructIncomingEdges() {
    graph.constructIncomingEdges();

    determineThreadRangesIn();
    determineThreadRangesMasterIn();
    initializeSpecificRangesIn();
  }

  /**
   * Does nothing (in other words WithInEdges is false).
   */
  template <bool T = WithInEdges, typename std::enable_if<!T>::type* = nullptr>
  void constructIncomingEdges() {
    // do nothing since the incoming edges template argument is false
    return;
  }

  /**
   * Get data of a node.
   *
   * @param N node to get the data of
   * @param mflag access flag for node data
   * @returns A node data object
   */
  inline NodeTy&
  getData(GraphNode N,
          galois::MethodFlag mflag = galois::MethodFlag::UNPROTECTED) {
    auto& r = graph.getData(N, mflag);
    return r;
  }

  /**
   * Get the node data for a particular node in the graph.
   *
   * @param ni edge to get the data of
   * @param mflag access flag for edge data
   * @returns The edge data for the requested edge
   */
  inline typename GraphTy::edge_data_reference
  getEdgeData(edge_iterator ni,
              galois::MethodFlag mflag = galois::MethodFlag::UNPROTECTED) {
    auto& r = graph.getEdgeData(ni, mflag);
    return r;
  }

  /**
   * Directly gets in-edge data.
   *
   * @param ni in-edge to get the data of
   * @param mflag access flag for edge data
   * @returns The edge data for the requested in-edge
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline typename GraphTy::edge_data_reference getInEdgeData(
      edge_iterator ni,
      galois::MethodFlag mflag = galois::MethodFlag::UNPROTECTED) const {
    return graph.getInEdgeData(ni, mflag);
  }

  /**
   * Gets edge destination of edge ni.
   *
   * @param ni edge id to get destination of
   * @returns Local ID of destination of edge ni
   */
  GraphNode getEdgeDst(edge_iterator ni) { return graph.getEdgeDst(ni); }

  /**
   * Gets edge destination of in-edge ni.
   *
   * @param ni in-edge id to get destination of
   * @returns Local ID of destination of in-edge ni
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  GraphNode getInEdgeDst(edge_iterator ni) const {
    return graph.getInEdgeDst(ni);
  }

  /**
   * Gets the first edge of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to first edge of N
   */
  inline edge_iterator edge_begin(GraphNode N) {
    return graph.edge_begin(N, galois::MethodFlag::UNPROTECTED);
  }

  /**
   * Gets the end edge boundary of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to the end of the edges of node N, i.e. the first edge
   * of the next node (or an "end" iterator if there is no next node)
   */
  inline edge_iterator edge_end(GraphNode N) {
    return graph.edge_end(N, galois::MethodFlag::UNPROTECTED);
  }

  /**
   * Gets the first in-edge of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to first in-edge of N
   */
  // template<bool T = WithInEdges, typename std::enable_if<T == true>::type*>
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline edge_iterator in_edge_begin(GraphNode N) {
    return graph.in_edge_begin(N, galois::MethodFlag::UNPROTECTED);
  }

  /**
   * Gets the end in-edge boundary of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to the end of the in-edges of node N, i.e. the first
   * in-edge of the next node (or an "end" iterator if there is no next node)
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline edge_iterator in_edge_end(GraphNode N) {
    return graph.in_edge_end(N, galois::MethodFlag::UNPROTECTED);
  }

  /**
   * Returns an iterable object over the edges of a particular node in the
   * graph.
   *
   * @param N node to get edges iterator over
   */
  inline galois::runtime::iterable<galois::NoDerefIterator<edge_iterator>>
  edges(GraphNode N) {
    return galois::graphs::internal::make_no_deref_range(edge_begin(N),
                                                         edge_end(N));
  }

  /**
   * Returns an iterable object over the in-edges of a particular node in the
   * graph.
   *
   * @param N node to get edges iterator over
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline galois::runtime::iterable<galois::NoDerefIterator<edge_iterator>>
  in_edges(GraphNode N) {
    return galois::graphs::internal::make_no_deref_range(in_edge_begin(N),
                                                         in_edge_end(N));
  }

  /**
   * Gets number of nodes on this (local) graph.
   *
   * @returns number of nodes present in this (local) graph
   */
  inline size_t size() const { return graph.size(); }

  /**
   * Gets number of edges on this (local) graph.
   *
   * @returns number of edges present in this (local) graph
   */
  inline size_t sizeEdges() const { return graph.sizeEdges(); }

  /**
   * Gets number of nodes on the global unpartitioned graph.
   *
   * @returns number of nodes present in the global unpartitioned graph
   */
  inline size_t globalSize() const { return numGlobalNodes; }

  /**
   * Gets number of edges on the global unpartitioned graph.
   *
   * @returns number of edges present in the global unpartitioned graph
   */
  inline size_t globalSizeEdges() const { return numGlobalEdges; }

  /**
   * Returns a range object that encapsulates all nodes of the graph.
   *
   * @returns A range object that contains all the nodes in this graph
   */
  inline const NodeRangeType& allNodesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[0];
  }

  /**
   * Returns a range object that encapsulates only master nodes in this
   * graph.
   *
   * @returns A range object that contains the master nodes in this graph
   */
  inline const NodeRangeType& masterNodesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[1];
  }

  /**
   * Returns a range object that encapsulates master nodes and nodes
   * with edges in this graph.
   *
   * @returns A range object that contains the master nodes and the nodes
   * with outgoing edges in this graph
   */
  inline const NodeRangeType& allNodesWithEdgesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[2];
  }

  /**
   * Returns a range object that encapsulates all nodes of the graph;
   * specific range based on in-edge distribution.
   *
   * @returns A range object that contains all the nodes in this graph
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline const NodeRangeType& allNodesRangeIn() const {
    return specificRangesIn[0];
  }

  /**
   * Returns a range object that encapsulates only master nodes in this
   * graph; specific range based on in-edge distribution.
   *
   * @returns A range object that contains the master nodes in this graph
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline const NodeRangeType& masterNodesRangeIn() const {
    return specificRangesIn[1];
  }

  /**
   * Returns a range object that encapsulates master nodes and nodes
   * with in-edges in this graph; specific range based on in-edge distribution.
   *
   * @returns A range object that contains the master nodes and the nodes
   * with incoming edges in this graph (i.e. all nodes)
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline const NodeRangeType& allNodesWithEdgesRangeIn() const {
    return specificRangesIn[0];
  }

protected:
  /**
   * Uses a pre-computed prefix sum to determine division of nodes among
   * threads.
   *
   * The call uses binary search to determine the ranges.
   */
  inline void determineThreadRanges() {
    allNodesRanges = galois::graphs::determineUnitRangesFromPrefixSum(
        galois::runtime::activeThreads, graph.getEdgePrefixSum());
  }

  /**
   * Determines the thread ranges for master nodes only and saves them to
   * the object.
   *
   * Only call after graph is constructed + only call once
   */
  inline void determineThreadRangesMaster() {
    // make sure this hasn't been called before
    assert(masterRanges.size() == 0);

    // first check if we even need to do any work; if already calculated,
    // use already calculated vector
    if (beginMaster == 0 && (beginMaster + numOwned) == size()) {
      masterRanges = allNodesRanges;
    } else if (beginMaster == 0 &&
               (beginMaster + numOwned) == numNodesWithEdges &&
               withEdgeRanges.size() != 0) {
      masterRanges = withEdgeRanges;
    } else {
      galois::gDebug("Manually det. master thread ranges");
      masterRanges = galois::graphs::determineUnitRangesFromGraph(
          graph, galois::runtime::activeThreads, beginMaster,
          beginMaster + numOwned, nodeAlphaRanges);
    }
  }

  /**
   * Determines the thread ranges for nodes with edges only and saves them to
   * the object.
   *
   * Only call after graph is constructed + only call once
   */
  inline void determineThreadRangesWithEdges() {
    // make sure not called before
    assert(withEdgeRanges.size() == 0);

    // first check if we even need to do any work; if already calculated,
    // use already calculated vector
    if (numNodesWithEdges == size()) {
      withEdgeRanges = allNodesRanges;
    } else if (beginMaster == 0 &&
               (beginMaster + numOwned) == numNodesWithEdges &&
               masterRanges.size() != 0) {
      withEdgeRanges = masterRanges;
    } else {
      galois::gDebug("Manually det. with edges thread ranges");
      withEdgeRanges = galois::graphs::determineUnitRangesFromGraph(
          graph, galois::runtime::activeThreads, 0, numNodesWithEdges,
          nodeAlphaRanges);
    }
  }

  /**
   * Initializes the 3 range objects that a user can access to iterate
   * over the graph in different ways.
   */
  void initializeSpecificRanges() {
    assert(specificRanges.size() == 0);

    // TODO/FIXME assertion likely not safe if a host gets no nodes
    // make sure the thread ranges have already been calculated
    // for the 3 ranges
    assert(allNodesRanges.size() != 0);
    assert(masterRanges.size() != 0);
    assert(withEdgeRanges.size() != 0);

    // 0 is all nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(size()), allNodesRanges.data()));
    allNodesRanges.clear();

    // 1 is master nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(beginMaster),
        boost::counting_iterator<size_t>(beginMaster + numOwned),
        masterRanges.data()));
    masterRanges.clear();

    // 2 is with edge nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(numNodesWithEdges),
        withEdgeRanges.data()));
    withEdgeRanges.clear();

    assert(specificRanges.size() == 3);
  }

  /**
   * Specific range editor: makes the range for edges equivalent to the range
   * for masters.
   */
  void edgesEqualMasters() { specificRanges[2] = specificRanges[1]; }

  /**
   * Uses a pre-computed prefix sum to determine division of nodes among
   * threads using in-edges.
   *
   * The call uses binary search to determine the ranges.
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline void determineThreadRangesIn() {
    galois::gDebug("Determining thread ranges in");
    allNodesRangesIn = galois::graphs::determineUnitRangesFromPrefixSum(
        galois::runtime::activeThreads, graph.getInEdgePrefixSum());
  }

  /**
   * Determine the thread ranges (splitting of nodes among threads) of master
   * nodes using in-edges.
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  inline void determineThreadRangesMasterIn() {
    assert(masterRangesIn.size() == 0);
    galois::gDebug("Determining master thread ranges in");

    // first check if we even need to do any work; if already calculated,
    // use already calculated vector
    if (beginMaster == 0 && (beginMaster + numOwned) == graph.size()) {
      masterRangesIn = allNodesRangesIn;
    } else {
      galois::gDebug("Manually det. master thread in ranges");
      masterRangesIn = galois::graphs::determineUnitRangesFromPrefixSum(
          galois::runtime::activeThreads, graph.getInEdgePrefixSum(),
          beginMaster, beginMaster + numOwned, nodeAlphaRanges);
    }
  }

  /**
   * Initialize the thread ranges for in-edges.
   */
  template <bool T = WithInEdges, typename std::enable_if<T>::type* = nullptr>
  void initializeSpecificRangesIn() {
    assert(specificRangesIn.size() == 0);
    galois::gDebug("Initialize specific ranges in");

    // 0 is all nodes
    specificRangesIn.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(size()), allNodesRangesIn.data()));

    // 1 is master nodes
    specificRangesIn.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(beginMaster),
        boost::counting_iterator<size_t>(beginMaster + numOwned),
        masterRangesIn.data()));

    assert(specificRangesIn.size() == 2);
  }

private:
  /**
   * Let other hosts know about which host has what mirrors/masters;
   * used for later communication of mirrors/masters.
   */
  void exchange_info_init() {
    auto& net = galois::runtime::getSystemNetworkInterface();

    // send off the mirror nodes
    for (unsigned x = 0; x < numHosts; ++x) {
      if (x == id)
        continue;

      galois::runtime::SendBuffer b;
      gSerialize(b, mirrorNodes[x]);
      net.sendTagged(x, galois::runtime::evilPhase, b);
    }

    // receive the mirror nodes
    for (unsigned x = 0; x < numHosts; ++x) {
      if (x == id)
        continue;

      decltype(net.recieveTagged(galois::runtime::evilPhase, nullptr)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase, nullptr);
      } while (!p);

      galois::runtime::gDeserialize(p->second, masterNodes[p->first]);
    }
    increment_evilPhase();
  }

  /**
   * Reports master/mirror stats.
   * Assumes that communication has already occured so that the host
   * calling it actually has the info required.
   *
   * @param global_total_mirror_nodes number of mirror nodes on all hosts
   * @param global_total_owned_nodes number of "owned" nodes on all hosts
   */
  void report_master_mirror_stats(uint64_t global_total_mirror_nodes,
                                  uint64_t global_total_owned_nodes) {
    float replication_factor =
        (float)(global_total_mirror_nodes + numGlobalNodes) /
        (float)numGlobalNodes;
    galois::runtime::reportStat_Single(GRNAME, "ReplicationFactor",
                                       replication_factor);

    galois::runtime::reportStatCond_Single<MORE_DIST_STATS>(
        GRNAME, "TotalNodes", numGlobalNodes);
    galois::runtime::reportStatCond_Single<MORE_DIST_STATS>(
        GRNAME, "TotalGlobalMirrorNodes", global_total_mirror_nodes);
  }

  /**
   * Send statistics about master/mirror nodes to each host, and
   * report the statistics.
   */
  void send_info_to_host() {
    auto& net = galois::runtime::getSystemNetworkInterface();

    uint64_t global_total_mirror_nodes = size() - numOwned;
    uint64_t global_total_owned_nodes  = numOwned;

    // send info to host
    for (unsigned x = 0; x < numHosts; ++x) {
      if (x == id)
        continue;

      galois::runtime::SendBuffer b;
      gSerialize(b, global_total_mirror_nodes, global_total_owned_nodes);
      net.sendTagged(x, galois::runtime::evilPhase, b);
    }

    // receive
    for (unsigned x = 0; x < numHosts; ++x) {
      if (x == id)
        continue;

      decltype(net.recieveTagged(galois::runtime::evilPhase, nullptr)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase, nullptr);
      } while (!p);

      uint64_t total_mirror_nodes_from_others;
      uint64_t total_owned_nodes_from_others;
      galois::runtime::gDeserialize(p->second, total_mirror_nodes_from_others,
                                    total_owned_nodes_from_others);
      global_total_mirror_nodes += total_mirror_nodes_from_others;
      global_total_owned_nodes += total_owned_nodes_from_others;
    }
    increment_evilPhase();

    assert(numGlobalNodes == global_total_owned_nodes);
    // report stats
    if (net.ID == 0) {
      report_master_mirror_stats(global_total_mirror_nodes,
                                 global_total_owned_nodes);
    }
  }

  /**
   * Given a bitset, determine the indices of the bitset that are currently
   * set.
   *
   * @tparam syncType either reduce or broadcast; only used to name the timer
   *
   * @param loopName string used to name the timer for this function
   * @param bitset_comm the bitset to get the offsets of
   * @param offsets output: the offset vector that will contain indices into
   * the bitset that are set
   * @param bit_set_count output: will be set to the number of bits set in the
   * bitset
   */
  template <SyncType syncType>
  void get_offsets_from_bitset(const std::string& loopName,
                               const galois::DynamicBitSet& bitset_comm,
                               std::vector<unsigned int>& offsets,
                               size_t& bit_set_count) const {
    // timer creation
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string offsets_timer_str(syncTypeStr + "Offsets_" +
                                  get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Toffsets(offsets_timer_str.c_str(),
                                                    GRNAME);

    Toffsets.start();

    auto activeThreads = galois::getActiveThreads();
    std::vector<unsigned int> t_prefix_bit_counts(activeThreads);

    // count how many bits are set on each thread
    galois::on_each([&](unsigned tid, unsigned nthreads) {
      // TODO use block_range instead
      unsigned int block_size = bitset_comm.size() / nthreads;
      if ((bitset_comm.size() % nthreads) > 0)
        ++block_size;
      assert((block_size * nthreads) >= bitset_comm.size());

      unsigned int start = tid * block_size;
      unsigned int end   = (tid + 1) * block_size;
      if (end > bitset_comm.size())
        end = bitset_comm.size();

      unsigned int count = 0;
      for (unsigned int i = start; i < end; ++i) {
        if (bitset_comm.test(i))
          ++count;
      }

      t_prefix_bit_counts[tid] = count;
    });

    // calculate prefix sum of bits per thread
    for (unsigned int i = 1; i < activeThreads; ++i) {
      t_prefix_bit_counts[i] += t_prefix_bit_counts[i - 1];
    }
    // total num of set bits
    bit_set_count = t_prefix_bit_counts[activeThreads - 1];

    // calculate the indices of the set bits and save them to the offset
    // vector
    if (bit_set_count > 0) {
      offsets.resize(bit_set_count);
      galois::on_each([&](unsigned tid, unsigned nthreads) {
        // TODO use block_range instead
        // TODO this is same calculation as above; maybe refactor it
        // into function?
        unsigned int block_size = bitset_comm.size() / nthreads;
        if ((bitset_comm.size() % nthreads) > 0)
          ++block_size;
        assert((block_size * nthreads) >= bitset_comm.size());

        unsigned int start = tid * block_size;
        unsigned int end   = (tid + 1) * block_size;
        if (end > bitset_comm.size())
          end = bitset_comm.size();

        unsigned int count = 0;
        unsigned int t_prefix_bit_count;
        if (tid == 0) {
          t_prefix_bit_count = 0;
        } else {
          t_prefix_bit_count = t_prefix_bit_counts[tid - 1];
        }

        for (unsigned int i = start; i < end; ++i) {
          if (bitset_comm.test(i)) {
            offsets[t_prefix_bit_count + count] = i;
            ++count;
          }
        }
      });
    }
    Toffsets.stop();
  }

  /**
   * Determine what data needs to be synchronized based on the passed in
   * bitset_compute and returns information regarding these need-to-be-sync'd
   * nodes.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done;
   * only used to get the size of the type being synchronized in this function
   * @tparam syncType type of synchronization this function is being called
   * for; only used to name a timer
   *
   * @param loopName loopname used to name the timer for the function
   * @param indices A vector that contains the local ids of the nodes that
   * you want to potentially synchronize
   * @param bitset_compute Contains the full bitset of all nodes in this
   * graph
   * @param bitset_comm OUTPUT: bitset that marks which indices in the passed
   * in indices array need to be synchronized
   * @param offsets OUTPUT: contains indices into bitset_comm that are set
   * @param bit_set_count OUTPUT: contains number of bits set in bitset_comm
   * @param data_mode OUTPUT: the way that this data should be communicated
   * based on how much data needs to be sent out
   */
  template <typename FnTy, SyncType syncType>
  void get_bitset_and_offsets(const std::string& loopName,
                              const std::vector<size_t>& indices,
                              const galois::DynamicBitSet& bitset_compute,
                              galois::DynamicBitSet& bitset_comm,
                              std::vector<unsigned int>& offsets,
                              size_t& bit_set_count,
                              DataCommMode& data_mode) const {
    if (enforce_data_mode != onlyData) {
      bitset_comm.reset();
      std::string syncTypeStr =
          (syncType == syncReduce) ? "Reduce" : "Broadcast";
      std::string doall_str(syncTypeStr + "Bitset_" + loopName);

      // determine which local nodes in the indices array need to be
      // sychronized
      galois::do_all(galois::iterate(0ul, indices.size()),
                     [&](unsigned int n) {
                       // assumes each lid is unique as test is not thread safe
                       size_t lid = indices[n];
                       if (bitset_compute.test(lid)) {
                         bitset_comm.set(n);
                       }
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());

      // get the number of set bits and the offsets into the comm bitset
      get_offsets_from_bitset<syncType>(loopName, bitset_comm, offsets,
                                        bit_set_count);
    }

    data_mode =
        get_data_mode<typename FnTy::ValTy>(bit_set_count, indices.size());
  }

  /**
   * Extracts data at provided lid.
   *
   * This version (reduce) resets the value after extract.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; determines if reset is
   * necessary
   *
   * @param lid local id of node to get data from
   * @returns data (specified by FnTy) of node with local id lid
   */
  /* Reduction extract resets the value afterwards */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline typename FnTy::ValTy extract_wrapper(size_t lid) {
#ifdef __GALOIS_HET_OPENCL__
    CLNodeDataWrapper d = clGraph.getDataW(lid);
    auto val            = FnTy::extract(lid, getData(lid, d));
    FnTy::reset(lid, d);
#else
    auto val = FnTy::extract(lid, getData(lid));
    FnTy::reset(lid, getData(lid));
#endif
    return val;
  }

  /**
   * Extracts data at provided lid; uses vecIndex to get the correct element
   * from the vector.
   *
   * This version (reduce) resets the value after extract.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; determines if reset is
   * necessary
   *
   * @param lid local id of node to get data from
   * @param vecIndex index to grab from vector in node
   * @returns data (specified by FnTy) of node with local id lid
   */
  /* Reduction extract resets the value afterwards */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline typename FnTy::ValTy extract_wrapper(size_t lid, unsigned vecIndex) {
    auto val = FnTy::extract(lid, getData(lid), vecIndex);
    FnTy::reset(lid, getData(lid), vecIndex);
    return val;
  }

  /**
   * Extracts data at provided lid.
   *
   * This version (broadcast) does not reset the value after extract.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; determines if reset is
   * necessary
   *
   * @param lid local id of node to get data from
   * @returns data (specified by FnTy) of node with local id lid
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline typename FnTy::ValTy extract_wrapper(size_t lid) {
#ifdef __GALOIS_HET_OPENCL__
    CLNodeDataWrapper d = clGraph.getDataW(lid);
    return FnTy::extract(lid, getData(lid, d));
#else
    return FnTy::extract(lid, getData(lid));
#endif
  }

  /**
   * Extracts data at provided lid; uses vecIndex to get the correct element
   * from the vector in the node.
   *
   * This version (broadcast) does not reset the value after extract.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; determines if reset is
   * necessary
   *
   * @param lid local id of node to get data from
   * @param vecIndex index to grab from vector in node
   * @returns data (specified by FnTy) of node with local id lid
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline typename FnTy::ValTy extract_wrapper(size_t lid, unsigned vecIndex) {
    return FnTy::extract(lid, getData(lid), vecIndex);
  }

  /**
   * Based on provided arguments, extracts the data that we are interested
   * in sending into val_vec.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; used to determine if reseting
   * the extracted field is necessary
   * @tparam identity_offsets If this is true, then ignore the offsets
   * array and just grab directly from indices (i.e. don't pick out
   * particular elements, just grab contiguous chunk)
   * @tparam parallelize Determines if parallelizing the extraction is done or
   * not
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param size Number of elements to extract
   * @param offsets Holds offsets into "indices" of the data that we are
   * interested in
   * @param val_vec OUTPUT: holds the extracted data
   * @param start Offset into val_vec to start saving data to
   */
  template <typename FnTy, SyncType syncType, bool identity_offsets = false,
            bool parallelize = true>
  void extract_subset(const std::string& loopName,
                      const std::vector<size_t>& indices, size_t size,
                      const std::vector<unsigned int>& offsets,
                      std::vector<typename FnTy::ValTy>& val_vec,
                      size_t start = 0) {
    if (parallelize) {
      std::string syncTypeStr =
          (syncType == syncReduce) ? "Reduce" : "Broadcast";
      std::string doall_str(syncTypeStr + "ExtractVal_" + loopName);

      galois::do_all(galois::iterate(start, start + size),
                     [&](unsigned int n) {
                       unsigned int offset;
                       if (identity_offsets)
                         offset = n;
                       else
                         offset = offsets[n];

                       size_t lid = indices[offset];
                       val_vec[n - start] =
                           extract_wrapper<FnTy, syncType>(lid);
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());
    } else { // non-parallel version
      for (unsigned n = start; n < start + size; ++n) {
        unsigned int offset;
        if (identity_offsets)
          offset = n;
        else
          offset = offsets[n];

        size_t lid         = indices[offset];
        val_vec[n - start] = extract_wrapper<FnTy, syncType>(lid);
      }
    }
  }

  /**
   * Based on provided arguments, extracts the data that we are interested
   * in sending into val_vec. Same as above, except it has the vecIndex
   * arguments and requires vecSync to be true
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType either reduce or broadcast; used to determine if reseting
   * the extracted field is necessary
   * @tparam identity_offsets If this is true, then ignore the offsets
   * array and just grab directly from indices (i.e. don't pick out
   * particular elements, just grab contiguous chunk)
   * @tparam parallelize Determines if parallelizing the extraction is done or
   * not
   * @tparam vecSync Only set to true if the field being synchronized is a
   * vector and synchronization is occuring element by element. MUST BE SET
   * TO TRUE IN ORDER FOR THIS FUNCTION TO COMPILE.
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param size Number of elements to extract
   * @param offsets Holds offsets into "indices" of the data that we are
   * interested in
   * @param val_vec OUTPUT: holds the extracted data
   * @param vecIndex which element of the vector to extract from node
   * @param start Offset into val_vec to start saving data to
   */
  // TODO find a better way to have this variant without code duplication
  template <typename FnTy, SyncType syncType, bool identity_offsets = false,
            bool parallelize = true, bool vecSync = false,
            typename std::enable_if<vecSync>::type* = nullptr>
  void extract_subset(const std::string& loopName,
                      const std::vector<size_t>& indices, size_t size,
                      const std::vector<unsigned int>& offsets,
                      std::vector<typename FnTy::ValTy>& val_vec,
                      unsigned vecIndex, size_t start = 0) {
    val_vec.resize(size); // resive val vec for this vecIndex

    if (parallelize) {
      std::string syncTypeStr =
          (syncType == syncReduce) ? "Reduce" : "Broadcast";
      std::string doall_str(syncTypeStr + "ExtractValVector_" + loopName);

      galois::do_all(galois::iterate(start, start + size),
                     [&](unsigned int n) {
                       unsigned int offset;
                       if (identity_offsets)
                         offset = n;
                       else
                         offset = offsets[n];

                       size_t lid = indices[offset];
                       val_vec[n - start] =
                           extract_wrapper<FnTy, syncType>(lid, vecIndex);
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());
    } else { // non-parallel version
      for (unsigned n = start; n < start + size; ++n) {
        unsigned int offset;
        if (identity_offsets)
          offset = n;
        else
          offset = offsets[n];

        size_t lid         = indices[offset];
        val_vec[n - start] = extract_wrapper<FnTy, syncType>(lid, vecIndex);
      }
    }
  }

  /**
   * Based on provided arguments, extracts the data that we are interested
   * in sending into a send buffer. Lazy serialize variant that works with
   * certain SeqTy.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SeqTy Type of sequence that we are getting data from
   * @tparam syncType either reduce or broadcast; used to determine if reseting
   * the extracted field is necessary
   * @tparam identity_offsets If this is true, then ignore the offsets
   * array and just grab directly from indices (i.e. don't pick out
   * particular elements, just grab contiguous chunk)
   * @tparam parallelize Determines if parallelizing the extraction is done or
   * not
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param size Number of elements to extract
   * @param offsets Holds offsets into "indices" of the data that we are
   * interested in
   * @param b send buffer to extract data into
   * @param lseq sequence to get data from
   * @param start Offset into send buffer to start saving data to
   */
  template <typename FnTy, typename SeqTy, SyncType syncType,
            bool identity_offsets = false, bool parallelize = true>
  void extract_subset(const std::string& loopName,
                      const std::vector<size_t>& indices, size_t size,
                      const std::vector<unsigned int>& offsets,
                      galois::runtime::SendBuffer& b, SeqTy lseq,
                      size_t start = 0) {
    if (parallelize) {
      std::string syncTypeStr =
          (syncType == syncReduce) ? "Reduce" : "Broadcast";
      std::string doall_str(syncTypeStr + "ExtractVal_" + loopName);

      galois::do_all(galois::iterate(start, start + size),
                     [&](unsigned int n) {
                       unsigned int offset;

                       if (identity_offsets)
                         offset = n;
                       else
                         offset = offsets[n];

                       size_t lid = indices[offset];
                       gSerializeLazy(b, lseq, n - start,
                                      extract_wrapper<FnTy, syncType>(lid));
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());
    } else {
      for (unsigned int n = start; n < start + size; ++n) {
        unsigned int offset;

        if (identity_offsets)
          offset = n;
        else
          offset = offsets[n];

        size_t lid = indices[offset];
        gSerializeLazy(b, lseq, n - start,
                       extract_wrapper<FnTy, syncType>(lid));
      }
    }
  }

  /**
   * GPU wrap function: extracts data from nodes and resets them to the
   * reduction identity value as specified by the sync structure. (Reduce only)
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be reduce
   *
   * @param x node id to extract from
   * @param v vector to extract data to
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline bool extract_batch_wrapper(unsigned x,
                                    std::vector<typename FnTy::ValTy>& v) {
    return FnTy::extract_reset_batch(x, v.data());
  }

  /**
   * GPU wrap function: extracts data from nodes. (Broadcast only)
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be broadcast
   *
   * @param x node id to extract from
   * @param v vector to extract data to
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline bool extract_batch_wrapper(unsigned x,
                                    std::vector<typename FnTy::ValTy>& v) {
    return FnTy::extract_batch(x, v.data());
  }

  /**
   * GPU wrap function: extracts data from nodes and resets them to the
   * reduction identity value as specified by the sync structure. (Reduce only)
   *
   * This version specifies more arguments.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be reduce
   *
   * @param x node id to extract from
   * @param b
   * @param o
   * @param v
   * @param s
   * @param data_mode
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline bool extract_batch_wrapper(unsigned x, galois::DynamicBitSet& b,
                                    std::vector<unsigned int>& o,
                                    std::vector<typename FnTy::ValTy>& v,
                                    size_t& s, DataCommMode& data_mode) {
    return FnTy::extract_reset_batch(x, (uint64_t*)b.get_vec().data(), o.data(),
                                     v.data(), &s, &data_mode);
  }

  /**
   * GPU wrap function: extracts data from nodes (Broadcast only)
   *
   * This version specifies more arguments.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be broadcast
   *
   * @param x node id to extract from
   * @param b
   * @param o
   * @param v
   * @param s
   * @param data_mode
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline bool extract_batch_wrapper(unsigned x, galois::DynamicBitSet& b,
                                    std::vector<unsigned int>& o,
                                    std::vector<typename FnTy::ValTy>& v,
                                    size_t& s, DataCommMode& data_mode) const {
    return FnTy::extract_batch(x, (uint64_t*)b.get_vec().data(), o.data(),
                               v.data(), &s, &data_mode);
  }

  /**
   * Reduce variant. Takes a value and reduces it according to the sync
   * structure provided to the function.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType Reduce sync or broadcast sync
   *
   * @param lid local id of node to reduce to
   * @param val value to reduce to
   * @param bit_set_compute bitset indicating which nodes have changed; updated
   * if reduction causes a change
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline void set_wrapper(size_t lid, typename FnTy::ValTy val,
                          galois::DynamicBitSet& bit_set_compute) {
#ifdef __GALOIS_HET_OPENCL__
    CLNodeDataWrapper d = clGraph.getDataW(lid);
    FnTy::reduce(lid, d, val);
#else
    if (FnTy::reduce(lid, getData(lid), val)) {
      if (bit_set_compute.size() != 0)
        bit_set_compute.set(lid);
    }
#endif
  }

  /**
   * VECTOR VARIANT.
   *
   * Reduce variant. Takes a value and reduces it according to the sync
   * structure provided to the function. Only reduces the element at a
   * particular index of the vector field being sychronized.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType Reduce sync or broadcast sync
   *
   * @param lid local id of node to reduce to
   * @param val value to reduce to
   * @param bit_set_compute bitset indicating which nodes have changed; updated
   * if reduction causes a change
   * @param vecIndex which element of the vector to reduce in the node
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline void set_wrapper(size_t lid, typename FnTy::ValTy val,
                          galois::DynamicBitSet& bit_set_compute,
                          unsigned vecIndex) {
    if (FnTy::reduce(lid, getData(lid), val, vecIndex)) {
      if (bit_set_compute.size() != 0)
        bit_set_compute.set(lid);
    }
  }

  /**
   * Broadcast variant. Takes a value and sets it according to the sync
   * structure provided to the function.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType Reduce sync or broadcast sync
   *
   * @param lid local id of node to reduce to
   * @param val value to reduce to
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline void set_wrapper(size_t lid, typename FnTy::ValTy val,
                          galois::DynamicBitSet&) {
#ifdef __GALOIS_HET_OPENCL__
    CLNodeDataWrapper d = clGraph.getDataW(lid);
    FnTy::setVal(lid, d, val_vec[n]);
#else
    FnTy::setVal(lid, getData(lid), val);
#endif
  }

  /**
   * VECTOR VARIANT.
   *
   * Broadcast variant. Takes a value and sets it according to the sync
   * structure provided to the function. Only sets the element at the specified
   * index of the vector in the node.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam syncType Reduce sync or broadcast sync
   *
   * @param lid local id of node to reduce to
   * @param val value to reduce to
   * @param vecIndex which element of the vector to set in the node
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline void set_wrapper(size_t lid, typename FnTy::ValTy val,
                          galois::DynamicBitSet&, unsigned vecIndex) {
    FnTy::setVal(lid, getData(lid), val, vecIndex);
  }

  /**
   * Given data received from another host and information on which nodes
   * to update, do the reduce/set of the received data to update local nodes.
   *
   * Complement function, in some sense, of extract_subset.
   *
   * @tparam VecTy type of indices variable
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Reduce or broadcast
   * @tparam identity_offsets If this is true, then ignore the offsets
   * array and just grab directly from indices (i.e. don't pick out
   * particular elements, just grab contiguous chunk)
   * @tparam parallelize True if updates to nodes are to be parallelized
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param size Number of elements to set
   * @param offsets Holds offsets into "indices" of the data that we are
   * interested in
   * @param val_vec holds data we will use to set
   * @param bit_set_compute bitset indicating which nodes have changed
   * @param start Offset into val_vec to start saving data to
   */
  template <typename VecTy, typename FnTy, SyncType syncType,
            bool identity_offsets = false, bool parallelize = true>
  void set_subset(const std::string& loopName, const VecTy& indices,
                  size_t size, const std::vector<unsigned int>& offsets,
                  std::vector<typename FnTy::ValTy>& val_vec,
                  galois::DynamicBitSet& bit_set_compute, size_t start = 0) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string doall_str(syncTypeStr + "SetVal_" +
                          get_run_identifier(loopName));

    if (parallelize) {
      galois::do_all(galois::iterate(start, start + size),
                     [&](unsigned int n) {
                       unsigned int offset;

                       if (identity_offsets)
                         offset = n;
                       else
                         offset = offsets[n];

                       auto lid = indices[offset];
                       set_wrapper<FnTy, syncType>(lid, val_vec[n - start],
                                                   bit_set_compute);
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());
    } else {
      for (unsigned int n = start; n < start + size; ++n) {
        unsigned int offset;

        if (identity_offsets)
          offset = n;
        else
          offset = offsets[n];

        auto lid = indices[offset];
        set_wrapper<FnTy, syncType>(lid, val_vec[n - start], bit_set_compute);
      }
    }
  }

  /**
   * VECTOR BITSET VARIANT.
   *
   * Given data received from another host and information on which nodes
   * to update, do the reduce/set of the received data to update local nodes.
   * It will only update a single index of the vector specified by the
   * sync structures at a time.
   *
   * Complement function, in some sense, of extract_subset, vector bitset
   * variant.
   *
   * @tparam VecTy type of indices variable
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Reduce or broadcast
   * @tparam identity_offsets If this is true, then ignore the offsets
   * array and just grab directly from indices (i.e. don't pick out
   * particular elements, just grab contiguous chunk)
   * @tparam parallelize True if updates to nodes are to be parallelized
   * @tparam vecSync Only set to true if the field being synchronized is a
   * vector. MUST BE SET TO TRUE FOR THIS FUNCTION TO COMPILE
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param size Number of elements to set
   * @param offsets Holds offsets into "indices" of the data that we are
   * interested in
   * @param val_vec holds data we will use to set
   * @param bit_set_compute bitset indicating which nodes have changed
   * @param vecIndex which element of the vector to set in the node
   * @param start Offset into val_vec to start saving data to
   */
  // TODO find a better way to have this variant without code duplication
  template <typename VecTy, typename FnTy, SyncType syncType,
            bool identity_offsets = false, bool parallelize = true,
            bool vecSync                            = false,
            typename std::enable_if<vecSync>::type* = nullptr>
  void set_subset(const std::string& loopName, const VecTy& indices,
                  size_t size, const std::vector<unsigned int>& offsets,
                  std::vector<typename FnTy::ValTy>& val_vec,
                  galois::DynamicBitSet& bit_set_compute, unsigned vecIndex,
                  size_t start = 0) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string doall_str(syncTypeStr + "SetValVector_" +
                          get_run_identifier(loopName));

    if (parallelize) {
      galois::do_all(galois::iterate(start, start + size),
                     [&](unsigned int n) {
                       unsigned int offset;

                       if (identity_offsets)
                         offset = n;
                       else
                         offset = offsets[n];

                       auto lid = indices[offset];
                       set_wrapper<FnTy, syncType>(lid, val_vec[n - start],
                                                   bit_set_compute, vecIndex);
                     },
#if MORE_COMM_STATS
                     galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                     galois::no_stats());
    } else {
      for (unsigned int n = start; n < start + size; ++n) {
        unsigned int offset;

        if (identity_offsets)
          offset = n;
        else
          offset = offsets[n];

        auto lid = indices[offset];
        set_wrapper<FnTy, syncType>(lid, val_vec[n - start], bit_set_compute,
                                    vecIndex);
      }
    }
  }

  /**
   * GPU wrapper function to reduce multiple nodes at once.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be reduce
   *
   * @param x node id to set
   * @param v
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline bool set_batch_wrapper(unsigned x,
                                std::vector<typename FnTy::ValTy>& v) {
    return FnTy::reduce_batch(x, v.data());
  }

  /**
   * GPU wrapper function to set multiple nodes at once.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be broadcast
   *
   * @param x node id to set
   * @param v
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline bool set_batch_wrapper(unsigned x,
                                std::vector<typename FnTy::ValTy>& v) {
    return FnTy::setVal_batch(x, v.data());
  }

  /**
   * GPU wrapper function to reduce multiple nodes at once. More detailed
   * arguments.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be reduce
   *
   * @param x node id to set
   * @param b
   * @param o
   * @param v
   * @param s
   * @param data_mode
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncReduce>::type* = nullptr>
  inline bool set_batch_wrapper(unsigned x, galois::DynamicBitSet& b,
                                std::vector<unsigned int>& o,
                                std::vector<typename FnTy::ValTy>& v, size_t& s,
                                DataCommMode& data_mode) {
    return FnTy::reduce_batch(x, (uint64_t*)b.get_vec().data(), o.data(),
                              v.data(), s, data_mode);
  }

  /**
   * GPU wrapper function to set multiple nodes at once. More detailed
   * arguments.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   * @tparam SyncType Must be broadcast
   *
   * @param x node id to set
   * @param b
   * @param o
   * @param v
   * @param s
   * @param data_mode
   *
   * @returns true if called on GPU device
   */
  template <typename FnTy, SyncType syncType,
            typename std::enable_if<syncType == syncBroadcast>::type* = nullptr>
  inline bool set_batch_wrapper(unsigned x, galois::DynamicBitSet& b,
                                std::vector<unsigned int>& o,
                                std::vector<typename FnTy::ValTy>& v, size_t& s,
                                DataCommMode& data_mode) {
    return FnTy::setVal_batch(x, (uint64_t*)b.get_vec().data(), o.data(),
                              v.data(), s, data_mode);
  }

  /**
   * Converts LIDs of nodes we are interested in into GIDs.
   *
   * @tparam syncType either reduce or broadcast; only used to name the timer
   *
   * @param loopName name of loop used to name timer
   * @param indices Local ids of nodes that we are interested in
   * @param offsets INPUT/OUTPUT holds offsets into "indices" that we should
   * use; after function completion, holds global ids of nodes we are interested
   * in
   */
  template <SyncType syncType>
  void convert_lid_to_gid(const std::string& loopName,
                          const std::vector<size_t>& indices,
                          std::vector<unsigned int>& offsets) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string doall_str(syncTypeStr + "LID2GID_" +
                          get_run_identifier(loopName));
    galois::do_all(galois::iterate(0ul, offsets.size()),
                   [&](unsigned int n) {
                     offsets[n] =
                         static_cast<uint32_t>(getGID(indices[offsets[n]]));
                   },
#if MORE_COMM_STATS
                   galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                   galois::no_stats());
  }

  /**
   * Converts a vector of GIDs into local ids.
   *
   * @tparam syncType either reduce or broadcast; only used to name the timer
   *
   * @param loopName name of loop used to name timer
   * @param offsets holds GIDs to convert to LIDs
   */
  template <SyncType syncType>
  void convert_gid_to_lid(const std::string& loopName,
                          std::vector<unsigned int>& offsets) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string doall_str(syncTypeStr + "GID2LID_" +
                          get_run_identifier(loopName));

    galois::do_all(galois::iterate(0ul, offsets.size()),
                   [&](unsigned int n) { offsets[n] = getLID(offsets[n]); },
#if MORE_COMM_STATS
                   galois::loopname(get_run_identifier(doall_str).c_str()),
#endif
                   galois::no_stats());
  }

  /**
   * Non-bitset extract that uses serializelazy to copy data over to the
   * buffer. REQUIRES that the ValTy be memory copyable.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam syncFnTy struct that has info on how to do synchronization
   *
   * @param loopName loop name used for timers
   * @param from_id
   * @param indices Vector that contains node ids of nodes that we will
   * potentially send things to
   * @param b OUTPUT: buffer that will be sent over the network; contains data
   * based on set bits in bitset
   */
  template <SyncType syncType, typename SyncFnTy,
            typename std::enable_if<galois::runtime::is_memory_copyable<
                typename SyncFnTy::ValTy>::value>::type* = nullptr>
  void syncExtract(std::string loopName, unsigned from_id,
                   std::vector<size_t>& indices,
                   galois::runtime::SendBuffer& b) {
    uint32_t num = indices.size();
    static std::vector<typename SyncFnTy::ValTy> val_vec; // sometimes wasteful
    static std::vector<unsigned int> offsets;
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string extract_timer_str(syncTypeStr + "Extract_" +
                                  get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textract(extract_timer_str.c_str(),
                                                    GRNAME);
    std::string extract_batch_timer_str(syncTypeStr + "ExtractBatch_" +
                                        get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textractbatch(
        extract_batch_timer_str.c_str(), GRNAME);

    DataCommMode data_mode;

    Textract.start();

    if (num > 0) {
      data_mode = onlyData;
      val_vec.resize(num);

      Textractbatch.start();
      bool batch_succeeded =
          extract_batch_wrapper<SyncFnTy, syncType>(from_id, val_vec);
      Textractbatch.stop();

      if (!batch_succeeded) {
        gSerialize(b, onlyData);
        auto lseq = gSerializeLazySeq(
            b, num, (std::vector<typename SyncFnTy::ValTy>*)nullptr);
        extract_subset<SyncFnTy, decltype(lseq), syncType, true, true>(
            loopName, indices, num, offsets, b, lseq);
      } else {
        gSerialize(b, onlyData, val_vec);
      }
    } else {
      data_mode = noData;
      gSerialize(b, noData);
    }

    Textract.stop();

    std::string metadata_str(syncTypeStr + "MetadataMode_" +
                             std::to_string(data_mode) + "_" +
                             get_run_identifier(loopName));
    galois::runtime::reportStatCond_Single<MORE_DIST_STATS>(GRNAME,
                                                            metadata_str, 1);
  }

  /**
   * Non-bitset extract for when the type of the item being sync'd isn't
   * memory copyable.
   *
   * Extracts all of the data for all nodes in indices and saves it into
   * a send buffer for return.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam syncFnTy struct that has info on how to do synchronization
   *
   * @param loopName loop name used for timers
   * @param from_id
   * @param indices Vector that contains node ids of nodes that we will
   * potentially send things to
   * @param b OUTPUT: buffer that will be sent over the network; contains data
   * based on set bits in bitset
   */
  template <SyncType syncType, typename SyncFnTy,
            typename std::enable_if<!galois::runtime::is_memory_copyable<
                typename SyncFnTy::ValTy>::value>::type* = nullptr>
  void syncExtract(std::string loopName, unsigned from_id,
                   std::vector<size_t>& indices,
                   galois::runtime::SendBuffer& b) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string extract_timer_str(syncTypeStr + "Extract_" +
                                  get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textract(extract_timer_str.c_str(),
                                                    GRNAME);
    std::string extract_batch_timer_str(syncTypeStr + "ExtractBatch_" +
                                        get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textractbatch(
        extract_batch_timer_str.c_str(), GRNAME);

    DataCommMode data_mode;

    uint32_t num = indices.size();
    static std::vector<typename SyncFnTy::ValTy> val_vec;
    static std::vector<unsigned int> dummyVector;

    Textract.start();

    if (num > 0) {
      data_mode = onlyData;
      val_vec.resize(num);

      Textractbatch.start();
      bool batch_succeeded =
          extract_batch_wrapper<SyncFnTy, syncType>(from_id, val_vec);
      Textractbatch.stop();

      if (!batch_succeeded) {
        // get everything (note I pass in "indices" as offsets as it won't
        // even get used anyways)
        extract_subset<SyncFnTy, syncType, true, true>(loopName, indices, num,
                                                       dummyVector, val_vec);
      }

      gSerialize(b, onlyData, val_vec);
    } else {
      data_mode = noData;
      gSerialize(b, noData);
    }

    Textract.stop();

    std::string metadata_str(syncTypeStr + "MetadataMode_" +
                             std::to_string(data_mode) + "_" +
                             get_run_identifier(loopName));
    galois::runtime::reportStatCond_Single<MORE_DIST_STATS>(GRNAME,
                                                            metadata_str, 1);
  }

  /**
   * Reports bytes saved by using the bitset to only selectively load data
   * to send.
   *
   * @tparam SyncFnTy synchronization structure with info needed to synchronize;
   * used for size calculation
   *
   * @param loopName loop name used for timers
   * @param syncTypeStr String used to name timers
   * @param totalToSend Total amount of nodes that are potentially sent (not
   * necessarily all nodees will be sent)
   * @param bitSetCount Number of nodes that will actually be sent
   * @param bitSetComm bitset used to send data
   */
  template <typename SyncFnTy>
  void reportRedundantSize(std::string loopName, std::string syncTypeStr,
                           uint32_t totalToSend, size_t bitSetCount,
                           const galois::DynamicBitSet& bitSetComm) {
    size_t redundant_size =
        (totalToSend - bitSetCount) * sizeof(typename SyncFnTy::ValTy);
    size_t bit_set_size = (bitSetComm.get_vec().size() * sizeof(uint64_t));

    if (redundant_size > bit_set_size) {
      std::string statSavedBytes_str(syncTypeStr + "SavedBytes_" +
                                     get_run_identifier(loopName));

      galois::runtime::reportStatCond_Tsum<MORE_DIST_STATS>(
          GRNAME, statSavedBytes_str, (redundant_size - bit_set_size));
    }
  }

  /**
   * Given data to serialize in val_vec, serialize it into the send buffer
   * depending on the mode of data communication selected for the data.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam VecType type of val_vec, which stores the data to send
   *
   * @param loopName loop name used for timers
   * @param data_mode the way that the data should be communicated
   * @param bit_set_count the number of items we are sending in this message
   * @param indices list of all nodes that we are potentially interested in
   * sending things to
   * @param offsets contains indicies into "indices" that we are interested in
   * @param val_vec contains the data that we are serializing to send
   * @param b the buffer in which to serialize the message we are sending
   * to
   */
  template <SyncType syncType, typename VecType>
  void serializeMessage(std::string loopName, DataCommMode data_mode,
                        size_t bit_set_count, std::vector<size_t>& indices,
                        std::vector<unsigned int>& offsets,
                        galois::DynamicBitSet& bit_set_comm, VecType& val_vec,
                        galois::runtime::SendBuffer& b) {
    if (data_mode == noData) {
      gSerialize(b, data_mode);
    } else if (data_mode == gidsData) {
      offsets.resize(bit_set_count);
      convert_lid_to_gid<syncType>(loopName, indices, offsets);
      val_vec.resize(bit_set_count);
      gSerialize(b, data_mode, bit_set_count, offsets, val_vec);
    } else if (data_mode == offsetsData) {
      offsets.resize(bit_set_count);
      val_vec.resize(bit_set_count);
      gSerialize(b, data_mode, bit_set_count, offsets, val_vec);
    } else if (data_mode == bitsetData) {
      val_vec.resize(bit_set_count);
      gSerialize(b, data_mode, bit_set_count, bit_set_comm, val_vec);
    } else { // onlyData
      gSerialize(b, data_mode, val_vec);
    }
  }

  /**
   * Extracts the data that will be sent to a host in this round of
   * synchronization based on the passed in bitset and saves it to a
   * send buffer.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam syncFnTy struct that has info on how to do synchronization
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   * being used for the extraction
   *
   * @param loopName loop name used for timers
   * @param from_id
   * @param indices Vector that contains node ids of nodes that we will
   * potentially send things to
   * @param b OUTPUT: buffer that will be sent over the network; contains data
   * based on set bits in bitset
   */
  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<!BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  void syncExtract(std::string loopName, unsigned from_id,
                   std::vector<size_t>& indices,
                   galois::runtime::SendBuffer& b) {
    uint32_t num = indices.size();
    static galois::DynamicBitSet bit_set_comm;
    static std::vector<typename SyncFnTy::ValTy> val_vec;
    static std::vector<unsigned int> offsets;

    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string extract_timer_str(syncTypeStr + "Extract_" +
                                  get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textract(extract_timer_str.c_str(),
                                                    GRNAME);
    std::string extract_batch_timer_str(syncTypeStr + "ExtractBatch_" +
                                        get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textractbatch(
        extract_batch_timer_str.c_str(), GRNAME);

    DataCommMode data_mode;

    Textract.start();

    if (num > 0) {
      bit_set_comm.resize(num);
      offsets.resize(num);
      val_vec.resize(num);
      size_t bit_set_count = 0;

      Textractbatch.start();
      bool batch_succeeded = extract_batch_wrapper<SyncFnTy, syncType>(
          from_id, bit_set_comm, offsets, val_vec, bit_set_count, data_mode);
      Textractbatch.stop();

      // GPUs have a batch function they can use; CPUs do not; therefore,
      // CPUS always enter this if block
      if (!batch_succeeded) {
        const galois::DynamicBitSet& bit_set_compute = BitsetFnTy::get();

        get_bitset_and_offsets<SyncFnTy, syncType>(
            loopName, indices, bit_set_compute, bit_set_comm, offsets,
            bit_set_count, data_mode);

        if (data_mode == onlyData) {
          bit_set_count = indices.size();
          extract_subset<SyncFnTy, syncType, true, true>(
              loopName, indices, bit_set_count, offsets, val_vec);
        } else if (data_mode !=
                   noData) { // bitsetData or offsetsData or gidsData
          extract_subset<SyncFnTy, syncType, false, true>(
              loopName, indices, bit_set_count, offsets, val_vec);
        }
      }

      reportRedundantSize<SyncFnTy>(loopName, syncTypeStr, num, bit_set_count,
                                    bit_set_comm);
      serializeMessage<syncType>(loopName, data_mode, bit_set_count, indices,
                                 offsets, bit_set_comm, val_vec, b);
    } else {
      data_mode = noData;
      gSerialize(b, noData);
    }

    Textract.stop();

    std::string metadata_str(syncTypeStr + "MetadataMode_" +
                             std::to_string(data_mode) + "_" +
                             get_run_identifier(loopName));
    galois::runtime::reportStatCond_Single<MORE_DIST_STATS>(GRNAME,
                                                            metadata_str, 1);
  }

  /**
   * Vector bitset variant.
   *
   * Extracts the data that will be sent to a host in this round of
   * synchronization based on the passed in bitset and saves it to a
   * send buffer. Unlike other variants, this will extract an entire
   * vector element by element.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam syncFnTy struct that has info on how to do synchronization
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   * being used for the extraction. MUST BE A VECTOR BITSET
   *
   * @param loopName loop name used for timers
   * @param from_id
   * @param indices Vector that contains node ids of nodes that we will
   * potentially send things to
   * @param b OUTPUT: buffer that will be sent over the network; contains data
   * based on set bits in bitset
   */
  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  void syncExtract(std::string loopName, unsigned from_id,
                   std::vector<size_t>& indices,
                   galois::runtime::SendBuffer& b) {
    uint32_t num = indices.size();
    static galois::DynamicBitSet bit_set_comm;
    static std::vector<typename SyncFnTy::ValTy> val_vec;
    static std::vector<unsigned int> offsets;

    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string extract_timer_str(syncTypeStr + "ExtractVector_" +
                                  get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Textract(extract_timer_str.c_str(),
                                                    GRNAME);

    Textract.start();

    if (num > 0) {
      bit_set_comm.resize(num);
      offsets.resize(num);
      val_vec.resize(num);
    }

    DataCommMode data_mode;
    // loop over all bitsets in the vector of bitsets; each one corresponds to
    // a different index in the vector field we are synchronizing
    for (unsigned i = 0; i < BitsetFnTy::numBitsets(); i++) {
      if (num > 0) {
        bit_set_comm.reset();

        size_t bit_set_count = 0;

        // No GPU support currently
        const galois::DynamicBitSet& bit_set_compute = BitsetFnTy::get(i);

        get_bitset_and_offsets<SyncFnTy, syncType>(
            loopName, indices, bit_set_compute, bit_set_comm, offsets,
            bit_set_count, data_mode);

        // note the extra template argument which specifies that this is a
        // vector extract, i.e. get element i of the vector (i passed in as
        // argument as well)
        if (data_mode == onlyData) {
          // galois::gInfo(id, " node ", i, " has data to send");
          bit_set_count = indices.size();
          extract_subset<SyncFnTy, syncType, true, true, true>(
              loopName, indices, bit_set_count, offsets, val_vec, i);
        } else if (data_mode !=
                   noData) { // bitsetData or offsetsData or gidsData
          // galois::gInfo(id, " node ", i, " has data to send");
          extract_subset<SyncFnTy, syncType, false, true, true>(
              loopName, indices, bit_set_count, offsets, val_vec, i);
        }

        reportRedundantSize<SyncFnTy>(loopName, syncTypeStr, num, bit_set_count,
                                      bit_set_comm);
        serializeMessage<syncType>(loopName, data_mode, bit_set_count, indices,
                                   offsets, bit_set_comm, val_vec, b);
      } else {
        // append noData for however many bitsets there are
        gSerialize(b, noData);
      }
    }

    Textract.stop();

    // FIXME report metadata mode for the different bitsets?
    // std::string metadata_str(syncTypeStr + "_METADATA_MODE" +
    //                         std::to_string(data_mode) +
    //                         get_run_identifier(loopName));
    // galois::runtime::reportStat_Single(GRNAME, metadata_str, 1);
  }

  /**
   * Get data that is going to be sent for synchronization and returns
   * it in a send buffer.
   *
   * @tparam syncType synchronization type
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has information needed to access bitset
   *
   * @param loopName Name to give timer
   * @param x Host to send to
   * @param b OUTPUT: Buffer that will hold data to send
   */
  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<!BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  void get_send_buffer(std::string loopName, unsigned x,
                       galois::runtime::SendBuffer& b) {
    auto& sharedNodes = (syncType == syncReduce) ? mirrorNodes : masterNodes;

    if (BitsetFnTy::is_valid()) {
      syncExtract<syncType, SyncFnTy, BitsetFnTy>(loopName, x, sharedNodes[x],
                                                  b);
    } else {
      syncExtract<syncType, SyncFnTy>(loopName, x, sharedNodes[x], b);
    }

    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string statSendBytes_str(syncTypeStr + "SendBytes_" +
                                  get_run_identifier(loopName));

    galois::runtime::reportStatCond_Tsum<MORE_COMM_STATS>(
        GRNAME, statSendBytes_str, b.size());
  }

  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  void get_send_buffer(std::string loopName, unsigned x,
                       galois::runtime::SendBuffer& b) {
    auto& sharedNodes = (syncType == syncReduce) ? mirrorNodes : masterNodes;

    syncExtract<syncType, SyncFnTy, BitsetFnTy>(loopName, x, sharedNodes[x], b);

    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string statSendBytes_str(syncTypeStr + "SendBytesVector_" +
                                  get_run_identifier(loopName));

    galois::runtime::reportStatCond_Tsum<MORE_COMM_STATS>(
        GRNAME, statSendBytes_str, b.size());
  }

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
  /**
   * Sync using MPI instead of network layer.
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_mpi_send(std::string loopName) {
    static std::vector<galois::runtime::SendBuffer> b;
    static std::vector<MPI_Request> request;
    b.resize(numHosts);
    request.resize(numHosts, MPI_REQUEST_NULL);

    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + h) % numHosts;

      if (nothingToSend(x, syncType, writeLocation, readLocation))
        continue;

      int ready = 0;
      MPI_Test(&request[x], &ready, MPI_STATUS_IGNORE);
      if (!ready) {
        assert(b[x].size() > 0);
        MPI_Wait(&request[x], MPI_STATUS_IGNORE);
      }
      if (b[x].size() > 0) {
        b[x].getVec().clear();
      }

      get_send_buffer<syncType, SyncFnTy, BitsetFnTy>(loopName, x, b[x]);

      MPI_Isend((uint8_t*)b[x].linearData(), b[x].size(), MPI_BYTE, x, 32767,
                MPI_COMM_WORLD, &request[x]);
    }

    if (BitsetFnTy::is_valid()) {
      reset_bitset(syncType, &BitsetFnTy::reset_range);
    }
  }

  /**
   * Sync put using MPI instead of network layer
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_mpi_put(std::string loopName, const MPI_Group& mpi_access_group,
                    const std::vector<MPI_Win>& window) {

    MPI_Win_start(mpi_access_group, 0, window[id]);

    std::vector<galois::runtime::SendBuffer> b(numHosts);
    std::vector<size_t> size(numHosts);
    uint64_t send_buffers_size = 0;

    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + h) % numHosts;

      if (nothingToSend(x, syncType, writeLocation, readLocation))
        continue;

      get_send_buffer<syncType, SyncFnTy, BitsetFnTy>(loopName, x, b[x]);

      size[x] = b[x].size();
      send_buffers_size += size[x];
      MPI_Put((uint8_t*)&size[x], sizeof(size_t), MPI_BYTE, x, 0,
              sizeof(size_t), MPI_BYTE, window[id]);
      MPI_Put((uint8_t*)b[x].linearData(), size[x], MPI_BYTE, x, sizeof(size_t),
              size[x], MPI_BYTE, window[id]);
    }

    auto& net = galois::runtime::getSystemNetworkInterface();
    net.incrementMemUsage(send_buffers_size);

    MPI_Win_complete(window[id]);
    net.decrementMemUsage(send_buffers_size);

    if (BitsetFnTy::is_valid()) {
      reset_bitset(syncType, &BitsetFnTy::reset_range);
    }
  }
#endif

  /**
   * Sends data to all hosts (if there is anything that needs to be sent
   * to that particular host) and adjusts bitset according to sync type.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has information needed to access bitset
   *
   * @param loopName used to name timers created by this sync send
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_net_send(std::string loopName) {
    auto& net = galois::runtime::getSystemNetworkInterface();

    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + h) % numHosts;

      if (nothingToSend(x, syncType, writeLocation, readLocation))
        continue;

      galois::runtime::SendBuffer b;

      get_send_buffer<syncType, SyncFnTy, BitsetFnTy>(loopName, x, b);

      net.sendTagged(x, galois::runtime::evilPhase, b);
    }
    // Will force all messages to be processed before continuing
    net.flush();

    if (BitsetFnTy::is_valid()) {
      reset_bitset(syncType, &BitsetFnTy::reset_range);
    }
  }

  /**
   * Sends data over the network to other hosts based on the provided template
   * arguments.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_send(std::string loopName) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    galois::CondStatTimer<MORE_COMM_STATS> TSendTime(
        (syncTypeStr + "Send_" + get_run_identifier(loopName)).c_str(), GRNAME);

    TSendTime.start();
    sync_net_send<writeLocation, readLocation, syncType, SyncFnTy, BitsetFnTy>(
        loopName);
    TSendTime.stop();
  }

  /**
   * Given the data mode, deserialize the rest of a message in a Receive Buffer.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam VecType type of val_vec, which data will be deserialized into
   *
   * @param loopName used to name timers for statistics
   * @param data_mode data mode with which the original message was sent;
   * determines how to deserialize the rest of the message
   * @param buf buffer which contains the received message to deserialize
   *
   * The rest of the arguments are output arguments (they are passed by
   * reference)
   *
   * @param bit_set_count Var that holds number of bits set (i.e. number of
   * node changed) after deserialization
   * @param offsets holds offsets data after deserialization if data mode is
   * offsets + data
   * @param bit_set_comm holds the bitset representing changed nodes after
   * deserialization of data mode is bitset + data
   * @param buf_start
   * @param retval
   * @param val_vec The data proper will be deserialized into this vector
   */
  template <SyncType syncType, typename VecType>
  void deserializeData(std::string loopName, DataCommMode data_mode,
                       uint32_t num, galois::runtime::RecvBuffer& buf,
                       size_t& bit_set_count,
                       std::vector<unsigned int>& offsets,
                       galois::DynamicBitSet& bit_set_comm, size_t& buf_start,
                       size_t& retval, VecType& val_vec) {
    // get other metadata associated with message if mode isn't OnlyData
    if (data_mode != onlyData) {
      galois::runtime::gDeserialize(buf, bit_set_count);

      if (data_mode == gidsData) {
        galois::runtime::gDeserialize(buf, offsets);
        convert_gid_to_lid<syncType>(loopName, offsets);
      } else if (data_mode == offsetsData) {
        galois::runtime::gDeserialize(buf, offsets);
      } else if (data_mode == bitsetData) {
        bit_set_comm.resize(num);
        galois::runtime::gDeserialize(buf, bit_set_comm);
      } else if (data_mode == dataSplit) {
        galois::runtime::gDeserialize(buf, buf_start);
      } else if (data_mode == dataSplitFirst) {
        galois::runtime::gDeserialize(buf, retval);
      }
    }

    // get data itself
    galois::runtime::gDeserialize(buf, val_vec);
  }

  /**
   * Deserializes messages from other hosts and applies them to update local
   * data based on the provided sync structures.
   *
   * Complement of syncExtract.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param from_id ID of host which the message we are processing was received
   * from
   * @param buf Buffer that contains received message from other host
   * @param loopName used to name timers for statistics
   */
  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<!BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  size_t syncRecvApply(uint32_t from_id, galois::runtime::RecvBuffer& buf,
                       std::string loopName) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string set_timer_str(syncTypeStr + "Set_" +
                              get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Tset(set_timer_str.c_str(), GRNAME);
    std::string set_batch_timer_str(syncTypeStr + "SetBatch_" +
                                    get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Tsetbatch(
        set_batch_timer_str.c_str(), GRNAME);

    static galois::DynamicBitSet bit_set_comm;
    static std::vector<typename SyncFnTy::ValTy> val_vec;
    static std::vector<unsigned int> offsets;

    auto& sharedNodes = (syncType == syncReduce) ? masterNodes : mirrorNodes;
    uint32_t num      = sharedNodes[from_id].size();
    size_t retval     = 0;

    Tset.start();

    if (num > 0) { // only enter if we expect message from that host
      DataCommMode data_mode;
      // 1st deserialize gets data mode
      galois::runtime::gDeserialize(buf, data_mode);

      if (data_mode != noData) {
        size_t bit_set_count = num;
        size_t buf_start     = 0;

        // deserialize the rest of the data in the buffer depending on the data
        // mode; arguments passed in here are mostly output vars
        deserializeData<syncType>(loopName, data_mode, num, buf, bit_set_count,
                                  offsets, bit_set_comm, buf_start, retval,
                                  val_vec);

        // GPU update call
        Tsetbatch.start();
        bool batch_succeeded = set_batch_wrapper<SyncFnTy, syncType>(
            from_id, bit_set_comm, offsets, val_vec, bit_set_count, data_mode);
        Tsetbatch.stop();

        // cpu always enters this block
        if (!batch_succeeded) {
          galois::DynamicBitSet& bit_set_compute = BitsetFnTy::get();

          if (data_mode == bitsetData) {
            size_t bit_set_count2;
            get_offsets_from_bitset<syncType>(loopName, bit_set_comm, offsets,
                                              bit_set_count2);
            assert(bit_set_count == bit_set_count2);
          }

          if (data_mode == onlyData) {
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType, true,
                       true>(loopName, sharedNodes[from_id], bit_set_count,
                             offsets, val_vec, bit_set_compute);
          } else if (data_mode == dataSplit || data_mode == dataSplitFirst) {
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType, true,
                       true>(loopName, sharedNodes[from_id], bit_set_count,
                             offsets, val_vec, bit_set_compute, buf_start);
          } else if (data_mode == gidsData) {
            set_subset<decltype(offsets), SyncFnTy, syncType, true, true>(
                loopName, offsets, bit_set_count, offsets, val_vec,
                bit_set_compute);
          } else { // bitsetData or offsetsData
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType,
                       false, true>(loopName, sharedNodes[from_id],
                                    bit_set_count, offsets, val_vec,
                                    bit_set_compute);
          }
          // TODO: reduce could update the bitset, so it needs to be copied
          // back to the device
        }
      }
    }

    Tset.stop();

    return retval;
  }

  /**
   * VECTOR BITSET VARIANT.
   *
   * Deserializes messages from other hosts and applies them to update local
   * data based on the provided sync structures. Each message will contain
   * a series of messages that must be deserialized (the number of such
   * messages corresponds to the size of the vector that is being synchronized).
   *
   * Complement of syncExtract, vector bitset version.
   *
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   * MUST BE VECTOR BITSET
   *
   * @param from_id ID of host which the message we are processing was received
   * from
   * @param buf Buffer that contains received message from other host
   * @param loopName used to name timers for statistics
   */
  template <
      SyncType syncType, typename SyncFnTy, typename BitsetFnTy,
      typename std::enable_if<BitsetFnTy::is_vector_bitset()>::type* = nullptr>
  size_t syncRecvApply(uint32_t from_id, galois::runtime::RecvBuffer& buf,
                       std::string loopName) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    std::string set_timer_str(syncTypeStr + "SetVector_" +
                              get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> Tset(set_timer_str.c_str(), GRNAME);

    static galois::DynamicBitSet bit_set_comm;
    static std::vector<typename SyncFnTy::ValTy> val_vec;
    static std::vector<unsigned int> offsets;

    auto& sharedNodes = (syncType == syncReduce) ? masterNodes : mirrorNodes;
    uint32_t num      = sharedNodes[from_id].size();
    size_t retval     = 0;

    Tset.start();

    if (num > 0) { // only enter if we expect message from that host
      for (unsigned i = 0; i < BitsetFnTy::numBitsets(); i++) {
        DataCommMode data_mode;
        // 1st deserialize gets data mode
        galois::runtime::gDeserialize(buf, data_mode);

        if (data_mode != noData) {
          size_t bit_set_count = num;
          size_t buf_start     = 0;

          // deserialize the rest of the data in the buffer depending on the
          // data mode; arguments passed in here are mostly output vars
          deserializeData<syncType>(loopName, data_mode, num, buf,
                                    bit_set_count, offsets, bit_set_comm,
                                    buf_start, retval, val_vec);

          galois::DynamicBitSet& bit_set_compute = BitsetFnTy::get(i);

          if (data_mode == bitsetData) {
            size_t bit_set_count2;
            get_offsets_from_bitset<syncType>(loopName, bit_set_comm, offsets,
                                              bit_set_count2);
            assert(bit_set_count == bit_set_count2);
          }

          // Note the extra template argument and i argument which cause
          // execution to deal with a particular element of the vector field
          // we are synchronizing
          if (data_mode == onlyData) {
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType, true,
                       true, true>(loopName, sharedNodes[from_id],
                                   bit_set_count, offsets, val_vec,
                                   bit_set_compute, i);
          } else if (data_mode == dataSplit || data_mode == dataSplitFirst) {
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType, true,
                       true, true>(loopName, sharedNodes[from_id],
                                   bit_set_count, offsets, val_vec,
                                   bit_set_compute, i, buf_start);
          } else if (data_mode == gidsData) {
            set_subset<decltype(offsets), SyncFnTy, syncType, true, true, true>(
                loopName, offsets, bit_set_count, offsets, val_vec,
                bit_set_compute, i);
          } else { // bitsetData or offsetsData
            set_subset<decltype(sharedNodes[from_id]), SyncFnTy, syncType,
                       false, true, true>(loopName, sharedNodes[from_id],
                                          bit_set_count, offsets, val_vec,
                                          bit_set_compute, i);
          }
        }
      }
    }

    Tset.stop();

    return retval;
  }

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
  /**
   * MPI Irecv wrapper for sync
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_mpi_recv_post(std::string loopName,
                          std::vector<MPI_Request>& request,
                          const std::vector<std::vector<uint8_t>>& rb) {
    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + numHosts - h) % numHosts;
      if (nothingToRecv(x, syncType, writeLocation, readLocation))
        continue;

      MPI_Irecv((uint8_t*)rb[x].data(), rb[x].size(), MPI_BYTE, x, 32767,
                MPI_COMM_WORLD, &request[x]);
    }
  }

  /**
   * MPI receive wrapper for sync
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_mpi_recv_wait(std::string loopName,
                          std::vector<MPI_Request>& request,
                          const std::vector<std::vector<uint8_t>>& rb) {
    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + numHosts - h) % numHosts;
      if (nothingToRecv(x, syncType, writeLocation, readLocation))
        continue;

      MPI_Status status;
      MPI_Wait(&request[x], &status);

      int size = 0;
      MPI_Get_count(&status, MPI_BYTE, &size);

      galois::runtime::RecvBuffer rbuf(rb[x].begin(), rb[x].begin() + size);

      syncRecvApply<syncType, SyncFnTy, BitsetFnTy>(x, rbuf, loopName);
    }
  }

  /**
   * MPI get wrapper for sync
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_mpi_get(std::string loopName, const std::vector<MPI_Win>& window,
                    const std::vector<std::vector<uint8_t>>& rb) {
    for (unsigned h = 1; h < numHosts; ++h) {
      unsigned x = (id + numHosts - h) % numHosts;
      if (nothingToRecv(x, syncType, writeLocation, readLocation))
        continue;

      MPI_Win_wait(window[x]);

      size_t size = 0;
      memcpy(&size, rb[x].data(), sizeof(size_t));

      galois::runtime::RecvBuffer rbuf(rb[x].begin() + sizeof(size_t),
                                       rb[x].begin() + sizeof(size_t) + size);

      MPI_Win_post(mpi_identity_groups[x], 0, window[x]);

      syncRecvApply<syncType, SyncFnTy, BitsetFnTy>(x, rbuf, loopName);
    }
  }
#endif

  /**
   * Determines if there is anything to receive from a host and receives/applies
   * the messages.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_net_recv(std::string loopName) {
    auto& net = galois::runtime::getSystemNetworkInterface();

    for (unsigned x = 0; x < numHosts; ++x) {
      if (x == id)
        continue;
      if (nothingToRecv(x, syncType, writeLocation, readLocation))
        continue;

      decltype(net.recieveTagged(galois::runtime::evilPhase, nullptr)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase, nullptr);
      } while (!p);

      syncRecvApply<syncType, SyncFnTy, BitsetFnTy>(p->first, p->second,
                                                    loopName);
    }
    increment_evilPhase();
  }

  /**
   * Receives messages from all other hosts and "applies" the message (reduce
   * or set) based on the sync structure provided.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam syncType either reduce or broadcast
   * @tparam SyncFnTy synchronization structure with info needed to synchronize
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_recv(std::string loopName) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    galois::CondStatTimer<MORE_COMM_STATS> TRecvTime(
        (syncTypeStr + "Recv_" + get_run_identifier(loopName)).c_str(), GRNAME);

    TRecvTime.start();
    sync_net_recv<writeLocation, readLocation, syncType, SyncFnTy, BitsetFnTy>(
        loopName);
    TRecvTime.stop();
  }

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
  /**
   * Nonblocking MPI sync
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_nonblocking_mpi(std::string loopName,
                            bool use_bitset_to_send = true) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    galois::CondStatTimer<MORE_COMM_STATS> TSendTime(
        (syncTypeStr + "Send_" + get_run_identifier(loopName)).c_str(), GRNAME);
    galois::CondStatTimer<MORE_COMM_STATS> TRecvTime(
        (syncTypeStr + "Recv_" + get_run_identifier(loopName)).c_str(), GRNAME);

    static std::vector<std::vector<uint8_t>> rb;
    static std::vector<MPI_Request> request;

    if (rb.size() == 0) { // create the receive buffers
      TRecvTime.start();
      auto& sharedNodes = (syncType == syncReduce) ? masterNodes : mirrorNodes;
      rb.resize(numHosts);
      request.resize(numHosts, MPI_REQUEST_NULL);

      for (unsigned h = 1; h < numHosts; ++h) {
        unsigned x = (id + numHosts - h) % numHosts;
        if (nothingToRecv(x, syncType, writeLocation, readLocation))
          continue;

        size_t size =
            (sharedNodes[x].size() * sizeof(typename SyncFnTy::ValTy));
        size += sizeof(size_t);       // vector size
        size += sizeof(DataCommMode); // data mode

        rb[x].resize(size);
      }
      TRecvTime.stop();
    }

    TRecvTime.start();
    sync_mpi_recv_post<writeLocation, readLocation, syncType, SyncFnTy,
                       BitsetFnTy>(loopName, request, rb);
    TRecvTime.stop();

    TSendTime.start();
    if (use_bitset_to_send) {
      sync_mpi_send<writeLocation, readLocation, syncType, SyncFnTy,
                    BitsetFnTy>(loopName);
    } else {
      sync_mpi_send<writeLocation, readLocation, syncType, SyncFnTy,
                    galois::InvalidBitsetFnTy>(loopName);
    }
    TSendTime.stop();

    TRecvTime.start();
    sync_mpi_recv_wait<writeLocation, readLocation, syncType, SyncFnTy,
                       BitsetFnTy>(loopName, request, rb);
    TRecvTime.stop();
  }

  /**
   * Onesided MPI sync
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            SyncType syncType, typename SyncFnTy, typename BitsetFnTy>
  void sync_onesided_mpi(std::string loopName, bool use_bitset_to_send = true) {
    std::string syncTypeStr = (syncType == syncReduce) ? "Reduce" : "Broadcast";
    galois::CondStatTimer<MORE_COMM_STATS> TSendTime(
        (syncTypeStr + "Send_" + get_run_identifier(loopName)).c_str(), GRNAME);
    galois::CondStatTimer<MORE_COMM_STATS> TRecvTime(
        (syncTypeStr + "Recv_" + get_run_identifier(loopName)).c_str(), GRNAME);

    static std::vector<MPI_Win> window;
    static MPI_Group mpi_access_group;
    static std::vector<std::vector<uint8_t>> rb;

    if (window.size() == 0) { // create the windows
      TRecvTime.start();
      auto& sharedNodes = (syncType == syncReduce) ? masterNodes : mirrorNodes;
      window.resize(numHosts);
      rb.resize(numHosts);

      uint64_t recv_buffers_size = 0;
      for (unsigned x = 0; x < numHosts; ++x) {
        size_t size =
            (sharedNodes[x].size() * sizeof(typename SyncFnTy::ValTy));
        size += sizeof(size_t);       // vector size
        size += sizeof(DataCommMode); // data mode
        size += sizeof(size_t);       // buffer size
        recv_buffers_size += size;

        rb[x].resize(size);

        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Info_set(info, "no_locks", "true");
        MPI_Info_set(info, "same_disp_unit", "true");

        MPI_Win_create(rb[x].data(), size, 1, info, MPI_COMM_WORLD, &window[x]);

        MPI_Info_free(&info);
      }
      auto& net = galois::runtime::getSystemNetworkInterface();
      net.incrementMemUsage(recv_buffers_size);

      for (unsigned h = 1; h < numHosts; ++h) {
        unsigned x = (id + numHosts - h) % numHosts;
        if (nothingToRecv(x, syncType, writeLocation, readLocation))
          continue;
        // exposure group of each window is same as identity group of that
        // window
        MPI_Win_post(mpi_identity_groups[x], 0, window[x]);
      }
      TRecvTime.stop();

      TSendTime.start();
      std::vector<int> access_hosts;
      for (unsigned h = 1; h < numHosts; ++h) {
        unsigned x = (id + h) % numHosts;

        if (nothingToSend(x, syncType, writeLocation, readLocation))
          continue;

        access_hosts.push_back(x);
      }
      MPI_Group world_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      // access group for only one window since only one window is accessed
      MPI_Group_incl(world_group, access_hosts.size(), access_hosts.data(),
                     &mpi_access_group);
      TSendTime.stop();
    }

    TSendTime.start();
    if (use_bitset_to_send) {
      sync_mpi_put<writeLocation, readLocation, syncType, SyncFnTy, BitsetFnTy>(
          loopName, mpi_access_group, window);
    } else {
      sync_mpi_put<writeLocation, readLocation, syncType, SyncFnTy,
                   galois::InvalidBitsetFnTy>(loopName, mpi_access_group,
                                              window);
    }
    TSendTime.stop();

    TRecvTime.start();
    sync_mpi_get<writeLocation, readLocation, syncType, SyncFnTy, BitsetFnTy>(
        loopName, window, rb);
    TRecvTime.stop();
  }
#endif

  /**
   * Does a reduction of data from mirror nodes to master nodes.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            typename ReduceFnTy, typename BitsetFnTy>
  inline void reduce(std::string loopName) {
    std::string timer_str("Reduce_" + get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> TsyncReduce(timer_str.c_str(),
                                                       GRNAME);
    TsyncReduce.start();

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
    switch (bare_mpi) {
    case noBareMPI:
#endif
      sync_send<writeLocation, readLocation, syncReduce, ReduceFnTy,
                BitsetFnTy>(loopName);
      sync_recv<writeLocation, readLocation, syncReduce, ReduceFnTy,
                BitsetFnTy>(loopName);
#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
      break;
    case nonBlockingBareMPI:
      sync_nonblocking_mpi<writeLocation, readLocation, syncReduce, ReduceFnTy,
                           BitsetFnTy>(loopName);
      break;
    case oneSidedBareMPI:
      sync_onesided_mpi<writeLocation, readLocation, syncReduce, ReduceFnTy,
                        BitsetFnTy>(loopName);
      break;
    default:
      GALOIS_DIE("Unsupported bare MPI");
    }
#endif

    TsyncReduce.stop();
  }

  /**
   * Does a broadcast of data from master to mirror nodes.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            typename BroadcastFnTy, typename BitsetFnTy>
  inline void broadcast(std::string loopName) {
    std::string timer_str("Broadcast_" + get_run_identifier(loopName));
    galois::CondStatTimer<MORE_COMM_STATS> TsyncBroadcast(timer_str.c_str(),
                                                          GRNAME);

    TsyncBroadcast.start();

    bool use_bitset = true;

    if (currentBVFlag != nullptr) {
      if (readLocation == readSource &&
          galois::runtime::src_invalid(*currentBVFlag)) {
        use_bitset     = false;
        *currentBVFlag = BITVECTOR_STATUS::NONE_INVALID;
        currentBVFlag  = nullptr;
      } else if (readLocation == readDestination &&
                 galois::runtime::dst_invalid(*currentBVFlag)) {
        use_bitset     = false;
        *currentBVFlag = BITVECTOR_STATUS::NONE_INVALID;
        currentBVFlag  = nullptr;
      } else if (readLocation == readAny &&
                 *currentBVFlag != BITVECTOR_STATUS::NONE_INVALID) {
        // the bitvector flag being non-null means this call came from
        // sync on demand; sync on demand will NEVER use readAny
        // if location is read Any + one of src or dst is invalid
        GALOIS_DIE("readAny + use of bitvector flag without none_invalid "
                   "should never happen");
      }
    }

#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
    switch (bare_mpi) {
    case noBareMPI:
#endif
      if (use_bitset) {
        sync_send<writeLocation, readLocation, syncBroadcast, BroadcastFnTy,
                  BitsetFnTy>(loopName);
      } else {
        sync_send<writeLocation, readLocation, syncBroadcast, BroadcastFnTy,
                  galois::InvalidBitsetFnTy>(loopName);
      }
      sync_recv<writeLocation, readLocation, syncBroadcast, BroadcastFnTy,
                BitsetFnTy>(loopName);
#ifdef __GALOIS_BARE_MPI_COMMUNICATION__
      break;
    case nonBlockingBareMPI:
      sync_nonblocking_mpi<writeLocation, readLocation, syncBroadcast,
                           BroadcastFnTy, BitsetFnTy>(loopName, use_bitset);
      break;
    case oneSidedBareMPI:
      sync_onesided_mpi<writeLocation, readLocation, syncBroadcast,
                        BroadcastFnTy, BitsetFnTy>(loopName, use_bitset);
      break;
    default:
      GALOIS_DIE("Unsupported bare MPI");
    }
#endif

    TsyncBroadcast.stop();
  }

  // OEC - outgoing edge-cut : source of any edge is master
  // IEC - incoming edge-cut : destination of any edge is master
  // CVC - cartesian vertex-cut : if source of an edge is mirror,
  //                              then destination is not, and vice-versa
  // UVC - unconstrained vertex-cut
  // Reduce - mirrors to master
  // Broadcast - master to mirrors

  /**
   * Do sync necessary for write source, read source.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_src_to_src(std::string loopName) {
    // do nothing for OEC
    // reduce and broadcast for IEC, CVC, UVC
    if (transposed || is_vertex_cut()) {
      reduce<writeSource, readSource, ReduceFnTy, BitsetFnTy>(loopName);
      broadcast<writeSource, readSource, BroadcastFnTy, BitsetFnTy>(loopName);
    }
  }

  /**
   * Do sync necessary for write source, read destination.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_src_to_dst(std::string loopName) {
    // only broadcast for OEC
    // only reduce for IEC
    // reduce and broadcast for CVC, UVC
    if (transposed) {
      reduce<writeSource, readDestination, ReduceFnTy, BitsetFnTy>(loopName);
      if (is_vertex_cut()) {
        broadcast<writeSource, readDestination, BroadcastFnTy, BitsetFnTy>(
            loopName);
      }
    } else {
      if (is_vertex_cut()) {
        reduce<writeSource, readDestination, ReduceFnTy, BitsetFnTy>(loopName);
      }
      broadcast<writeSource, readDestination, BroadcastFnTy, BitsetFnTy>(
          loopName);
    }
  }

  /**
   * Do sync necessary for write source, read any.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_src_to_any(std::string loopName) {
    // only broadcast for OEC
    // reduce and broadcast for IEC, CVC, UVC
    if (transposed || is_vertex_cut()) {
      reduce<writeSource, readAny, ReduceFnTy, BitsetFnTy>(loopName);
    }
    broadcast<writeSource, readAny, BroadcastFnTy, BitsetFnTy>(loopName);
  }

  /**
   * Do sync necessary for write dest, read source.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_dst_to_src(std::string loopName) {
    // only reduce for OEC
    // only broadcast for IEC
    // reduce and broadcast for CVC, UVC
    if (transposed) {
      if (is_vertex_cut()) {
        reduce<writeDestination, readSource, ReduceFnTy, BitsetFnTy>(loopName);
      }
      broadcast<writeDestination, readSource, BroadcastFnTy, BitsetFnTy>(
          loopName);
    } else {
      reduce<writeDestination, readSource, ReduceFnTy, BitsetFnTy>(loopName);
      if (is_vertex_cut()) {
        broadcast<writeDestination, readSource, BroadcastFnTy, BitsetFnTy>(
            loopName);
      }
    }
  }

  /**
   * Do sync necessary for write dest, read dest.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_dst_to_dst(std::string loopName) {
    // do nothing for IEC
    // reduce and broadcast for OEC, CVC, UVC
    if (!transposed || is_vertex_cut()) {
      reduce<writeDestination, readDestination, ReduceFnTy, BitsetFnTy>(
          loopName);
      broadcast<writeDestination, readDestination, BroadcastFnTy, BitsetFnTy>(
          loopName);
    }
  }

  /**
   * Do sync necessary for write dest, read any.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_dst_to_any(std::string loopName) {
    // only broadcast for IEC
    // reduce and broadcast for OEC, CVC, UVC
    if (!transposed || is_vertex_cut()) {
      reduce<writeDestination, readAny, ReduceFnTy, BitsetFnTy>(loopName);
    }
    broadcast<writeDestination, readAny, BroadcastFnTy, BitsetFnTy>(loopName);
  }

  /**
   * Do sync necessary for write any, read src.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_any_to_src(std::string loopName) {
    // only reduce for OEC
    // reduce and broadcast for IEC, CVC, UVC
    reduce<writeAny, readSource, ReduceFnTy, BitsetFnTy>(loopName);
    if (transposed || is_vertex_cut()) {
      broadcast<writeAny, readSource, BroadcastFnTy, BitsetFnTy>(loopName);
    }
  }

  /**
   * Do sync necessary for write any, read dst.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_any_to_dst(std::string loopName) {
    // only reduce for IEC
    // reduce and broadcast for OEC, CVC, UVC
    reduce<writeAny, readDestination, ReduceFnTy, BitsetFnTy>(loopName);

    if (!transposed || is_vertex_cut()) {
      broadcast<writeAny, readDestination, BroadcastFnTy, BitsetFnTy>(loopName);
    }
  }

  /**
   * Do sync necessary for write any, read any.
   *
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  inline void sync_any_to_any(std::string loopName) {
    // reduce and broadcast for OEC, IEC, CVC, UVC
    reduce<writeAny, readAny, ReduceFnTy, BitsetFnTy>(loopName);
    broadcast<writeAny, readAny, BroadcastFnTy, BitsetFnTy>(loopName);
  }

public:
  /**
   * Main sync call exposed to the user that calls the correct sync function
   * based on provided template arguments. Must provide information through
   * structures on how to do synchronization/which fields to synchronize.
   *
   * @tparam writeLocation Location data is written (src or dst)
   * @tparam readLocation Location data is read (src or dst)
   * @tparam ReduceFnTy specify how to do reductions
   * @tparam BroadcastFnTy specify how to do broadcasts
   * @tparam BitsetFnTy struct that has info on how to access the bitset
   *
   * @param loopName used to name timers for statistics
   */
  template <WriteLocation writeLocation, ReadLocation readLocation,
            typename ReduceFnTy, typename BroadcastFnTy,
            typename BitsetFnTy = galois::InvalidBitsetFnTy>
  inline void sync(std::string loopName) {
    std::string timer_str("Sync_" + loopName + "_" + get_run_identifier());
    galois::StatTimer Tsync(timer_str.c_str(), GRNAME);

    Tsync.start();

    if (partitionAgnostic) {
      sync_any_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
    } else {
      if (writeLocation == writeSource) {
        if (readLocation == readSource) {
          sync_src_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else if (readLocation == readDestination) {
          sync_src_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else { // readAny
          sync_src_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        }
      } else if (writeLocation == writeDestination) {
        if (readLocation == readSource) {
          sync_dst_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else if (readLocation == readDestination) {
          sync_dst_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else { // readAny
          sync_dst_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        }
      } else { // writeAny
        if (readLocation == readSource) {
          sync_any_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else if (readLocation == readDestination) {
          sync_any_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else { // readAny
          sync_any_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        }
      }
    }

    Tsync.stop();
  }

private:
  /**
   * Generic Sync on demand handler. Should NEVER get to this (hence
   * the galois die).
   */
  template <ReadLocation rl, typename ReduceFnTy, typename BroadcastFnTy,
            typename BitsetFnTy>
  struct SyncOnDemandHandler {
    // note this call function signature is diff. from specialized versions:
    // will cause compile time error if this struct is used (which is what
    // we want)
    void call() { GALOIS_DIE("Invalid read location for sync on demand"); }
  };

  /**
   * Sync on demand handler specialized for read source.
   *
   * @tparam ReduceFnTy specify how to do reductions
   * @tparam BroadcastFnTy specify how to do broadcasts
   * @tparam BitsetFnTy tells program what data needs to be sync'd
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  struct SyncOnDemandHandler<readSource, ReduceFnTy, BroadcastFnTy,
                             BitsetFnTy> {
    /**
     * Based on sync flags, handles syncs for cases when you need to read
     * at source
     *
     * @param g The graph to sync
     * @param fieldFlags the flags structure specifying what needs to be
     * sync'd
     * @param loopName loopname used to name timers
     * @param bvFlag Copy of the bitvector status (valid/invalid at particular
     * locations)
     */
    static inline void call(DistGraph* g,
                            galois::runtime::FieldFlags& fieldFlags,
                            std::string loopName,
                            const BITVECTOR_STATUS& bvFlag) {
      if (fieldFlags.src_to_src() && fieldFlags.dst_to_src()) {
        g->sync_any_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      } else if (fieldFlags.src_to_src()) {
        g->sync_src_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      } else if (fieldFlags.dst_to_src()) {
        g->sync_dst_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      }

      fieldFlags.clear_read_src();
    }
  };

  /**
   * Sync on demand handler specialized for read destination.
   *
   * @tparam ReduceFnTy specify how to do reductions
   * @tparam BroadcastFnTy specify how to do broadcasts
   * @tparam BitsetFnTy tells program what data needs to be sync'd
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  struct SyncOnDemandHandler<readDestination, ReduceFnTy, BroadcastFnTy,
                             BitsetFnTy> {
    /**
     * Based on sync flags, handles syncs for cases when you need to read
     * at destination
     *
     * @param g The graph to sync
     * @param fieldFlags the flags structure specifying what needs to be
     * sync'd
     * @param loopName loopname used to name timers
     * @param bvFlag Copy of the bitvector status (valid/invalid at particular
     * locations)
     */
    static inline void call(DistGraph* g,
                            galois::runtime::FieldFlags& fieldFlags,
                            std::string loopName,
                            const BITVECTOR_STATUS& bvFlag) {
      if (fieldFlags.src_to_dst() && fieldFlags.dst_to_dst()) {
        g->sync_any_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      } else if (fieldFlags.src_to_dst()) {
        g->sync_src_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      } else if (fieldFlags.dst_to_dst()) {
        g->sync_dst_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
      }

      fieldFlags.clear_read_dst();
    }
  };

  /**
   * Sync on demand handler specialized for read any.
   *
   * @tparam ReduceFnTy specify how to do reductions
   * @tparam BroadcastFnTy specify how to do broadcasts
   * @tparam BitsetFnTy tells program what data needs to be sync'd
   */
  template <typename ReduceFnTy, typename BroadcastFnTy, typename BitsetFnTy>
  struct SyncOnDemandHandler<readAny, ReduceFnTy, BroadcastFnTy, BitsetFnTy> {
    /**
     * Based on sync flags, handles syncs for cases when you need to read
     * at both source and destination
     *
     * @param g The graph to sync
     * @param fieldFlags the flags structure specifying what needs to be
     * sync'd
     * @param loopName loopname used to name timers
     * @param bvFlag Copy of the bitvector status (valid/invalid at particular
     * locations)
     */
    static inline void call(DistGraph* g,
                            galois::runtime::FieldFlags& fieldFlags,
                            std::string loopName,
                            const BITVECTOR_STATUS& bvFlag) {
      bool src_write = fieldFlags.src_to_src() || fieldFlags.src_to_dst();
      bool dst_write = fieldFlags.dst_to_src() || fieldFlags.dst_to_dst();

      if (!(src_write && dst_write)) {
        // src or dst write flags aren't set (potentially both are not set),
        // but it's NOT the case that both are set, meaning "any" isn't
        // required in the "from"; can work at granularity of just src
        // write or dst wrte

        if (src_write) {
          if (fieldFlags.src_to_src() && fieldFlags.src_to_dst()) {
            if (bvFlag == BITVECTOR_STATUS::NONE_INVALID) {
              g->sync_src_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else if (galois::runtime::src_invalid(bvFlag)) {
              // src invalid bitset; sync individually so it can be called
              // without bitset
              g->sync_src_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
              g->sync_src_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else if (galois::runtime::dst_invalid(bvFlag)) {
              // dst invalid bitset; sync individually so it can be called
              // without bitset
              g->sync_src_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
              g->sync_src_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else {
              GALOIS_DIE("Invalid bitvector flag setting in sync_on_demand");
            }
          } else if (fieldFlags.src_to_src()) {
            g->sync_src_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          } else { // src to dst is set
            g->sync_src_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          }
        } else if (dst_write) {
          if (fieldFlags.dst_to_src() && fieldFlags.dst_to_dst()) {
            if (bvFlag == BITVECTOR_STATUS::NONE_INVALID) {
              g->sync_dst_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else if (galois::runtime::src_invalid(bvFlag)) {
              g->sync_dst_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
              g->sync_dst_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else if (galois::runtime::dst_invalid(bvFlag)) {
              g->sync_dst_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
              g->sync_dst_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(
                  loopName);
            } else {
              GALOIS_DIE("Invalid bitvector flag setting in sync_on_demand");
            }
          } else if (fieldFlags.dst_to_src()) {
            g->sync_dst_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          } else { // dst to dst is set
            g->sync_dst_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          }
        }

        // note the "no flags are set" case will enter into this block
        // as well, and it is correctly handled by doing nothing since
        // both src/dst_write will be false
      } else {
        // it is the case that both src/dst write flags are set, so "any"
        // is required in the "from"; what remains to be determined is
        // the use of src, dst, or any for the destination of the sync
        bool src_read = fieldFlags.src_to_src() || fieldFlags.dst_to_src();
        bool dst_read = fieldFlags.src_to_dst() || fieldFlags.dst_to_dst();

        if (src_read && dst_read) {
          if (bvFlag == BITVECTOR_STATUS::NONE_INVALID) {
            g->sync_any_to_any<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          } else if (galois::runtime::src_invalid(bvFlag)) {
            g->sync_any_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
            g->sync_any_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          } else if (galois::runtime::dst_invalid(bvFlag)) {
            g->sync_any_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
            g->sync_any_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
          } else {
            GALOIS_DIE("Invalid bitvector flag setting in sync_on_demand");
          }
        } else if (src_read) {
          g->sync_any_to_src<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        } else { // dst_read
          g->sync_any_to_dst<ReduceFnTy, BroadcastFnTy, BitsetFnTy>(loopName);
        }
      }

      fieldFlags.clear_read_src();
      fieldFlags.clear_read_dst();
    }
  };

public:
  /**
   * Given a structure that contains flags signifying what needs to be
   * synchronized, sync_on_demand will synchronize what is necessary based
   * on the read location of the * field.
   *
   * @tparam readLocation Location in which field will need to be read
   * @tparam ReduceFnTy reduce sync structure for the field
   * @tparam BroadcastFnTy broadcast sync structure for the field
   * @tparam BitsetFnTy struct which holds a bitset which can be used
   * to control synchronization at a more fine grain level
   * @param fieldFlags structure for field you are syncing
   * @param loopName Name of loop this sync is for for naming timers
   */
  template <ReadLocation readLocation, typename ReduceFnTy,
            typename BroadcastFnTy,
            typename BitsetFnTy = galois::InvalidBitsetFnTy>
  inline void sync_on_demand(galois::runtime::FieldFlags& fieldFlags,
                             std::string loopName) {
    std::string timer_str("Sync_" + get_run_identifier(loopName));
    galois::StatTimer Tsync(timer_str.c_str(), GRNAME);
    Tsync.start();

    currentBVFlag = &(fieldFlags.bitvectorStatus);

    // call a template-specialized function depending on the read location
    SyncOnDemandHandler<readLocation, ReduceFnTy, BroadcastFnTy,
                        BitsetFnTy>::call(this, fieldFlags, loopName,
                                          *currentBVFlag);

    currentBVFlag = nullptr;

    Tsync.stop();
  }

#ifdef __GALOIS_CHECKPOINT__
public:
  /**
   * Checkpoint the complete structure on the node to disk
   */
  void checkpointSaveNodeData(std::string checkpointFileName = "checkpoint") {
    using namespace boost::archive;
    galois::StatTimer TimerSaveCheckPoint(
        get_run_identifier("TimerSaveCheckpoint").c_str(), GRNAME);

    TimerSaveCheckPoint.start();
    std::string checkpointFileName_local =
        checkpointFileName + "_" + std::to_string(id);

    std::ofstream outputStream(checkpointFileName_local, std::ios::binary);
    if (!outputStream.is_open()) {
      galois::gPrint("ERROR: Could not open ", checkpointFileName_local,
                     " to save checkpoint!!!\n");
    }
    galois::gPrint("[", id,
                   "] Saving local checkpoint to :", checkpointFileName_local,
                   "\n");

    boost::archive::binary_oarchive ar(outputStream, boost::archive::no_header);

    graph.serializeNodeData(ar);

    std::string statSendBytes_str("CheckpointBytesTotal");
    constexpr static const char* const RREGION = "RECOVERY";
    size_t cp_size                             = outputStream.tellp();
    galois::runtime::reportStat_Tsum(RREGION, statSendBytes_str, cp_size);

    outputStream.flush();
    outputStream.close();
    TimerSaveCheckPoint.stop();
  }

  /**
   * Load checkpointed data from disk.
   */
  void checkpointApplyNodeData(std::string checkpointFileName = "checkpoint") {
    using namespace boost::archive;
    galois::StatTimer TimerApplyCheckPoint(
        get_run_identifier("TimerApplyCheckpoint").c_str(), GRNAME);

    TimerApplyCheckPoint.start();
    std::string checkpointFileName_local =
        checkpointFileName + "_" + std::to_string(id);

    std::ifstream inputStream(checkpointFileName_local, std::ios::binary);

    if (!inputStream.is_open()) {
      galois::gPrint("ERROR: Could not open ", checkpointFileName_local,
                     " to read checkpoint!!!\n");
    }
    galois::gPrint("[", id, "] reading local checkpoint from: ",
                   checkpointFileName_local, "\n");

    boost::archive::binary_iarchive ar(inputStream, boost::archive::no_header);

    graph.deSerializeNodeData(ar);

    inputStream.close();
    TimerApplyCheckPoint.stop();
  }
#endif

public:
  /**
   * Converts a local node id into a global node id
   *
   * @param nodeID local node id
   * @returns global node id corresponding to the local one
   */
  inline uint64_t getGID(const uint32_t nodeID) const { return L2G(nodeID); }

  /**
   * Converts a global node id into a local node id
   *
   * @param nodeID global node id
   * @returns local node id corresponding to the global one
   */
  inline uint32_t getLID(const uint64_t nodeID) const { return G2L(nodeID); }

#ifdef __GALOIS_HET_CUDA__
private:
  // Code that handles getting the graph onto the GPU
  template <bool isVoidType,
            typename std::enable_if<isVoidType>::type* = nullptr>
  inline void setMarshalEdge(MarshalGraph& m, const size_t index,
                             const edge_iterator& e) {
    // do nothing
  }

  template <bool isVoidType,
            typename std::enable_if<!isVoidType>::type* = nullptr>
  inline void setMarshalEdge(MarshalGraph& m, const size_t index,
                             const edge_iterator& e) {
    m.edge_data[index] = getEdgeData(e);
  }

public:
  void getMarshalGraph(MarshalGraph& m) {
    m.nnodes = size();
    m.nedges = sizeEdges();
    assert(numOwned > 0);
    m.numOwned          = numOwned;
    m.beginMaster       = beginMaster;
    m.numNodesWithEdges = numNodesWithEdges;
    m.id                = id;
    m.numHosts          = masterNodes.size();
    m.row_start         = (index_type*)calloc(m.nnodes + 1, sizeof(index_type));
    m.edge_dst          = (index_type*)calloc(m.nedges, sizeof(index_type));
    m.node_data         = (index_type*)calloc(m.nnodes, sizeof(node_data_type));

    if (std::is_void<EdgeTy>::value) {
      m.edge_data = NULL;
    } else {
      if (!std::is_same<EdgeTy, edge_data_type>::value) {
        galois::gWarn("Edge data type mismatch between CPU and GPU\n");
      }

      m.edge_data = (edge_data_type*)calloc(m.nedges, sizeof(edge_data_type));
    }

    galois::do_all(
        galois::iterate(graph),
        [&](const typename GraphTy::GraphNode& nodeID) {
          // initialize node_data with localID-to-globalID mapping
          m.node_data[nodeID] = getGID(nodeID);
          m.row_start[nodeID] = *edge_begin(nodeID);
          for (auto e = edge_begin(nodeID); e != edge_end(nodeID); e++) {
            auto edgeID = *e;
            setMarshalEdge<std::is_void<EdgeTy>::value>(m, edgeID, e);
            m.edge_dst[edgeID] = getEdgeDst(e);
          }
        },
        galois::steal());
    m.row_start[m.nnodes] = m.nedges;

    // copy memoization meta-data
    m.num_master_nodes =
        (unsigned int*)calloc(masterNodes.size(), sizeof(unsigned int));
    ;
    m.master_nodes =
        (unsigned int**)calloc(masterNodes.size(), sizeof(unsigned int*));
    ;

    for (uint32_t h = 0; h < masterNodes.size(); ++h) {
      m.num_master_nodes[h] = masterNodes[h].size();

      if (masterNodes[h].size() > 0) {
        m.master_nodes[h] =
            (unsigned int*)calloc(masterNodes[h].size(), sizeof(unsigned int));
        ;
        std::copy(masterNodes[h].begin(), masterNodes[h].end(),
                  m.master_nodes[h]);
      } else {
        m.master_nodes[h] = NULL;
      }
    }

    m.num_mirror_nodes =
        (unsigned int*)calloc(mirrorNodes.size(), sizeof(unsigned int));
    ;
    m.mirror_nodes =
        (unsigned int**)calloc(mirrorNodes.size(), sizeof(unsigned int*));
    ;
    for (uint32_t h = 0; h < mirrorNodes.size(); ++h) {
      m.num_mirror_nodes[h] = mirrorNodes[h].size();

      if (mirrorNodes[h].size() > 0) {
        m.mirror_nodes[h] =
            (unsigned int*)calloc(mirrorNodes[h].size(), sizeof(unsigned int));
        ;
        std::copy(mirrorNodes[h].begin(), mirrorNodes[h].end(),
                  m.mirror_nodes[h]);
      } else {
        m.mirror_nodes[h] = NULL;
      }
    }

    graph.deallocate();
  }
#endif

#ifdef __GALOIS_HET_OPENCL__
  typedef galois::opencl::Graphs::CL_LC_Graph<NodeTy, EdgeTy> CLGraphType;
  typedef typename CLGraphType::NodeDataWrapper CLNodeDataWrapper;
  typedef typename CLGraphType::NodeIterator CLNodeIterator;
  CLGraphType clGraph;

  const cl_mem& device_ptr() const { return clGraph.device_ptr(); }
  CLNodeDataWrapper
  getDataW(GraphNode N,
           galois::MethodFlag mflag = galois::MethodFlag::UNPROTECTED) {
    return clGraph.getDataW(N);
  }
  const CLNodeDataWrapper
  getDataR(GraphNode N,
           galois::MethodFlag mflag = galois::MethodFlag::UNPROTECTED) {
    return clGraph.getDataR(N);
  }
#endif

  /**
   * Set the run number.
   *
   * @param runNum Number to set the run to
   */
  inline void set_num_run(const uint32_t runNum) { num_run = runNum; }

  /**
   * Get the set run number.
   *
   * @returns The set run number saved in the graph
   */
  inline uint32_t get_run_num() const { return num_run; }

  /**
   * Set the round number for use in the run identifier.
   *
   * @param round round number to set to
   */
  inline void set_num_round(const uint32_t round) { num_round = round; }

  /**
   * Get a run identifier using the set run and set round.
   *
   * @returns a string run identifier
   * @deprecated We want to move away from calling this by itself; use ones
   * that take an argument; will be removed once we eliminate all instances
   * of its use from code
   */
  inline std::string get_run_identifier() const {
#if DIST_PER_ROUND_TIMER
    return std::string(std::to_string(num_run) + "_" +
                       std::to_string(num_round));
#else
    return std::string(std::to_string(num_run));
#endif
  }

  /**
   * Get a run identifier using the set run and set round and
   * append to the passed in string.
   *
   * @param loop_name String to append the run identifier
   * @returns String with run identifier appended to passed in loop name
   */
  inline std::string get_run_identifier(std::string loop_name) const {
#if DIST_PER_ROUND_TIMER
    return std::string(std::string(loop_name) + "_" + std::to_string(num_run) +
                       "_" + std::to_string(num_round));
#else
    return std::string(std::string(loop_name) + "_" + std::to_string(num_run));
#endif
  }

  /**
   * Get a run identifier using the set run and set round and
   * append to the passed in string in addition to the number identifier passed
   * in.
   *
   * @param loop_name String to append the run identifier
   * @param alterID another ID with which to add to the timer name.
   *
   * @returns String with run identifier appended to passed in loop name +
   * alterID
   */
  inline std::string get_run_identifier(std::string loop_name,
                                        unsigned alterID) const {
#if DIST_PER_ROUND_TIMER
    return std::string(std::string(loop_name) + "_" + std::to_string(alterID) +
                       "_" + std::to_string(num_run) + "_" +
                       std::to_string(num_round));
#else
    return std::string(std::string(loop_name) + "_" + std::to_string(alterID) +
                       "_" + std::to_string(num_run));
#endif
  }

  /**
   * Write the local LC_CSR graph to the file on a disk.
   *
   * @param localGraphFileName file name to write local graph to.
   */
  void
  save_local_graph_to_file(std::string localGraphFileName = "local_graph") {
    using namespace boost::archive;
    galois::StatTimer dGraphTimerSaveLocalGraph("TimerSaveLocalGraph", GRNAME);
    dGraphTimerSaveLocalGraph.start();

    std::string fileName = localGraphFileName + "_" + std::to_string(id);

    galois::gDebug("[", id, "] inside save_local_graph_to_file \n");

    std::ofstream outputStream(fileName, std::ios::binary);
    if (!outputStream.is_open()) {
      std::cerr << "ERROR: Could not open " << fileName
                << " to save local graph!!!\n";
    }

    boost::archive::binary_oarchive ar(outputStream, boost::archive::no_header);

    // graph topology
    ar << graph;
    ar << numGlobalNodes;
    ar << numGlobalEdges;

    // bool
    ar << transposed;

    // Proxy information
    // TODO: Find better way to serialize vector of vectors in boost
    // serialization
    for (uint32_t i = 0; i < numHosts; ++i) {
      ar << masterNodes[i];
      ar << mirrorNodes[i];
    }

    ar << numOwned;
    ar << beginMaster;
    ar << numNodesWithEdges;
    ar << gid2host;

    // Serialize partitioning scheme specific data structures.
    boostSerializeLocalGraph(ar);

    outputStream.close();
    dGraphTimerSaveLocalGraph.stop();
  }

  /**
   * Read the local LC_CSR graph from the file on a disk.
   *
   * @param localGraphFileName file name to read local graph from.
   */
  void
  read_local_graph_from_file(std::string localGraphFileName = "local_graph") {
    using namespace boost::archive;
    galois::StatTimer dGraphTimerReadLocalGraph("TimerReadLocalGraph", GRNAME);
    dGraphTimerReadLocalGraph.start();

    std::string fileName = localGraphFileName + "_" + std::to_string(id);

    std::ifstream inputStream(fileName, std::ios::binary);
    if (!inputStream.is_open()) {
      std::cerr << "ERROR: Could not open " << fileName
                << " to read local graph!!!\n";
    }

    galois::gPrint("[", id, "] inside read_local_graph_from_file \n");

    boost::archive::binary_iarchive ar(inputStream, boost::archive::no_header);

    // Graph topology
    ar >> graph;
    ar >> numGlobalNodes;
    ar >> numGlobalEdges;

    // bool
    ar >> transposed;

    // Proxy information
    // TODO: Find better way to Deserialize vector of vectors in boost
    // serialization
    for (uint32_t i = 0; i < numHosts; ++i) {
      ar >> masterNodes[i];
      ar >> mirrorNodes[i];
    }

    ar >> numOwned;
    ar >> beginMaster;
    ar >> numNodesWithEdges;
    ar >> gid2host;

    // Serialize partitioning scheme specific data structures.
    boostDeSerializeLocalGraph(ar);

    allNodesRanges.clear();
    masterRanges.clear();
    withEdgeRanges.clear();
    specificRanges.clear();

    // find ranges again
    determineThreadRanges();
    determineThreadRangesMaster();
    determineThreadRangesWithEdges();
    initializeSpecificRanges();

    // Exchange information among hosts
    // send_info_to_host();

    inputStream.close();
    dGraphTimerReadLocalGraph.stop();
  }

  /**
   * Given a sync structure, reset the field specified by the structure
   * to the 0 of the reduction on mirrors.
   *
   * @tparam FnTy structure that specifies how synchronization is to be done
   */
  template <typename FnTy>
  void reset_mirrorField() {
    auto mirrorRanges = getMirrorRanges();
    for (auto r : mirrorRanges) {
      assert(r.first < r.second);

      // GPU call
      bool batch_succeeded = FnTy::reset_batch(r.first, r.second - 1);

      // CPU always enters this block
      if (!batch_succeeded) {
        galois::do_all(
            galois::iterate(r.first, r.second),
            [&](uint32_t lid) {
#ifdef __GALOIS_HET_OPENCL__
              CLNodeDataWrapper d = clGraph.getDataW(lid);
              FnTy::reset(lid, d);
#else
              FnTy::reset(lid, getData(lid));
#endif
            },
            galois::no_stats(),
            galois::loopname(get_run_identifier("RESET:MIRRORS").c_str()));
      }
    }
  }
};

template <typename NodeTy, typename EdgeTy, bool WithInEdges>
constexpr const char* const
    galois::graphs::DistGraph<NodeTy, EdgeTy, WithInEdges>::GRNAME;
} // end namespace graphs
} // end namespace galois

#endif //_GALOIS_DIST_HGRAPH_H
