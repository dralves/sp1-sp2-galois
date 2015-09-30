/** Kruskal MST -*- C++ -*-
 * @file
 * @section License
 *
 * This file is part of Galois.  Galoisis a gramework to exploit
 * amorphous data-parallelism in irregular programs.
 *
 * Galois is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Galois is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Galois.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @section Copyright
 *
 * Copyright (C) 2015, The University of Texas at Austin. All rights
 * reserved.
 *
 * @section Description
 *
 * Kruskal MST.
 *
 * @author <ahassaan@ices.utexas.edu>
 */

#ifndef KRUSKALSTRICTOBIM_H
#define KRUSKALSTRICTOBIM_H

#include "Galois/Graphs/Graph.h"
#include "Galois/Runtime/LevelExecutor.h"
#include "Galois/WorkList/WorkListExperimental.h"

#include "Kruskal.h"
#include "KruskalParallel.h"


namespace kruskal {


class KruskalStrictOBIM: public Kruskal {
  protected:

  typedef Galois::Graph::FirstGraph<void,void,false> Graph;
  typedef Graph::GraphNode Lockable;
  typedef std::vector<Lockable> VecLocks;

  virtual const std::string getVersion () const { return "Parallel Kruskal using Strict OBIM"; }
  
  struct FindLoopSpec {

    Graph& graph;
    VecLocks& locks;
    VecRep& repVec;
    Accumulator& findIter;


    FindLoopSpec (
        Graph& graph,
        VecLocks& locks,
        VecRep& repVec,
        Accumulator& findIter)
      :
        graph (graph),
        locks (locks),
        repVec (repVec),
        findIter (findIter)

    {}


    template <typename C>
    void operator () (const Edge& e, C& ctx) {
      int repSrc = kruskal::findPCiter_int (e.src, repVec);
      int repDst = kruskal::findPCiter_int (e.dst, repVec);
      // int repSrc = kruskal::getRep_int (e.src, repVec);
      // int repDst = kruskal::getRep_int (e.dst, repVec);
      

      if (repSrc != repDst) {
        graph.getData (locks[repSrc]);
        graph.getData (locks[repDst]);
      }

      findIter += 1;
    }
  };


  struct LinkUpLoopSpec {

    static const unsigned CHUNK_SIZE = 64;

    VecRep& repVec;
    Accumulator& mstSum;
    Accumulator& linkUpIter;

    LinkUpLoopSpec (
        VecRep& repVec,
        Accumulator& mstSum,
        Accumulator& linkUpIter)
      :
        repVec (repVec),
        mstSum (mstSum),
        linkUpIter (linkUpIter)

    {}


    template <typename C>
    void operator () (const Edge& e, C& ctx) {
      int repSrc = kruskal::findPCiter_int (e.src, repVec);
      int repDst = kruskal::findPCiter_int (e.dst, repVec);

      // int repSrc = kruskal::getRep_int (e.src, repVec);
      // int repDst = kruskal::getRep_int (e.dst, repVec);

      if (repSrc != repDst) {
        unionByRank_int (repSrc, repDst, repVec);
        linkUpIter += 1;
        mstSum += e.weight;
      }

    }
  };

  struct GetWeight: public std::unary_function<Edge,Weight_ty> {
    Weight_ty operator () (const Edge& e) const {
      return e.weight;
    }
  };

  virtual void runMST (const size_t numNodes, VecEdge& edges,
      size_t& mstWeight, size_t& totalIter) {

    Graph graph;
    VecLocks locks;
    locks.resize(numNodes);
    Galois::do_all(boost::counting_iterator<size_t>(0), boost::counting_iterator<size_t>(numNodes), [&](size_t i) {
      locks[i] = graph.createNode();
      graph.addNode(locks[i]);
    });
      
    VecRep repVec (numNodes, -1);
    Accumulator findIter;
    Accumulator linkUpIter;
    Accumulator mstSum;

    FindLoopSpec findLoop (graph, locks, repVec, findIter);
    LinkUpLoopSpec linkUpLoop (repVec, mstSum, linkUpIter);

    Galois::TimeAccumulator runningTime;

    runningTime.start ();
    typedef Galois::WorkList::OrderedByIntegerMetric<GetWeight, Galois::WorkList::ChunkedFIFO<128>>::with_barrier<true>::type WL;
    Galois::for_each(edges.begin(), edges.end(),
        [&](const Edge& e, Galois::UserContext<Edge>& ctx) {
          findLoop(e, ctx);
          linkUpLoop(e, ctx);
    }, Galois::wl<WL>());
    //}, Galois::wl<Galois::WorkList::StableIterator<>>());
    //Galois::Runtime::for_each_ordered_level (
    //    Galois::Runtime::makeStandardRange (edges.begin (), edges.end ()),
    //    GetWeight (), std::less<Weight_ty> (), findLoop, linkUpLoop);

    runningTime.stop ();

    mstWeight = mstSum.reduce ();
    totalIter = findIter.reduce ();

    std::cout << "Number of FindLoop iterations = " << findIter.reduce () << std::endl;
    std::cout << "Number of LinkUpLoop iterations = " << linkUpIter.reduce () << std::endl;

    std::cout << "MST running time without initialization/destruction: " << runningTime.get () << std::endl;
  }
};



}// end namespace kruskal




#endif //  KRUSKAL_LEVEL_EXEC_H
