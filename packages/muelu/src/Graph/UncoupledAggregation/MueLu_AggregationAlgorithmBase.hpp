// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_AggregationAlgorithmBase.hpp
 *
 *  Created on: Sep 17, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_AGGREGATIONALGORITHMBASE_HPP_
#define MUELU_AGGREGATIONALGORITHMBASE_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

//#include "MueLu_Graph_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_AggOptions.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

using namespace AggOptions; // necessary

/* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggregated. */
namespace NodeStats {
enum NodeState {
  READY   = 1,   /* indicates that a node is available to be */
  /* selected as a root node of an aggregate  */

  NOTSEL  = 2,   /* indicates that a node has been rejected  */
  /* as a root node. This could perhaps be    */
  /* because if this node had been selected a */
  /* small aggregate would have resulted.     */

  AGGREGATED = 3,   /* indicates that a node has been assigned  */
  /* to an aggregate.                         */

  ONEPT    = 4,  /* indicates that a node shall be preserved over all multigrid levels as 1 point aggregate */
  SMALLAGG = 5,   /* indicates that a node shall be aggregated separately from standard nodes with small aggregates (only neighbour nodes which are also marked with the SMALLAGG flag) */
  BOUNDARY = 6     // node is a Dirichlet node and should never be aggregated
};
} // namespace NodeStats



// TODO: dangerous: same definition as in CheapAggregationAlgorithm
class Aggregate {
public:
  int length;                   // current size of aggregate
  int maxLength;                // max size of aggregate
  int index;                    // local aggregate id
  std::vector<int> list;  // list of node ids in aggregate
};

/*!
     @class pure virtual base class for all aggregation algorithms.
     @brief Base class for MueLu aggregation algorithms

     @ingroup MueLuBaseClasses
 */
template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class AggregationAlgorithmBase
: public BaseClass
  {
#undef MUELU_AGGREGATIONALGORITHMBASE_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"
  public:

  //! @name Constructors/Destructors
  //@{

  //! Destructor.
  virtual ~AggregationAlgorithmBase() {}

  //@}

  //! @name Build routines
  //@{

  //! BuildAggregates routine.
  virtual LocalOrdinal BuildAggregates(Teuchos::ParameterList const & params, GraphBase const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat) const = 0;
  //@}

  //! @name Build routines
  //@{

  //! BuildAggregates routine.
  virtual void PrintAggregationInformation(const std::string phase, GraphBase const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat) const {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();
    const LocalOrdinal nRows = graph.GetNodeNumVertices();
    const LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates();

    if(IsPrint(Statistics1)) {
      LO localAggregated  = 0;
      GO globalAggregated = 0;
      GO globalNRows    = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] == NodeStats::AGGREGATED) localAggregated++;
      sumAll(comm, (GO)localAggregated, globalAggregated);
      sumAll(comm, (GO)nRows, globalNRows);
      GetOStream(Statistics1, 0) << "Aggregation (UC): " << phase << " Nodes aggregated = " << globalAggregated << " out of " << globalNRows << " nodes" << std::endl;
      GO nAggregatesGlobal = 0;
      sumAll(comm, (GO)nLocalAggregates, nAggregatesGlobal);
      GetOStream(Statistics1, 0) << "Aggregation (UC): " << phase << " Total aggregates = " << nAggregatesGlobal << std::endl;
    }
    if(IsPrint(Warnings0)) {
      GO localNotAggregated  = 0;
      GO globalNotAggregated = 0;
      for(LO i=0; i<nRows; ++i)
        if(aggStat[i] != NodeStats::AGGREGATED ) localNotAggregated++;
        sumAll(comm, (GO)localNotAggregated, globalNotAggregated);
        if(globalNotAggregated > 0)
          GetOStream(Warnings0,0) << "Aggregation (UC): " << phase << " (WARNING) " << globalNotAggregated << " unaggregated nodes left" << std::endl;
    }

  }

  //@}

  private:

  }; // class AggregationAlgorithmBase

} // namespace MueLu

#define MUELU_AGGREGATIONALGORITHMBASE_SHORT
#endif /* MUELU_AGGREGATIONALGORITHMBASE_HPP_ */
