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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_OnePtAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_InterfaceAggregationAlgorithm.hpp"

//#include "MueLu_Graph.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::InterfaceAggregationAlgorithm(RCP<const FactoryBase> const &graphFact)
{
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(Teuchos::ParameterList const & params, GraphBase const & graph, Aggregates & aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregates");

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal numLocalAggregates = aggregates.GetNumAggregates();    // number of local aggregates on current proc

  // main loop over all local rows of graph(A)
  for(int iNode1 = 0; iNode1 < nRows; ++iNode1) {

    if (aggStat[iNode1] == INTERFACE) {

      aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'agg'
      int aggIndex = numLocalAggregates;
      std::vector<int> aggList;
      aggList.push_back(iNode1);
      ArrayView<const LO> neighOfINode = graph.getNeighborVertices(iNode1);

      for(int j = 0; j < neighOfINode.size(); ++j) {
        LO neigh = neighOfINode[j];
        if(neigh != iNode1 && graph.isLocalNeighborVertex(neigh)) {
          if(aggStat[neigh] != AGGREGATED && aggStat[neigh] != INTERFACE &&
             aggStat[neigh] != IGNORED) {
            aggList.push_back(neigh);
          }
        }
      }

      for (size_t k = 0; k < aggList.size(); k++) {
        aggStat[aggList[k]]      = AGGREGATED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
      }
      ++numLocalAggregates;
      numNonAggregatedNodes -= aggList.size();
    }

  } // end for

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

} // end namespace


#endif /* MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_ */
