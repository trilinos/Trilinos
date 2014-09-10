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
 * MueLu_MaxLinkAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_

#include <algorithm>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_MaxLinkAggregationAlgorithm.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void MaxLinkAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    LO MaxNodesPerAggregate = params.get<LO>("aggregation: max agg size");

    // vertex ids for output
    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    ArrayRCP<LO> aggSizes     = aggregates.ComputeAggregateSizes(); // contains number of nodes in aggregate with given aggId

    const LO  nRows  = graph.GetNodeNumVertices();
    const int myRank = graph.GetComm()->getRank();

    size_t           aggSize = 0;
    std::vector<int> aggList(graph.getNodeMaxNumRowEntries());

    //bool recomputeAggregateSizes=false; // variable not used TODO remove it

    for (LO iNode = 0; iNode < nRows; iNode++) {
      if (aggStat[iNode] == AGGREGATED || aggStat[iNode] == IGNORED)
        continue;

      ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);

      aggSize = 0;
      for (LO j = 0; j < neighOfINode.size(); j++) {
        LO neigh = neighOfINode[j];

        // NOTE: we don't need the check (neigh != iNode), as we work only
        // if aggStat[neigh] == AGGREGATED, which we know is different from aggStat[iNode]
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == AGGREGATED)
          aggList[aggSize++] = vertex2AggId[neigh];
      }

      // Ideally, we would have a _fast_ hash table here.
      // But for the absense of that, sorting works just fine.
      std::sort(aggList.begin(), aggList.begin() + aggSize);

      // terminator
      aggList[aggSize] = -1;

      // Find an aggregate id with most connections to
      LO maxNumConnections =  0, curNumConnections = 0;
      LO selectedAggregate = -1;
      for (size_t i = 0; i < aggSize; i++) {
        curNumConnections++;
        if (aggList[i+1] != aggList[i]) {
          if (curNumConnections > maxNumConnections &&         // only select aggregate if it has more connections
              aggSizes[aggList[i]] < MaxNodesPerAggregate) {   // and if it is not too big (i.e. can have one more node)
            maxNumConnections = curNumConnections;
            selectedAggregate = aggList[i];
            //recomputeAggregateSizes=true;
          }
          curNumConnections = 0;
        }
      }

      // Add node iNode to aggregate
      if (selectedAggregate != -1) {
        aggStat[iNode]      = AGGREGATED;
        vertex2AggId[iNode] = selectedAggregate;
        procWinner[iNode]   = myRank;

        numNonAggregatedNodes--;
      }
    }
  }

} // end namespace


#endif /* MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_ */
