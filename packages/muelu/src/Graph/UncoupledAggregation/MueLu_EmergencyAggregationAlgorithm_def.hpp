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
 * MueLu_EmergencyAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_EMERGENCYAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_EMERGENCYAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_EmergencyAggregationAlgorithm.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void EmergencyAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    // vertex ids for output
    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    const LO  nRows  = graph.GetNodeNumVertices();
    const int myRank = graph.GetComm()->getRank();

    int              aggIndex = -1;
    size_t           aggSize  =  0;
    std::vector<int> aggList(graph.getNodeMaxNumRowEntries());

    LO nLocalAggregates = aggregates.GetNumAggregates();
    for (LO iNode = 0; iNode < nRows; iNode++) {
      if (aggStat[iNode] == AGGREGATED || aggStat[iNode] == IGNORED)
        continue;

      aggSize = 0;
      aggregates.SetIsRoot(iNode);

      aggList[aggSize++] = iNode;
      aggIndex = nLocalAggregates++;

      ArrayView<const LO> neighOfINode = graph.getNeighborVertices(iNode);

      for (LO j = 0; j < neighOfINode.size(); j++) {
        LO neigh = neighOfINode[j];

        if (neigh != iNode && graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY)
          aggList[aggSize++] = neigh;
      }

      // finalize aggregate
      for (size_t k = 0; k < aggSize; k++) {
        aggStat     [aggList[k]] = AGGREGATED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner  [aggList[k]] = myRank;
      }

      numNonAggregatedNodes -= aggSize;
    }

    aggregates.SetNumAggregates(nLocalAggregates);
  }

} // end namespace

#endif /* MUELU_EMERGENCYAGGREGATIONALGORITHM_DEF_HPP_ */
