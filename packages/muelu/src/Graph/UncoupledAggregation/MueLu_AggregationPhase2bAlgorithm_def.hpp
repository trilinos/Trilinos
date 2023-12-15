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
#ifndef MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2bAlgorithm_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2bAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregates");
  bool matchMLbehavior = params.get<bool>("aggregation: match ML phase2b");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  LO numLocalAggregates = aggregates.GetNumAggregates();

  const int defaultConnectWeight = 100;
  const int penaltyConnectWeight = 10;

  std::vector<int> aggWeight(numLocalAggregates, 0);
  std::vector<int> connectWeight(numRows, defaultConnectWeight);
  std::vector<int> aggPenalties(numRows, 0);

  // We do this cycle twice.
  // I don't know why, but ML does it too
  // taw: by running the aggregation routine more than once there is a chance that also
  // non-aggregated nodes with a node distance of two are added to existing aggregates.
  // Assuming that the aggregate size is 3 in each direction running the algorithm only twice
  // should be sufficient.
  for (int k = 0; k < 2; k++) {
    for (LO i = 0; i < numRows; i++) {
      if (aggStat[i] != READY)
        continue;

      ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(i);

      for (int j = 0; j < neighOfINode.size(); j++) {
        LO neigh = neighOfINode[j];

        // We don't check (neigh != i), as it is covered by checking (aggStat[neigh] == AGGREGATED)
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == AGGREGATED)
          aggWeight[vertex2AggId[neigh]] += connectWeight[neigh];
      }

      int bestScore   = -100000;
      int bestAggId   = -1;
      int bestConnect = -1;

      for (int j = 0; j < neighOfINode.size(); j++) {
        LO neigh  = neighOfINode[j];
        int aggId = vertex2AggId[neigh];

        // Note: The third condition is only relevant if the ML matching is enabled
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == AGGREGATED && (!matchMLbehavior || aggWeight[aggId] != 0)) {
          int score = aggWeight[aggId] - aggPenalties[aggId];

          if (score > bestScore) {
            bestAggId   = aggId;
            bestScore   = score;
            bestConnect = connectWeight[neigh];

          } else if (aggId == bestAggId && connectWeight[neigh] > bestConnect) {
            bestConnect = connectWeight[neigh];
          }

          // Reset the weights for the next loop
          aggWeight[aggId] = 0;
        }
      }

      if (bestScore >= 0) {
        aggStat[i]      = AGGREGATED;
        vertex2AggId[i] = bestAggId;
        procWinner[i]   = myRank;

        numNonAggregatedNodes--;

        aggPenalties[bestAggId]++;
        connectWeight[i] = bestConnect - penaltyConnectWeight;
      }
    }
  }
}

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_ */
