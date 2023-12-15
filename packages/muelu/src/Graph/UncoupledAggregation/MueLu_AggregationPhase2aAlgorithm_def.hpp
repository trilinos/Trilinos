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
#ifndef MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2aAlgorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2aAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregates");

  int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
  int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
  bool matchMLbehavior     = params.get<bool>("aggregation: match ML phase2a");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  LO numLocalAggregates = aggregates.GetNumAggregates();

  LO numLocalNodes      = procWinner.size();
  LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

  const double aggFactor = params.get<double>("aggregation: phase2a agg factor");
  double factor;

  if (matchMLbehavior) {
    // Note: ML uses global counts to set the factor
    // Passing  # of nonaggregated nodes and # of nodes via aggStat
    GO in_data[2] = {(GO)numNonAggregatedNodes, (GO)aggStat.size()};
    GO out_data[2];
    Teuchos::reduceAll(*graph.GetComm(), Teuchos::REDUCE_SUM, 2, in_data, out_data);
    GO phase_one_aggregated = out_data[1] - out_data[0];
    factor                  = as<double>(phase_one_aggregated) / (out_data[1] + 1);

    LO agg_stat_unaggregated = 0;
    LO agg_stat_aggregated   = 0;
    LO agg_stat_bdry         = 0;
    for (LO i = 0; i < (LO)aggStat.size(); i++) {
      if (aggStat[i] == AGGREGATED)
        agg_stat_aggregated++;
      else if (aggStat[i] == BOUNDARY)
        agg_stat_bdry++;
      else
        agg_stat_unaggregated++;
    }

    // NOTE: ML always uses 3 as minNodesPerAggregate
    minNodesPerAggregate = 3;

  } else {
    // MueLu defaults to using local counts to set the factor
    factor = as<double>(numLocalAggregated) / (numLocalNodes + 1);
  }

  // Now apply aggFactor
  factor = pow(factor, aggFactor);

  int aggIndex   = -1;
  size_t aggSize = 0;
  std::vector<int> aggList(graph.getLocalMaxNumRowEntries());

  for (LO rootCandidate = 0; rootCandidate < numRows; rootCandidate++) {
    if (aggStat[rootCandidate] != READY) {
      continue;
    }

    LO numNeighbors = 0;
    aggSize         = 0;
    if (matchMLbehavior) {
      aggList[aggSize++] = rootCandidate;
      numNeighbors++;
    }

    ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(rootCandidate);

    LO num_nonaggd_neighbors = 0, num_local_neighbors = 0;
    for (int j = 0; j < neighOfINode.size(); j++) {
      LO neigh = neighOfINode[j];
      if (graph.isLocalNeighborVertex(neigh))
        num_local_neighbors++;

      if (neigh != rootCandidate) {
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY) {
          // If aggregate size does not exceed max size, add node to the tentative aggregate
          // NOTE: We do not exit the loop over all neighbours since we have still
          //       to count all aggregated neighbour nodes for the aggregation criteria
          // NOTE: We check here for the maximum aggregation size. If we would do it below
          //       with all the other check too big aggregates would not be accepted at all.
          if (aggSize < as<size_t>(maxNodesPerAggregate))
            aggList[aggSize++] = neigh;
          num_nonaggd_neighbors++;
        }

        numNeighbors++;
      }
    }

    bool accept_aggregate;
    if (matchMLbehavior) {
      // ML does this calculation slightly differently than MueLu does by default, specifically it
      // uses the *local* number of neigbors, regardless of what they are.
      // NOTE: ML does zero compression here.  Not sure if it matters
      // NOTE: ML uses a hardcoded value 3 instead of minNodesPerAggregate.  This has been set above
      LO rowi_N = num_local_neighbors;
      num_nonaggd_neighbors++;  // ML counts the node itself as a nonaggd_neighbor
      accept_aggregate = (rowi_N > as<LO>(minNodesPerAggregate)) && (num_nonaggd_neighbors > (factor * rowi_N));
    } else {
      accept_aggregate = (aggSize > as<size_t>(minNodesPerAggregate)) && (aggSize > factor * numNeighbors);
    }

    if (accept_aggregate) {
      // Accept new aggregate
      // rootCandidate becomes the root of the newly formed aggregate
      aggregates.SetIsRoot(rootCandidate);
      aggIndex = numLocalAggregates++;

      for (size_t k = 0; k < aggSize; k++) {
        aggStat[aggList[k]]      = AGGREGATED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
      }

      numNonAggregatedNodes -= aggSize;
    }
  }

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_ */
