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
#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase1Algorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    AggOptions::Ordering ordering    = params.get<AggOptions::Ordering>("Ordering");
    LO MaxNeighAlreadySelected       = params.get<LO>                  ("MaxNeighAlreadySelected");
    LO MinNodesPerAggregate          = params.get<LO>                  ("MinNodesPerAggregate");
    LO MaxNodesPerAggregate          = params.get<LO>                  ("MaxNodesPerAggregate");

    TEUCHOS_TEST_FOR_EXCEPTION(MaxNodesPerAggregate < MinNodesPerAggregate, Exceptions::RuntimeError, "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: MinNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

    if (ordering != NATURAL && ordering != RANDOM && ordering != GRAPH)
      throw Exceptions::RuntimeError("UncoupledAggregation::BuildAggregates : bad aggregation ordering option");

    const LO  nRows  = graph.GetNodeNumVertices();
    const int myRank = graph.GetComm()->getRank();

    // vertex ids for output
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    // some internal variables
    LO nLocalAggregates = aggregates.GetNumAggregates();    // number of local aggregates on current proc
    std::queue<LO> graph_ordering_inodes; // inodes for graph ordering

    ArrayRCP<LO> randomVector;
    if (ordering == RANDOM) {
      randomVector = arcp<LO>(nRows);
      for (LO i = 0; i < nRows; i++)
        randomVector[i] = i;
      RandomReorder(randomVector);
    }

    int              aggIndex = -1;
    size_t           aggSize  =  0;
    std::vector<int> aggList(graph.getNodeMaxNumRowEntries());

    // Main loop over all local rows of graph(A)
    for (LO iNode2 = 0; iNode2 < nRows; iNode2++) {
      // Step 1: pick the next node to aggregate
      LO iNode1 = 0;
      if      (ordering == NATURAL) iNode1 = iNode2;
      else if (ordering == RANDOM)  iNode1 = randomVector[iNode2];
      else if (ordering == GRAPH) {

        if (graph_ordering_inodes.size() == 0) {
          // There are no nodes for graph ordering scheme,
          // add exactly one ready node for graph ordering aggregates
          for (LO jnode = 0; jnode < nRows; jnode++)
            if (aggStat[jnode] == NodeStats::READY) {
              graph_ordering_inodes.push(jnode);
              break;
            }
        }
        if (graph_ordering_inodes.size() == 0) {
          // There are no more ready nodes, end the phase
          break;
        }
        iNode1 = graph_ordering_inodes.front();   // take next node from graph ordering queue
        graph_ordering_inodes.pop();              // delete this node in list
      }

      if (aggStat[iNode1] == NodeStats::READY) {
        // Step 2: build tentative aggregate
        aggSize = 0;
        aggList[aggSize++] = iNode1;

        ArrayView<const LO> neighOfINode = graph.getNeighborVertices(iNode1);

        LO numAggregatedNeighbours = 0;

        // NOTE: if neighOfINode.size() < MinNodesPerAggregate, we could skip this loop,
        // but only for NATURAL and RANDOM (for GRAPH we still need the list of local neighbors)
        for (LO j = 0; j < neighOfINode.size(); j++) {
          LO neigh = neighOfINode[j];

          if (neigh != iNode1 && graph.isLocalNeighborVertex(neigh)) {

            if (aggStat[neigh] == NodeStats::READY || aggStat[neigh] == NodeStats::NOTSEL) {
              // Add neighbor node to tentative aggregate
              // but only if aggregate size is not exceeding maximum size
              // NOTE: We do not exit the loop over all neighbours since we have still
              //       to count all aggregated neighbour nodes for the aggregation criteria
              // NOTE: We check here for the maximum aggregation size. If we would do it below
              //       with all the other check too big aggregates would not be accepted at all.
              if (aggSize < as<size_t>(MaxNodesPerAggregate))
                aggList[aggSize++] = neigh;

            } else {
              numAggregatedNeighbours++;
            }
          }
        }

        // Step 3: check if tentative aggregate is acceptable
        if ((numAggregatedNeighbours <= MaxNeighAlreadySelected) &&   // too many connections to other aggregates
            (as<LO>(aggSize)         >= MinNodesPerAggregate)) {      // too few nodes in the tentative aggregate
          // Accept new aggregate
          // iNode1 becomes the root of the newly formed aggregate
          aggregates.SetIsRoot(iNode1);
          aggIndex = nLocalAggregates++;

          for (size_t k = 0; k < aggSize; k++) {
            aggStat     [aggList[k]] = NodeStats::AGGREGATED;
            vertex2AggId[aggList[k]] = aggIndex;
            procWinner  [aggList[k]] = myRank;

            if (ordering == GRAPH) {
              Teuchos::ArrayView<const LO> neighOfJNode = graph.getNeighborVertices(aggList[k]);
              for (int j = 0; j < neighOfJNode.size(); j++) {
                LO neigh = neighOfJNode[j];

                if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == NodeStats::READY)
                  graph_ordering_inodes.push(neigh);
              }
            }
          }

          numNonAggregatedNodes -= aggSize;

        } else {
          // Aggregate is not accepted
          aggStat[iNode1] = NodeStats::NOTSEL;

          if (ordering == GRAPH) {
            // Even though the aggregate around iNode1 is not perfect, we want to try
            // the neighbor nodes of iNode1
            for (int j = 0; j < neighOfINode.size(); j++) {
              LO neigh = neighOfINode[j];

              if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == NodeStats::READY)
                graph_ordering_inodes.push(neigh);
            }
          }
        }
      }
    }

    // update aggregate object
    aggregates.SetNumAggregates(nLocalAggregates);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(ArrayRCP<LO> list) const {
    //TODO: replace int
    int n = list.size();
    for(int i = 0; i < n-1; i++) {
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const {
    return min + Teuchos::as<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

} // end namespace


#endif /* MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_ */
