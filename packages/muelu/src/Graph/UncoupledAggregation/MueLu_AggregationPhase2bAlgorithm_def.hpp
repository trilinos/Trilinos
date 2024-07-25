// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2bAlgorithm_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2bAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");
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

      auto neighOfINode = graph.getNeighborVertices(i);

      for (int j = 0; j < neighOfINode.length; j++) {
        LO neigh = neighOfINode(j);

        // We don't check (neigh != i), as it is covered by checking (aggStat[neigh] == AGGREGATED)
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == AGGREGATED)
          aggWeight[vertex2AggId[neigh]] += connectWeight[neigh];
      }

      int bestScore   = -100000;
      int bestAggId   = -1;
      int bestConnect = -1;

      for (int j = 0; j < neighOfINode.length; j++) {
        LO neigh  = neighOfINode(j);
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

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2bAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const ParameterList& params,
                    const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                    LO& numNonAggregatedNodes) const {
  if (params.get<bool>("aggregation: deterministic")) {
    Monitor m(*this, "BuildAggregatesDeterministic");
    BuildAggregatesDeterministic(params, graph, aggregates, aggStat, numNonAggregatedNodes);
  } else {
    Monitor m(*this, "BuildAggregatesRandom");
    BuildAggregatesRandom(params, graph, aggregates, aggStat, numNonAggregatedNodes);
  }

}  // BuildAggregates

template <class LO, class GO, class Node>
void AggregationPhase2bAlgorithm<LO, GO, Node>::
    BuildAggregatesRandom(const ParameterList& params,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          typename AggregationAlgorithmBase<LO, GO, Node>::AggStatType& aggStat,
                          LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId           = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner             = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors                 = aggregates.GetGraphColors();
  const LO numColors          = aggregates.GetGraphNumColors();
  const LO numLocalAggregates = aggregates.GetNumAggregates();

  auto lclLWGraph = graph;

  const LO defaultConnectWeight = 100;
  const LO penaltyConnectWeight = 10;

  Kokkos::View<LO*, device_type> aggWeight(Kokkos::ViewAllocateWithoutInitializing("aggWeight"), numLocalAggregates);  // This gets re-initialized at the start of each "color" loop
  Kokkos::View<LO*, device_type> connectWeight(Kokkos::ViewAllocateWithoutInitializing("connectWeight"), numRows);
  Kokkos::View<LO*, device_type> aggPenalties("aggPenalties", numLocalAggregates);  // This gets initialized to zero here

  Kokkos::deep_copy(connectWeight, defaultConnectWeight);

  // taw: by running the aggregation routine more than once there is a chance that also
  // non-aggregated nodes with a node distance of two are added to existing aggregates.
  // Assuming that the aggregate size is 3 in each direction running the algorithm only twice
  // should be sufficient.
  // lbv: If the prior phase of aggregation where run without specifying an aggregate size,
  // the distance 2 coloring and phase 1 aggregation actually guarantee that only one iteration
  // is needed to reach distance 2 neighbors.
  int maxIters             = 2;
  int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
  if (maxNodesPerAggregate == std::numeric_limits<int>::max()) {
    maxIters = 1;
  }
  for (int iter = 0; iter < maxIters; ++iter) {
    for (LO color = 1; color <= numColors; ++color) {
      Kokkos::deep_copy(aggWeight, 0);

      // the reduce counts how many nodes are aggregated by this phase,
      // which will then be subtracted from numNonAggregatedNodes
      LO numAggregated = 0;
      Kokkos::parallel_reduce(
          "Aggregation Phase 2b: aggregates expansion",
          Kokkos::RangePolicy<execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LO i, LO& tmpNumAggregated) {
            if (aggStat(i) != READY || colors(i) != color)
              return;

            auto neighOfINode = lclLWGraph.getNeighborVertices(i);
            for (int j = 0; j < neighOfINode.length; j++) {
              LO neigh = neighOfINode(j);

              // We don't check (neigh != i), as it is covered by checking
              // (aggStat[neigh] == AGGREGATED)
              if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                  aggStat(neigh) == AGGREGATED)
                Kokkos::atomic_add(&aggWeight(vertex2AggId(neigh, 0)),
                                   connectWeight(neigh));
            }

            int bestScore   = -100000;
            int bestAggId   = -1;
            int bestConnect = -1;

            for (int j = 0; j < neighOfINode.length; j++) {
              LO neigh = neighOfINode(j);

              if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                  aggStat(neigh) == AGGREGATED) {
                auto aggId = vertex2AggId(neigh, 0);
                int score  = aggWeight(aggId) - aggPenalties(aggId);

                if (score > bestScore) {
                  bestAggId   = aggId;
                  bestScore   = score;
                  bestConnect = connectWeight(neigh);

                } else if (aggId == bestAggId &&
                           connectWeight(neigh) > bestConnect) {
                  bestConnect = connectWeight(neigh);
                }
              }
            }
            if (bestScore >= 0) {
              aggStat(i)         = AGGREGATED;
              vertex2AggId(i, 0) = bestAggId;
              procWinner(i, 0)   = myRank;

              Kokkos::atomic_add(&aggPenalties(bestAggId), 1);
              connectWeight(i) = bestConnect - penaltyConnectWeight;
              tmpNumAggregated++;
            }
          },
          numAggregated);  // parallel_for
      numNonAggregatedNodes -= numAggregated;
    }
  }  // loop over maxIters

}  // BuildAggregatesRandom

template <class LO, class GO, class Node>
void AggregationPhase2bAlgorithm<LO, GO, Node>::
    BuildAggregatesDeterministic(const ParameterList& params,
                                 const LWGraph_kokkos& graph,
                                 Aggregates& aggregates,
                                 typename AggregationAlgorithmBase<LO, GO, Node>::AggStatType& aggStat,
                                 LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId     = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner       = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors           = aggregates.GetGraphColors();
  const LO numColors    = aggregates.GetGraphNumColors();
  LO numLocalAggregates = aggregates.GetNumAggregates();

  auto lclLWGraph = graph;

  const int defaultConnectWeight = 100;
  const int penaltyConnectWeight = 10;

  Kokkos::View<int*, device_type> connectWeight(Kokkos::ViewAllocateWithoutInitializing("connectWeight"), numRows);
  Kokkos::View<int*, device_type> aggWeight(Kokkos::ViewAllocateWithoutInitializing("aggWeight"), numLocalAggregates);  // This gets re-initialized at the start of each "color" loop
  Kokkos::View<int*, device_type> aggPenaltyUpdates("aggPenaltyUpdates", numLocalAggregates);
  Kokkos::View<int*, device_type> aggPenalties("aggPenalties", numLocalAggregates);

  Kokkos::deep_copy(connectWeight, defaultConnectWeight);

  // We do this cycle twice.
  // I don't know why, but ML does it too
  // taw: by running the aggregation routine more than once there is a chance that also
  // non-aggregated nodes with a node distance of two are added to existing aggregates.
  // Assuming that the aggregate size is 3 in each direction running the algorithm only twice
  // should be sufficient.
  int maxIters             = 2;
  int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
  if (maxNodesPerAggregate == std::numeric_limits<int>::max()) {
    maxIters = 1;
  }
  for (int iter = 0; iter < maxIters; ++iter) {
    for (LO color = 1; color <= numColors; color++) {
      Kokkos::deep_copy(aggWeight, 0);

      // the reduce counts how many nodes are aggregated by this phase,
      // which will then be subtracted from numNonAggregatedNodes
      LO numAggregated = 0;
      Kokkos::parallel_for(
          "Aggregation Phase 2b: updating agg weights",
          Kokkos::RangePolicy<execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LO i) {
            if (aggStat(i) != READY || colors(i) != color)
              return;
            auto neighOfINode = lclLWGraph.getNeighborVertices(i);
            for (int j = 0; j < neighOfINode.length; j++) {
              LO neigh = neighOfINode(j);
              // We don't check (neigh != i), as it is covered by checking
              // (aggStat[neigh] == AGGREGATED)
              if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                  aggStat(neigh) == AGGREGATED)
                Kokkos::atomic_add(&aggWeight(vertex2AggId(neigh, 0)),
                                   connectWeight(neigh));
            }
          });

      Kokkos::parallel_reduce(
          "Aggregation Phase 2b: aggregates expansion",
          Kokkos::RangePolicy<execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LO i, LO& tmpNumAggregated) {
            if (aggStat(i) != READY || colors(i) != color)
              return;
            int bestScore   = -100000;
            int bestAggId   = -1;
            int bestConnect = -1;

            auto neighOfINode = lclLWGraph.getNeighborVertices(i);
            for (int j = 0; j < neighOfINode.length; j++) {
              LO neigh = neighOfINode(j);

              if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                  aggStat(neigh) == AGGREGATED) {
                auto aggId = vertex2AggId(neigh, 0);
                int score  = aggWeight(aggId) - aggPenalties(aggId);

                if (score > bestScore) {
                  bestAggId   = aggId;
                  bestScore   = score;
                  bestConnect = connectWeight(neigh);

                } else if (aggId == bestAggId &&
                           connectWeight(neigh) > bestConnect) {
                  bestConnect = connectWeight(neigh);
                }
              }
            }
            if (bestScore >= 0) {
              aggStat(i)         = AGGREGATED;
              vertex2AggId(i, 0) = bestAggId;
              procWinner(i, 0)   = myRank;

              Kokkos::atomic_add(&aggPenaltyUpdates(bestAggId), 1);
              connectWeight(i) = bestConnect - penaltyConnectWeight;
              tmpNumAggregated++;
            }
          },
          numAggregated);  // parallel_reduce

      Kokkos::parallel_for(
          "Aggregation Phase 2b: updating agg penalties",
          Kokkos::RangePolicy<execution_space>(0, numLocalAggregates),
          KOKKOS_LAMBDA(const LO agg) {
            aggPenalties(agg) += aggPenaltyUpdates(agg);
            aggPenaltyUpdates(agg) = 0;
          });
      numNonAggregatedNodes -= numAggregated;
    }
  }  // loop over k
}  // BuildAggregatesDeterministic

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE2BALGORITHM_DEF_HPP_ */
