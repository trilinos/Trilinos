// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2aAlgorithm_decl.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

#include "Kokkos_Sort.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2aAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

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

    auto neighOfINode = graph.getNeighborVertices(rootCandidate);

    LO num_nonaggd_neighbors = 0, num_local_neighbors = 0;
    for (int j = 0; j < neighOfINode.length; j++) {
      LO neigh = neighOfINode(j);
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

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase2aAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
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
void AggregationPhase2aAlgorithm<LO, GO, Node>::
    BuildAggregatesRandom(const ParameterList& params,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          typename AggregationAlgorithmBase<LO, GO, Node>::AggStatType& aggStat,
                          LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
  const int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
  bool matchMLbehavior           = params.get<bool>("aggregation: match ML phase2a");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner    = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors        = aggregates.GetGraphColors();
  const LO numColors = aggregates.GetGraphNumColors();

  auto lclLWGraph = graph;

  LO numLocalNodes      = numRows;
  LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

  const double aggFactor = 0.5;
  double factor          = static_cast<double>(numLocalAggregated) / (numLocalNodes + 1);
  factor                 = pow(factor, aggFactor);

  // LBV on Sept 12, 2019: this looks a little heavy handed,
  // I'm not sure a view is needed to perform atomic updates.
  // If we can avoid this and use a simple LO that would be
  // simpler for later maintenance.
  Kokkos::View<LO, device_type> numLocalAggregates("numLocalAggregates");
  typename Kokkos::View<LO, device_type>::HostMirror h_numLocalAggregates =
      Kokkos::create_mirror_view(numLocalAggregates);
  h_numLocalAggregates() = aggregates.GetNumAggregates();
  Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);

  // Now we create new aggregates using root nodes in all colors other than the first color,
  // as the first color was already exhausted in Phase 1.
  for (int color = 2; color < numColors + 1; ++color) {
    LO tmpNumNonAggregatedNodes = 0;
    Kokkos::parallel_reduce(
        "Aggregation Phase 2a: loop over each individual color",
        Kokkos::RangePolicy<execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO rootCandidate, LO& lNumNonAggregatedNodes) {
          if (aggStat(rootCandidate) == READY &&
              colors(rootCandidate) == color) {
            LO numNeighbors = 0;
            LO aggSize      = 0;
            if (matchMLbehavior) {
              aggSize += 1;
              numNeighbors += 1;
            }

            auto neighbors = lclLWGraph.getNeighborVertices(rootCandidate);

            // Loop over neighbors to count how many nodes could join
            // the new aggregate

            for (int j = 0; j < neighbors.length; ++j) {
              LO neigh = neighbors(j);
              if (neigh != rootCandidate) {
                if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                    (aggStat(neigh) == READY) &&
                    (aggSize < maxNodesPerAggregate)) {
                  ++aggSize;
                }
                ++numNeighbors;
              }
            }

            // If a sufficient number of nodes can join the new aggregate
            // then we actually create the aggregate.
            if (aggSize > minNodesPerAggregate &&
                (aggSize > factor * numNeighbors)) {
              // aggregates.SetIsRoot(rootCandidate);
              LO aggIndex = Kokkos::
                  atomic_fetch_add(&numLocalAggregates(), 1);

              LO numAggregated = 0;

              if (matchMLbehavior) {
                // Add the root.
                aggStat(rootCandidate)         = AGGREGATED;
                vertex2AggId(rootCandidate, 0) = aggIndex;
                procWinner(rootCandidate, 0)   = myRank;
                ++numAggregated;
                --lNumNonAggregatedNodes;
              }

              for (int neighIdx = 0; neighIdx < neighbors.length; ++neighIdx) {
                LO neigh = neighbors(neighIdx);
                if (neigh != rootCandidate) {
                  if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                      (aggStat(neigh) == READY) &&
                      (numAggregated < aggSize)) {
                    aggStat(neigh)         = AGGREGATED;
                    vertex2AggId(neigh, 0) = aggIndex;
                    procWinner(neigh, 0)   = myRank;

                    ++numAggregated;
                    --lNumNonAggregatedNodes;
                  }
                }
              }
            }
          }
        },
        tmpNumNonAggregatedNodes);
    numNonAggregatedNodes += tmpNumNonAggregatedNodes;
  }

  // update aggregate object
  Kokkos::deep_copy(h_numLocalAggregates, numLocalAggregates);
  aggregates.SetNumAggregates(h_numLocalAggregates());
}  // BuildAggregatesRandom

template <class LO, class GO, class Node>
void AggregationPhase2aAlgorithm<LO, GO, Node>::
    BuildAggregatesDeterministic(const ParameterList& params,
                                 const LWGraph_kokkos& graph,
                                 Aggregates& aggregates,
                                 typename AggregationAlgorithmBase<LO, GO, Node>::AggStatType& aggStat,
                                 LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
  const int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner    = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors        = aggregates.GetGraphColors();
  const LO numColors = aggregates.GetGraphNumColors();

  auto lclLWGraph = graph;

  LO numLocalNodes      = procWinner.size();
  LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

  const double aggFactor = 0.5;
  double factor          = as<double>(numLocalAggregated) / (numLocalNodes + 1);
  factor                 = pow(factor, aggFactor);

  Kokkos::View<LO, device_type> numLocalAggregates("numLocalAggregates");
  typename Kokkos::View<LO, device_type>::HostMirror h_numLocalAggregates =
      Kokkos::create_mirror_view(numLocalAggregates);
  h_numLocalAggregates() = aggregates.GetNumAggregates();
  Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);

  // Now we create new aggregates using root nodes in all colors other than the first color,
  // as the first color was already exhausted in Phase 1.
  //
  // In the deterministic version, exactly the same set of aggregates will be created
  // (as the nondeterministic version)
  // because no vertex V can be a neighbor of two vertices of the same color, so two root
  // candidates can't fight over V
  //
  // But, the precise values in vertex2AggId need to match exactly, so just sort the new
  // roots of each color before assigning aggregate IDs

  // numNonAggregatedNodes is the best available upper bound for the number of aggregates
  // which may be created in this phase, so use it for the size of newRoots
  Kokkos::View<LO*, device_type> newRoots("New root LIDs", numNonAggregatedNodes);
  Kokkos::View<LO, device_type> numNewRoots("Number of new aggregates of current color");
  auto h_numNewRoots = Kokkos::create_mirror_view(numNewRoots);
  for (int color = 1; color < numColors + 1; ++color) {
    h_numNewRoots() = 0;
    Kokkos::deep_copy(numNewRoots, h_numNewRoots);
    Kokkos::parallel_for(
        "Aggregation Phase 2a: determining new roots of current color",
        Kokkos::RangePolicy<execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO rootCandidate) {
          if (aggStat(rootCandidate) == READY &&
              colors(rootCandidate) == color) {
            LO aggSize     = 0;
            auto neighbors = lclLWGraph.getNeighborVertices(rootCandidate);
            // Loop over neighbors to count how many nodes could join
            // the new aggregate
            LO numNeighbors = 0;
            for (int j = 0; j < neighbors.length; ++j) {
              LO neigh = neighbors(j);
              if (neigh != rootCandidate) {
                if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                    aggStat(neigh) == READY &&
                    aggSize < maxNodesPerAggregate) {
                  ++aggSize;
                }
                ++numNeighbors;
              }
            }
            // If a sufficient number of nodes can join the new aggregate
            // then we mark rootCandidate as a future root.
            if (aggSize > minNodesPerAggregate && aggSize > factor * numNeighbors) {
              LO newRootIndex        = Kokkos::atomic_fetch_add(&numNewRoots(), 1);
              newRoots(newRootIndex) = rootCandidate;
            }
          }
        });
    Kokkos::deep_copy(h_numNewRoots, numNewRoots);

    if (h_numNewRoots() > 0) {
      // sort the new root indices
      Kokkos::sort(newRoots, 0, h_numNewRoots());
      // now, loop over all new roots again and actually create the aggregates
      LO tmpNumNonAggregatedNodes = 0;
      // First, just find the set of color vertices which will become aggregate roots
      Kokkos::parallel_reduce(
          "Aggregation Phase 2a: create new aggregates",
          Kokkos::RangePolicy<execution_space>(0, h_numNewRoots()),
          KOKKOS_LAMBDA(const LO newRootIndex, LO& lNumNonAggregatedNodes) {
            LO root        = newRoots(newRootIndex);
            LO newAggID    = numLocalAggregates() + newRootIndex;
            auto neighbors = lclLWGraph.getNeighborVertices(root);
            // Loop over neighbors and add them to new aggregate
            aggStat(root)         = AGGREGATED;
            vertex2AggId(root, 0) = newAggID;
            LO aggSize            = 1;
            for (int j = 0; j < neighbors.length; ++j) {
              LO neigh = neighbors(j);
              if (neigh != root) {
                if (lclLWGraph.isLocalNeighborVertex(neigh) &&
                    aggStat(neigh) == READY &&
                    aggSize < maxNodesPerAggregate) {
                  aggStat(neigh)         = AGGREGATED;
                  vertex2AggId(neigh, 0) = newAggID;
                  procWinner(neigh, 0)   = myRank;
                  aggSize++;
                }
              }
            }
            lNumNonAggregatedNodes -= aggSize;
          },
          tmpNumNonAggregatedNodes);
      numNonAggregatedNodes += tmpNumNonAggregatedNodes;
      h_numLocalAggregates() += h_numNewRoots();
      Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);
    }
  }
  aggregates.SetNumAggregates(h_numLocalAggregates());
}

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE2AALGORITHM_DEF_HPP_ */
