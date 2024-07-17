// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE3ALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE3ALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase3Algorithm_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big. Otherwise, make a new aggregate
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase3Algorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

  bool makeNonAdjAggs    = false;
  bool error_on_isolated = false;
  if (params.isParameter("aggregation: error on nodes with no on-rank neighbors"))
    error_on_isolated = params.get<bool>("aggregation: error on nodes with no on-rank neighbors");
  if (params.isParameter("aggregation: phase3 avoid singletons"))
    makeNonAdjAggs = params.get<bool>("aggregation: phase3 avoid singletons");

  size_t numSingletons = 0;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  LO numLocalAggregates = aggregates.GetNumAggregates();

  for (LO i = 0; i < numRows; i++) {
    if (aggStat[i] == AGGREGATED || aggStat[i] == IGNORED)
      continue;

    auto neighOfINode = graph.getNeighborVertices(i);

    // We don't want a singleton. So lets see if there is an unaggregated
    // neighbor that we can also put with this point.
    bool isNewAggregate    = false;
    bool failedToAggregate = true;
    for (int j = 0; j < neighOfINode.length; j++) {
      LO neigh = neighOfINode(j);

      if (neigh != i && graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY) {
        isNewAggregate = true;

        aggStat[neigh]      = AGGREGATED;
        vertex2AggId[neigh] = numLocalAggregates;
        procWinner[neigh]   = myRank;

        numNonAggregatedNodes--;
      }
    }

    if (isNewAggregate) {
      // Create new aggregate (not singleton)
      aggStat[i]    = AGGREGATED;
      procWinner[i] = myRank;
      numNonAggregatedNodes--;
      aggregates.SetIsRoot(i);
      vertex2AggId[i] = numLocalAggregates++;

      failedToAggregate = false;
    } else {
      // We do not want a singleton, but there are no non-aggregated
      // neighbors. Lets see if we can connect to any other aggregates
      // NOTE: This is very similar to phase 2b, but simplier: we stop with
      // the first found aggregate
      int j = 0;
      for (; j < neighOfINode.length; j++) {
        LO neigh = neighOfINode(j);

        // We don't check (neigh != rootCandidate), as it is covered by checking (aggStat[neigh] == AGGREGATED)
        if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == AGGREGATED)
          break;
      }

      if (j < neighOfINode.length) {
        // Assign to an adjacent aggregate
        vertex2AggId[i] = vertex2AggId[neighOfINode(j)];
        numNonAggregatedNodes--;
        failedToAggregate = false;
      }
    }

    if (failedToAggregate && makeNonAdjAggs) {
      //  it we are still didn't find an aggregate home for i (i.e., we have
      //  a potential singleton), we are desperate. Basically, we seek to
      //  group i with any other local point to form an aggregate (even if
      //  it is not a neighbor of i. Either we find a vertex that is already
      //  aggregated or not aggregated.
      //    1) if found vertex is aggregated, then assign i to this aggregate
      //    2) if found vertex is not aggregated, create new aggregate

      for (LO ii = 0; ii < numRows; ii++) {  // look for anyone else
        if ((ii != i) && (aggStat[ii] != IGNORED)) {
          failedToAggregate = false;       // found someone so start
          aggStat[i]        = AGGREGATED;  // marking i as aggregated
          procWinner[i]     = myRank;

          if (aggStat[ii] == AGGREGATED)
            vertex2AggId[i] = vertex2AggId[ii];
          else {
            vertex2AggId[i]  = numLocalAggregates;
            vertex2AggId[ii] = numLocalAggregates;
            aggStat[ii]      = AGGREGATED;
            procWinner[ii]   = myRank;
            numNonAggregatedNodes--;  // acounts for ii now being aggregated
            aggregates.SetIsRoot(i);
            numLocalAggregates++;
          }
          numNonAggregatedNodes--;  // accounts for i now being aggregated
          break;
        }  // if ( (ii != i) && (aggStat[ii] != IGNORED ...
      }    // for (LO ii = 0; ...
    }
    if (failedToAggregate) {
      if (error_on_isolated) {
        // Error on this isolated node, as the user has requested
        std::ostringstream oss;
        oss << "MueLu::AggregationPhase3Algorithm::BuildAggregatesNonKokkos: MueLu has detected a non-Dirichlet node that has no on-rank neighbors and is terminating (by user request). " << std::endl;
        oss << "If this error is being generated at level 0, this is due to an initial partitioning problem in your matrix." << std::endl;
        oss << "If this error is being generated at any other level, try turning on repartitioning, which may fix this problem." << std::endl;
        throw Exceptions::RuntimeError(oss.str());
      } else {
        // Create new aggregate (singleton)
        //          this->GetOStream(Warnings1) << "Found singleton: " << i << std::endl;
        numSingletons++;

        aggregates.SetIsRoot(i);
        vertex2AggId[i] = numLocalAggregates++;
        numNonAggregatedNodes--;
      }
    }

    // One way or another, the node is aggregated (possibly into a singleton)
    aggStat[i]    = AGGREGATED;
    procWinner[i] = myRank;

  }  // loop over numRows

  if (numSingletons > 0)
    this->GetOStream(Runtime0) << "  WARNING Rank " << myRank << " singletons :" << numSingletons << " (phase)" << std::endl;

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big. Otherwise, make a new aggregate
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase3Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const ParameterList& params,
                    const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                    LO& numNonAggregatedNodes) const {
  // So far we only have the non-deterministic version of the algorithm...
  if (params.get<bool>("aggregation: deterministic")) {
    Monitor m(*this, "BuildAggregatesDeterministic");
    BuildAggregatesRandom(params, graph, aggregates, aggStat, numNonAggregatedNodes);
  } else {
    Monitor m(*this, "BuildAggregatesRandom");
    BuildAggregatesRandom(params, graph, aggregates, aggStat, numNonAggregatedNodes);
  }
}

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big. Otherwise, make a new aggregate
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase3Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesRandom(const ParameterList& params,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                          LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  bool error_on_isolated = params.get<bool>("aggregation: error on nodes with no on-rank neighbors");
  bool makeNonAdjAggs    = params.get<bool>("aggregation: phase3 avoid singletons");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner    = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors        = aggregates.GetGraphColors();
  const LO numColors = aggregates.GetGraphNumColors();

  auto lclLWGraph = graph;

  Kokkos::View<LO, device_type> numAggregates("numAggregates");
  Kokkos::deep_copy(numAggregates, aggregates.GetNumAggregates());

  Kokkos::View<unsigned*, device_type> aggStatOld(Kokkos::ViewAllocateWithoutInitializing("Initial aggregation status"), aggStat.extent(0));
  Kokkos::deep_copy(aggStatOld, aggStat);
  Kokkos::View<LO, device_type> numNonAggregated("numNonAggregated");
  Kokkos::deep_copy(numNonAggregated, numNonAggregatedNodes);
  for (int color = 1; color < numColors + 1; ++color) {
    Kokkos::parallel_for(
        "Aggregation Phase 3: aggregates clean-up",
        Kokkos::RangePolicy<execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO nodeIdx) {
          // Check if node has already been treated?
          if ((colors(nodeIdx) != color) ||
              (aggStatOld(nodeIdx) == AGGREGATED) ||
              (aggStatOld(nodeIdx) == IGNORED)) {
            return;
          }

          // Grab node neighbors
          auto neighbors = lclLWGraph.getNeighborVertices(nodeIdx);
          LO neighIdx;

          // We don't want a singleton.
          // So lets see if any neighbors can be used to form a new aggregate?
          bool isNewAggregate = false;
          for (int neigh = 0; neigh < neighbors.length; ++neigh) {
            neighIdx = neighbors(neigh);

            if ((neighIdx != nodeIdx) &&
                lclLWGraph.isLocalNeighborVertex(neighIdx) &&
                (aggStatOld(neighIdx) == READY)) {
              isNewAggregate = true;
              break;
            }
          }

          // We can form a new non singleton aggregate!
          if (isNewAggregate) {
            // If this is the aggregate root
            // we need to process the nodes in the aggregate
            const LO aggId           = Kokkos::atomic_fetch_add(&numAggregates(), 1);
            aggStat(nodeIdx)         = AGGREGATED;
            procWinner(nodeIdx, 0)   = myRank;
            vertex2AggId(nodeIdx, 0) = aggId;
            // aggregates.SetIsRoot(nodeIdx);
            Kokkos::atomic_decrement(&numNonAggregated());
            for (int neigh = 0; neigh < neighbors.length; ++neigh) {
              neighIdx = neighbors(neigh);
              if ((neighIdx != nodeIdx) &&
                  lclLWGraph.isLocalNeighborVertex(neighIdx) &&
                  (aggStatOld(neighIdx) == READY)) {
                aggStat(neighIdx)         = AGGREGATED;
                procWinner(neighIdx, 0)   = myRank;
                vertex2AggId(neighIdx, 0) = aggId;
                Kokkos::atomic_decrement(&numNonAggregated());
              }
            }
            return;
          }

          // Getting a little desperate!
          // Let us try to aggregate into a neighboring aggregate
          for (int neigh = 0; neigh < neighbors.length; ++neigh) {
            neighIdx = neighbors(neigh);
            if (lclLWGraph.isLocalNeighborVertex(neighIdx) &&
                (aggStatOld(neighIdx) == AGGREGATED)) {
              aggStat(nodeIdx)         = AGGREGATED;
              procWinner(nodeIdx, 0)   = myRank;
              vertex2AggId(nodeIdx, 0) = vertex2AggId(neighIdx, 0);
              Kokkos::atomic_decrement(&numNonAggregated());
              return;
            }
          }

          // Getting quite desperate!
          // Let us try to make a non contiguous aggregate
          if (makeNonAdjAggs) {
            for (LO otherNodeIdx = 0; otherNodeIdx < numRows; ++otherNodeIdx) {
              if ((otherNodeIdx != nodeIdx) &&
                  (aggStatOld(otherNodeIdx) == AGGREGATED)) {
                aggStat(nodeIdx)         = AGGREGATED;
                procWinner(nodeIdx, 0)   = myRank;
                vertex2AggId(nodeIdx, 0) = vertex2AggId(otherNodeIdx, 0);
                Kokkos::atomic_decrement(&numNonAggregated());
                return;
              }
            }
          }

          // Total deperation!
          // Let us make a singleton
          if (!error_on_isolated) {
            const LO aggId           = Kokkos::atomic_fetch_add(&numAggregates(), 1);
            aggStat(nodeIdx)         = AGGREGATED;
            procWinner(nodeIdx, 0)   = myRank;
            vertex2AggId(nodeIdx, 0) = aggId;
            Kokkos::atomic_decrement(&numNonAggregated());
          }
        });
    // LBV on 09/27/19: here we could copy numNonAggregated to host
    // and check for it to be equal to 0 in which case we can stop
    // looping over the different colors...
    Kokkos::deep_copy(aggStatOld, aggStat);
  }  // loop over colors

  auto numNonAggregated_h = Kokkos::create_mirror_view(numNonAggregated);
  Kokkos::deep_copy(numNonAggregated_h, numNonAggregated);
  numNonAggregatedNodes = numNonAggregated_h();
  if ((error_on_isolated) && (numNonAggregatedNodes > 0)) {
    // Error on this isolated node, as the user has requested
    std::ostringstream oss;
    oss << "MueLu::AggregationPhase3Algorithm::BuildAggregates: MueLu has detected a non-Dirichlet node that has no on-rank neighbors and is terminating (by user request). " << std::endl;
    oss << "If this error is being generated at level 0, this is due to an initial partitioning problem in your matrix." << std::endl;
    oss << "If this error is being generated at any other level, try turning on repartitioning, which may fix this problem." << std::endl;
    throw Exceptions::RuntimeError(oss.str());
  }

  // update aggregate object
  auto numAggregates_h = Kokkos::create_mirror_view(numAggregates);
  Kokkos::deep_copy(numAggregates_h, numAggregates);
  aggregates.SetNumAggregates(numAggregates_h());
}

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE3ALGORITHM_DEF_HPP_ */
