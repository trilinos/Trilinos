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
#ifndef MUELU_AGGREGATIONPHASE3ALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONPHASE3ALGORITHM_KOKKOS_DEF_HPP

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase3Algorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

// #include "Kokkos_Core.hpp"

namespace MueLu {

// Try to stick unaggregated nodes into a neighboring aggregate if they are
// not already too big. Otherwise, make a new aggregate
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase3Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const ParameterList& params,
                    const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
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
void AggregationPhase3Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesRandom(const ParameterList& params,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                          LO& numNonAggregatedNodes) const {
  bool error_on_isolated = params.get<bool>("aggregation: error on nodes with no on-rank neighbors");
  bool makeNonAdjAggs    = params.get<bool>("aggregation: phase3 avoid singletons");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner    = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors        = aggregates.GetGraphColors();
  const LO numColors = aggregates.GetGraphNumColors();

  auto lclLWGraph = graph.getLocalLWGraph();

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

#endif  // MUELU_AGGREGATIONPHASE3ALGORITHM_KOKKOS_DEF_HPP
