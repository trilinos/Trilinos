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
#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP

#include <queue>
#include <vector>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase1Algorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

#include "Kokkos_Sort.hpp"
#include <Kokkos_ScatterView.hpp>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const Teuchos::ParameterList& params,
                    const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                    LO& numNonAggregatedNodes) const {
  int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
  int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");

  TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate,
                             Exceptions::RuntimeError,
                             "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

  // Distance-2 gives less control than serial uncoupled phase 1
  // no custom row reordering because would require making deep copy
  // of local matrix entries and permuting it can only enforce
  // max aggregate size
  {
    if (params.get<bool>("aggregation: deterministic")) {
      Monitor m(*this, "BuildAggregatesDeterministic");
      BuildAggregatesDeterministic(maxNodesPerAggregate, graph,
                                   aggregates, aggStat, numNonAggregatedNodes);
    } else {
      Monitor m(*this, "BuildAggregatesRandom");
      BuildAggregatesRandom(maxNodesPerAggregate, graph,
                            aggregates, aggStat, numNonAggregatedNodes);
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesRandom(const LO maxAggSize,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                          LO& numNonAggregatedNodes) const {
  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // Extract data from aggregates
  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner   = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors       = aggregates.GetGraphColors();

  auto lclLWGraph = graph.getLocalLWGraph();

  LO numAggregatedNodes = 0;
  LO numLocalAggregates = aggregates.GetNumAggregates();
  Kokkos::View<LO, device_type> aggCount("aggCount");
  Kokkos::deep_copy(aggCount, numLocalAggregates);
  Kokkos::parallel_for(
      "Aggregation Phase 1: initial reduction over color == 1",
      Kokkos::RangePolicy<LO, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO nodeIdx) {
        if (colors(nodeIdx) == 1 && aggStat(nodeIdx) == READY) {
          const LO aggIdx          = Kokkos::atomic_fetch_add(&aggCount(), 1);
          vertex2AggId(nodeIdx, 0) = aggIdx;
          aggStat(nodeIdx)         = AGGREGATED;
          procWinner(nodeIdx, 0)   = myRank;
        }
      });
  // Truely we wish to compute: numAggregatedNodes = aggCount - numLocalAggregates
  // before updating the value of numLocalAggregates.
  // But since we also do not want to create a host mirror of aggCount we do some trickery...
  numAggregatedNodes -= numLocalAggregates;
  Kokkos::deep_copy(numLocalAggregates, aggCount);
  numAggregatedNodes += numLocalAggregates;

  // Compute the initial size of the aggregates.
  // Note lbv 12-21-17: I am pretty sure that the aggregates will always be of size 1
  //                    at this point so we could simplify the code below a lot if this
  //                    assumption is correct...
  Kokkos::View<LO*, device_type> aggSizesView("aggSizes", numLocalAggregates);
  {
    // Here there is a possibility that two vertices assigned to two different threads contribute
    // to the same aggregate if somethings happened before phase 1?
    auto aggSizesScatterView = Kokkos::Experimental::create_scatter_view(aggSizesView);
    Kokkos::parallel_for(
        "Aggregation Phase 1: compute initial aggregates size",
        Kokkos::RangePolicy<LO, execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO nodeIdx) {
          auto aggSizesScatterViewAccess = aggSizesScatterView.access();
          if (vertex2AggId(nodeIdx, 0) >= 0)
            aggSizesScatterViewAccess(vertex2AggId(nodeIdx, 0)) += 1;
        });
    Kokkos::Experimental::contribute(aggSizesView, aggSizesScatterView);
  }

  LO tmpNumAggregatedNodes = 0;
  Kokkos::parallel_reduce(
      "Aggregation Phase 1: main parallel_reduce over aggSizes",
      Kokkos::RangePolicy<size_t, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const size_t nodeIdx, LO& lNumAggregatedNodes) {
        if (colors(nodeIdx) != 1 && (aggStat(nodeIdx) == READY || aggStat(nodeIdx) == NOTSEL)) {
          // Get neighbors of vertex i and look for local, aggregated,
          // color 1 neighbor (valid root).
          auto neighbors = lclLWGraph.getNeighborVertices(nodeIdx);
          for (LO j = 0; j < neighbors.length; ++j) {
            auto nei = neighbors.colidx(j);
            if (lclLWGraph.isLocalNeighborVertex(nei) && colors(nei) == 1 && aggStat(nei) == AGGREGATED) {
              // This atomic guarentees that any other node trying to
              // join aggregate agg has the correct size.
              LO agg           = vertex2AggId(nei, 0);
              const LO aggSize = Kokkos::atomic_fetch_add(&aggSizesView(agg),
                                                          1);
              if (aggSize < maxAggSize) {
                // assign vertex i to aggregate with root j
                vertex2AggId(nodeIdx, 0) = agg;
                procWinner(nodeIdx, 0)   = myRank;
                aggStat(nodeIdx)         = AGGREGATED;
                ++lNumAggregatedNodes;
                break;
              } else {
                // Decrement back the value of aggSizesView(agg)
                Kokkos::atomic_decrement(&aggSizesView(agg));
              }
            }
          }
        }
        // if(aggStat(nodeIdx) != AGGREGATED) {
        //   lNumNonAggregatedNodes++;
        if (aggStat(nodeIdx) == NOTSEL) {
          aggStat(nodeIdx) = READY;
        }
        // }
      },
      tmpNumAggregatedNodes);
  numAggregatedNodes += tmpNumAggregatedNodes;
  numNonAggregatedNodes -= numAggregatedNodes;

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesDeterministic(const LO maxAggSize,
                                 const LWGraph_kokkos& graph,
                                 Aggregates& aggregates,
                                 Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                                 LO& numNonAggregatedNodes) const {
  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner   = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors       = aggregates.GetGraphColors();

  auto lclLWGraph = graph.getLocalLWGraph();

  LO numLocalAggregates = aggregates.GetNumAggregates();
  Kokkos::View<LO, device_type> numLocalAggregatesView("Num aggregates");
  {
    auto h_nla = Kokkos::create_mirror_view(numLocalAggregatesView);
    h_nla()    = numLocalAggregates;
    Kokkos::deep_copy(numLocalAggregatesView, h_nla);
  }

  Kokkos::View<LO*, device_type> newRoots("New root LIDs", numNonAggregatedNodes);
  Kokkos::View<LO, device_type> numNewRoots("Number of new aggregates of current color");
  auto h_numNewRoots = Kokkos::create_mirror_view(numNewRoots);

  // first loop build the set of new roots
  Kokkos::parallel_for(
      "Aggregation Phase 1: building list of new roots",
      Kokkos::RangePolicy<execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO i) {
        if (colors(i) == 1 && aggStat(i) == READY) {
          // i will become a root
          newRoots(Kokkos::atomic_fetch_add(&numNewRoots(), 1)) = i;
        }
      });
  Kokkos::deep_copy(h_numNewRoots, numNewRoots);
  // sort new roots by LID to guarantee determinism in agg IDs
  Kokkos::sort(newRoots, 0, h_numNewRoots());
  LO numAggregated = 0;
  Kokkos::parallel_reduce(
      "Aggregation Phase 1: aggregating nodes",
      Kokkos::RangePolicy<execution_space>(0, h_numNewRoots()),
      KOKKOS_LAMBDA(const LO rootIndex, LO& lnumAggregated) {
        LO root               = newRoots(rootIndex);
        LO aggID              = numLocalAggregatesView() + rootIndex;
        LO aggSize            = 1;
        vertex2AggId(root, 0) = aggID;
        procWinner(root, 0)   = myRank;
        aggStat(root)         = AGGREGATED;
        auto neighOfRoot      = lclLWGraph.getNeighborVertices(root);
        for (LO n = 0; n < neighOfRoot.length; n++) {
          LO neigh = neighOfRoot(n);
          if (lclLWGraph.isLocalNeighborVertex(neigh) && aggStat(neigh) == READY) {
            // add neigh to aggregate
            vertex2AggId(neigh, 0) = aggID;
            procWinner(neigh, 0)   = myRank;
            aggStat(neigh)         = AGGREGATED;
            aggSize++;
            if (aggSize == maxAggSize) {
              // can't add any more nodes
              break;
            }
          }
        }
        lnumAggregated += aggSize;
      },
      numAggregated);
  numNonAggregatedNodes -= numAggregated;
  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates + h_numNewRoots());
}

}  // namespace MueLu

#endif  // MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP
