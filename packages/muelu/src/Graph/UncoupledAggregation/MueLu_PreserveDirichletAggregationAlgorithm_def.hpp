// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_PreserveDirichletAggregationAlgorithm_decl.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void PreserveDirichletAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(Teuchos::ParameterList const& params, LWGraph const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

  bool preserve = params.get<bool>("aggregation: preserve Dirichlet points");

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  LO numLocalAggregates = aggregates.GetNumAggregates();

  for (LO i = 0; i < numRows; i++)
    if (aggStat[i] == BOUNDARY) {
      aggStat[i] = IGNORED;
      numNonAggregatedNodes--;

      if (preserve) {
        aggregates.SetIsRoot(i);

        vertex2AggId[i] = numLocalAggregates++;
        procWinner[i]   = myRank;
      }
    }

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void PreserveDirichletAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(Teuchos::ParameterList const& params,
                    LWGraph_kokkos const& graph,
                    Aggregates& aggregates,
                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                    LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  Monitor m(*this, "BuildAggregates");
  using local_ordinal_type = typename LWGraph_kokkos::local_ordinal_type;

  // Extract parameters and data from:
  // 1) the parameter list
  const bool preserve = params.get<bool>("aggregation: preserve Dirichlet points");

  // 2) the amalgamated graph
  const LO numNodes = graph.GetNodeNumVertices();
  const int myRank  = graph.GetComm()->getRank();

  // 3) the aggregates
  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner   = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);

  // A view is needed to count on the fly the current number of local aggregates
  Kokkos::View<LO, device_type> aggCount("aggCount");
  if (preserve) {
    Kokkos::deep_copy(aggCount, aggregates.GetNumAggregates());
  }
  Kokkos::parallel_for(
      "MueLu - PreserveDirichlet: tagging ignored nodes",
      Kokkos::RangePolicy<local_ordinal_type, execution_space>(0, numNodes),
      KOKKOS_LAMBDA(const local_ordinal_type nodeIdx) {
        if (aggStat(nodeIdx) == BOUNDARY) {
          aggStat(nodeIdx) = IGNORED;
          const LO aggIdx  = Kokkos::atomic_fetch_add(&aggCount(), 1);

          if (preserve) {
            // aggregates.SetIsRoot(nodeIdx);

            vertex2AggId(nodeIdx, 0) = aggIdx;
            procWinner(nodeIdx, 0)   = myRank;
          }
        }
      });
  typename Kokkos::View<LO, device_type>::HostMirror aggCount_h = Kokkos::create_mirror_view(aggCount);
  Kokkos::deep_copy(aggCount_h, aggCount);
  // In this phase the number of new aggregates is the same
  // as the number of newly aggregated nodes.
  numNonAggregatedNodes -= (aggCount_h() - aggregates.GetNumAggregates());

  // update aggregate object
  if (preserve) {
    aggregates.SetNumAggregates(aggCount_h());
  }
}

}  // namespace MueLu

#endif /* MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_DEF_HPP_ */
