// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_OnePtAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_OnePtAggregationAlgorithm_decl.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::OnePtAggregationAlgorithm(RCP<const FactoryBase> const& /* graphFact */) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(Teuchos::ParameterList const& /* params */, LWGraph const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank         = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates();  // number of local aggregates on current proc
  LocalOrdinal iNode1           = 0;                              // current node

  // main loop over all local rows of graph(A)
  while (iNode1 < nRows) {
    if (aggStat[iNode1] == ONEPT) {
      aggregates.SetIsRoot(iNode1);  // mark iNode1 as root node for new aggregate 'ag'
      std::vector<int> aggList;
      aggList.push_back(iNode1);
      int aggIndex = nLocalAggregates++;

      for (size_t k = 0; k < aggList.size(); k++) {
        aggStat[aggList[k]]      = IGNORED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
      }
      numNonAggregatedNodes -= aggList.size();
    }

    iNode1++;
  }  // end while

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(Teuchos::ParameterList const& /* params */,
                    LWGraph_kokkos const& graph,
                    Aggregates& aggregates,
                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggstat,
                    LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  Monitor m(*this, "BuildAggregates");

  typename Kokkos::View<unsigned*, device_type>::HostMirror aggstatHost = Kokkos::create_mirror(aggstat);
  Kokkos::deep_copy(aggstatHost, aggstat);
  std::vector<unsigned> aggStat;
  aggStat.resize(aggstatHost.extent(0));
  for (size_t idx = 0; idx < aggstatHost.extent(0); ++idx) {
    aggStat[idx] = aggstatHost(idx);
  }

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank         = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates();  // number of local aggregates on current proc
  LocalOrdinal iNode1           = 0;                              // current node

  // main loop over all local rows of graph(A)
  while (iNode1 < nRows) {
    if (aggStat[iNode1] == ONEPT) {
      aggregates.SetIsRoot(iNode1);  // mark iNode1 as root node for new aggregate 'ag'
      std::vector<int> aggList;
      aggList.push_back(iNode1);
      int aggIndex = nLocalAggregates++;

      // finalize aggregate
      for (size_t k = 0; k < aggList.size(); k++) {
        aggStat[aggList[k]]      = IGNORED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
      }
      numNonAggregatedNodes -= aggList.size();
    }

    iNode1++;
  }  // end while

  for (size_t idx = 0; idx < aggstatHost.extent(0); ++idx) {
    aggstatHost(idx) = aggStat[idx];
  }
  Kokkos::deep_copy(aggstat, aggstatHost);

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);
}

}  // namespace MueLu

#endif /* MUELU_ONEPTAGGREGATIONALGORITHM_DEF_HPP_ */
