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

#ifndef MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_InterfaceAggregationAlgorithm_decl.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::InterfaceAggregationAlgorithm(RCP<const FactoryBase> const& /* graphFact */) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregatesNonKokkos(Teuchos::ParameterList const& /* params */, LWGraph const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank         = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal numLocalAggregates = aggregates.GetNumAggregates();  // number of local aggregates on current proc

  // main loop over all local rows of graph(A)
  for (int iNode1 = 0; iNode1 < nRows; ++iNode1) {
    if (aggStat[iNode1] == INTERFACE) {
      aggregates.SetIsRoot(iNode1);  // mark iNode1 as root node for new aggregate 'agg'
      int aggIndex = numLocalAggregates;
      std::vector<int> aggList;
      aggList.push_back(iNode1);
      auto neighOfINode = graph.getNeighborVertices(iNode1);

      for (int j = 0; j < neighOfINode.length; ++j) {
        LO neigh = neighOfINode(j);
        if (neigh != iNode1 && graph.isLocalNeighborVertex(neigh)) {
          if (aggStat[neigh] != AGGREGATED && aggStat[neigh] != INTERFACE &&
              aggStat[neigh] != IGNORED) {
            aggList.push_back(neigh);
          }
        }
      }

      for (size_t k = 0; k < aggList.size(); k++) {
        aggStat[aggList[k]]      = AGGREGATED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
      }
      ++numLocalAggregates;
      numNonAggregatedNodes -= aggList.size();
    }

  }  // end for

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(Teuchos::ParameterList const& /* params */, LWGraph_kokkos const& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat, LO& numNonAggregatedNodes) const {
  TEUCHOS_ASSERT(false);
}

}  // namespace MueLu

#endif /* MUELU_INTERFACEAGGREGATIONALGORITHM_DEF_HPP_ */
