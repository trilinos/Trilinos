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
#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_KOKKOS_DEF_HPP

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_IndexManager_kokkos.hpp"
#include "MueLu_AggregationStructuredAlgorithm_kokkos_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const Teuchos::ParameterList& /* params */, const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    Kokkos::View<unsigned*, device_type>& aggStat,
                    LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregates");

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  RCP<IndexManager_kokkos> geoData = aggregates.GetIndexManagerKokkos();
  const LO numLocalFineNodes       = geoData->getNumLocalFineNodes();
  const LO numCoarseNodes          = geoData->getNumCoarseNodes();
  LOVectorView vertex2AggId        = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  LOVectorView procWinner          = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);

  *out << "Loop over fine nodes and assign them to an aggregate and a rank" << std::endl;
  LO numAggregatedNodes;
  fillAggregatesFunctor fillAggregates(geoData,
                                       graph.GetComm()->getRank(),
                                       aggStat,
                                       vertex2AggId,
                                       procWinner);
  Kokkos::parallel_reduce("StructuredAggregation: fill aggregates data",
                          Kokkos::RangePolicy<execution_space>(0, numLocalFineNodes),
                          fillAggregates,
                          numAggregatedNodes);

  *out << "numCoarseNodes= " << numCoarseNodes
       << ", numAggregatedNodes= " << numAggregatedNodes << std::endl;
  numNonAggregatedNodes = numNonAggregatedNodes - numAggregatedNodes;

}  // BuildAggregates()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    BuildGraph(const LWGraph_kokkos& graph, RCP<IndexManager_kokkos>& geoData, const LO dofsPerNode,
               RCP<CrsGraph>& myGraph) const {
  Monitor m(*this, "BuildGraphP");

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  // Compute the number of coarse points needed to interpolate quantities to a fine point
  int numInterpolationPoints = 0;
  if (geoData->getInterpolationOrder() == 0) {
    numInterpolationPoints = 1;
  } else if (geoData->getInterpolationOrder() == 1) {
    // Compute 2^numDimensions using bit logic to avoid round-off errors from std::pow()
    numInterpolationPoints = 1 << geoData->getNumDimensions();
  }
  *out << "numInterpolationPoints=" << numInterpolationPoints << std::endl;

  const LO numLocalFineNodes = geoData->getNumLocalFineNodes();
  const LO numCoarseNodes    = geoData->getNumCoarseNodes();
  const LO numNnzEntries     = dofsPerNode * (numCoarseNodes + numInterpolationPoints * (numLocalFineNodes - numCoarseNodes));

  non_const_row_map_type rowPtr("Prolongator graph, rowPtr", dofsPerNode * (numLocalFineNodes + 1));
  entries_type colIndex("Prolongator graph, colIndices", numNnzEntries);

  *out << "Compute prolongatorGraph data" << std::endl;
  if (geoData->getInterpolationOrder() == 0) {
    computeGraphDataConstantFunctor computeGraphData(geoData,
                                                     numCoarseNodes,
                                                     dofsPerNode,
                                                     geoData->getCoarseningRates(),
                                                     geoData->getCoarseningEndRates(),
                                                     geoData->getLocalFineNodesPerDir(),
                                                     rowPtr,
                                                     colIndex);
    Kokkos::parallel_for("Structured Aggregation: compute loca graph data",
                         Kokkos::RangePolicy<execution_space>(0, numLocalFineNodes),
                         computeGraphData);
  } else if (geoData->getInterpolationOrder() == 1) {
    // Note, lbv 2018-11-08: in the piece-wise linear case I am computing the rowPtr
    // using a parallel scan, it might be possible to do something faster than that
    // by including this calculation in computeGraphDataLinearFunctor but at the moment
    // all the ideas I have include a bunch of if statements which I would like to avoid.
    computeGraphRowPtrFunctor computeGraphRowPtr(geoData,
                                                 dofsPerNode,
                                                 numInterpolationPoints,
                                                 numLocalFineNodes,
                                                 geoData->getCoarseningRates(),
                                                 geoData->getLocalFineNodesPerDir(),
                                                 rowPtr);
    Kokkos::parallel_scan("Structured Aggregation: compute rowPtr for prolongator graph",
                          Kokkos::RangePolicy<execution_space>(0, numLocalFineNodes + 1),
                          computeGraphRowPtr);

    computeGraphDataLinearFunctor computeGraphData(geoData,
                                                   geoData->getNumDimensions(),
                                                   numCoarseNodes,
                                                   dofsPerNode,
                                                   numInterpolationPoints,
                                                   geoData->getCoarseningRates(),
                                                   geoData->getCoarseningEndRates(),
                                                   geoData->getLocalFineNodesPerDir(),
                                                   geoData->getCoarseNodesPerDir(),
                                                   rowPtr,
                                                   colIndex);
    Kokkos::parallel_for("Structured Aggregation: compute loca graph data",
                         Kokkos::RangePolicy<execution_space>(0, numLocalFineNodes),
                         computeGraphData);
  }

  local_graph_type myLocalGraph(colIndex, rowPtr);

  // Compute graph's colMap and domainMap
  RCP<Map> colMap, domainMap;
  *out << "Compute domain and column maps of the CrsGraph" << std::endl;
  colMap    = MapFactory::Build(graph.GetDomainMap()->lib(),
                                Teuchos::OrdinalTraits<GO>::invalid(),
                                numCoarseNodes,
                                graph.GetDomainMap()->getIndexBase(),
                                graph.GetDomainMap()->getComm());
  domainMap = colMap;

  myGraph = CrsGraphFactory::Build(myLocalGraph, graph.GetDomainMap(), colMap,
                                   colMap, graph.GetDomainMap());

}  // BuildGraph()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    fillAggregatesFunctor::fillAggregatesFunctor(RCP<IndexManager_kokkos> geoData,
                                                 const int myRank,
                                                 Kokkos::View<unsigned*, device_type> aggStat,
                                                 LOVectorView vertex2AggID,
                                                 LOVectorView procWinner)
  : geoData_(*geoData)
  , myRank_(myRank)
  , aggStat_(aggStat)
  , vertex2AggID_(vertex2AggID)
  , procWinner_(procWinner) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
KOKKOS_INLINE_FUNCTION void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    fillAggregatesFunctor::operator()(const LO nodeIdx, LO& lNumAggregatedNodes) const {
  // Compute coarse ID associated with fine LID
  LO rem, rate;
  LO coarseNodeCoarseLID;
  LO nodeFineTuple[3], coarseIdx[3];
  auto coarseRate       = geoData_.getCoarseningRates();
  auto endRate          = geoData_.getCoarseningEndRates();
  auto lFineNodesPerDir = geoData_.getLocalFineNodesPerDir();
  // Compute coarse ID associated with fine LID
  geoData_.getFineLID2FineTuple(nodeIdx, nodeFineTuple);

  for (int dim = 0; dim < 3; ++dim) {
    coarseIdx[dim] = nodeFineTuple[dim] / coarseRate(dim);
    rem            = nodeFineTuple[dim] % coarseRate(dim);
    rate           = (nodeFineTuple[dim] < lFineNodesPerDir(dim) - endRate(dim)) ? coarseRate(dim) : endRate(dim);
    if (rem > (rate / 2)) {
      ++coarseIdx[dim];
    }
  }

  geoData_.getCoarseTuple2CoarseLID(coarseIdx[0], coarseIdx[1], coarseIdx[2],
                                    coarseNodeCoarseLID);

  vertex2AggID_(nodeIdx, 0) = coarseNodeCoarseLID;
  procWinner_(nodeIdx, 0)   = myRank_;
  aggStat_(nodeIdx)         = AGGREGATED;
  ++lNumAggregatedNodes;

}  // fillAggregatesFunctor::operator()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphDataConstantFunctor::
        computeGraphDataConstantFunctor(RCP<IndexManager_kokkos> geoData,
                                        const LO NumGhostedNodes,
                                        const LO dofsPerNode,
                                        constIntTupleView coarseRate,
                                        constIntTupleView endRate,
                                        constLOTupleView lFineNodesPerDir,
                                        non_const_row_map_type rowPtr,
                                        entries_type colIndex)
  : geoData_(*geoData)
  , numGhostedNodes_(NumGhostedNodes)
  , dofsPerNode_(dofsPerNode)
  , coarseRate_(coarseRate)
  , endRate_(endRate)
  , lFineNodesPerDir_(lFineNodesPerDir)
  , rowPtr_(rowPtr)
  , colIndex_(colIndex) {
}  // computeGraphDataConstantFunctor()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
KOKKOS_INLINE_FUNCTION void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphDataConstantFunctor::operator()(const LO nodeIdx) const {
  LO nodeFineTuple[3]   = {0, 0, 0};
  LO nodeCoarseTuple[3] = {0, 0, 0};

  // Compute ghosted tuple associated with fine LID
  geoData_.getFineLID2FineTuple(nodeIdx, nodeFineTuple);

  // Compute coarse tuple associated with fine point
  // then overwrite it with tuple associated with aggregate
  LO rem, rate, coarseNodeCoarseLID;
  for (int dim = 0; dim < 3; ++dim) {
    nodeCoarseTuple[dim] = nodeFineTuple[dim] / coarseRate_(dim);
    rem                  = nodeFineTuple[dim] % coarseRate_(dim);
    if (nodeFineTuple[dim] < (lFineNodesPerDir_(dim) - endRate_(dim))) {
      rate = coarseRate_(dim);
    } else {
      rate = endRate_(dim);
    }
    if (rem > (rate / 2)) {
      ++nodeCoarseTuple[dim];
    }
  }

  // get LID associted with aggregate
  geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1], nodeCoarseTuple[2],
                                    coarseNodeCoarseLID);

  // store data into CrsGraph taking care of multiple dofs case
  for (LO dof = 0; dof < dofsPerNode_; ++dof) {
    rowPtr_(nodeIdx * dofsPerNode_ + dof + 1) = nodeIdx * dofsPerNode_ + dof + 1;
    colIndex_(nodeIdx * dofsPerNode_ + dof)   = coarseNodeCoarseLID * dofsPerNode_ + dof;
  }

}  // computeGraphDataConstantFunctor::operator()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphRowPtrFunctor::computeGraphRowPtrFunctor(RCP<IndexManager_kokkos> geoData,
                                                         const LO dofsPerNode,
                                                         const int numInterpolationPoints,
                                                         const LO numLocalRows,
                                                         constIntTupleView coarseRate,
                                                         constLOTupleView lFineNodesPerDir,
                                                         non_const_row_map_type rowPtr)
  : geoData_(*geoData)
  , dofsPerNode_(dofsPerNode)
  , numInterpolationPoints_(numInterpolationPoints)
  , numLocalRows_(numLocalRows)
  , coarseRate_(coarseRate)
  , lFineNodesPerDir_(lFineNodesPerDir)
  , rowPtr_(rowPtr) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
KOKKOS_INLINE_FUNCTION void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphRowPtrFunctor::operator()(const LO rowIdx, GO& update, const bool final) const {
  if (final) {
    // Kokkos uses a multipass algorithm to implement scan.
    // Only update the array on the final pass.  Updating the
    // array before changing 'update' means that we do an
    // exclusive scan.  Update the array after for an inclusive
    // scan.
    rowPtr_(rowIdx) = update;
  }
  if (rowIdx < numLocalRows_) {
    LO nodeIdx          = rowIdx / dofsPerNode_;
    bool allCoarse      = true;
    LO nodeFineTuple[3] = {0, 0, 0};
    geoData_.getFineLID2FineTuple(nodeIdx, nodeFineTuple);
    for (int dim = 0; dim < 3; ++dim) {
      const LO rem = nodeFineTuple[dim] % coarseRate_(dim);

      // Check if Fine node lies on Coarse Node
      allCoarse = (allCoarse && ((rem == 0) || (nodeFineTuple[dim] == lFineNodesPerDir_(dim) - 1)));
    }
    update += (allCoarse ? 1 : numInterpolationPoints_);
  }
}  // computeGraphRowPtrFunctor::operator()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphDataLinearFunctor::computeGraphDataLinearFunctor(RCP<IndexManager_kokkos> geoData,
                                                                 const int numDimensions,
                                                                 const LO numGhostedNodes,
                                                                 const LO dofsPerNode,
                                                                 const int numInterpolationPoints,
                                                                 constIntTupleView coarseRate,
                                                                 constIntTupleView endRate,
                                                                 constLOTupleView lFineNodesPerDir,
                                                                 constLOTupleView ghostedNodesPerDir,
                                                                 non_const_row_map_type rowPtr,
                                                                 entries_type colIndex)
  : geoData_(*geoData)
  , numDimensions_(numDimensions)
  , numGhostedNodes_(numGhostedNodes)
  , dofsPerNode_(dofsPerNode)
  , numInterpolationPoints_(numInterpolationPoints)
  , coarseRate_(coarseRate)
  , endRate_(endRate)
  , lFineNodesPerDir_(lFineNodesPerDir)
  , ghostedNodesPerDir_(ghostedNodesPerDir)
  , rowPtr_(rowPtr)
  , colIndex_(colIndex) {
}  // computeGraphDataLinearFunctor()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
KOKKOS_INLINE_FUNCTION void AggregationStructuredAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    computeGraphDataLinearFunctor::operator()(const LO nodeIdx) const {
  LO nodeFineTuple[3]   = {0, 0, 0};
  LO nodeCoarseTuple[3] = {0, 0, 0};

  // Compute coarse ID associated with fine LID
  geoData_.getFineLID2FineTuple(nodeIdx, nodeFineTuple);

  LO coarseNodeCoarseLID;
  bool allCoarse = false;
  for (int dim = 0; dim < 3; ++dim) {
    nodeCoarseTuple[dim] = nodeFineTuple[dim] / coarseRate_(dim);
  }
  if (rowPtr_(nodeIdx + 1) == rowPtr_(nodeIdx) + 1) {
    allCoarse = true;
  }

  geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1], nodeCoarseTuple[2],
                                    coarseNodeCoarseLID);

  if (allCoarse) {
    // Fine node lies on Coarse node, easy case, we only need the LID of the coarse node.
    for (LO dof = 0; dof < dofsPerNode_; ++dof) {
      colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof)) = coarseNodeCoarseLID * dofsPerNode_ + dof;
    }
  } else {
    for (int dim = 0; dim < numDimensions_; ++dim) {
      if (nodeCoarseTuple[dim] == ghostedNodesPerDir_(dim) - 1) {
        --nodeCoarseTuple[dim];
      }
    }
    // Compute Coarse Node LID
    // Note lbv 10-06-2018: it is likely benefitial to remove the two if statments and somehow
    // find out the number of dimensions before calling the opertor() of the functor.
    for (LO dof = 0; dof < dofsPerNode_; ++dof) {
      geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1], nodeCoarseTuple[2], colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 0));
      geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0] + 1, nodeCoarseTuple[1], nodeCoarseTuple[2], colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 1));
      if (numDimensions_ > 1) {
        geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1] + 1, nodeCoarseTuple[2], colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 2));
        geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0] + 1, nodeCoarseTuple[1] + 1, nodeCoarseTuple[2], colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 3));
        if (numDimensions_ > 2) {
          geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1], nodeCoarseTuple[2] + 1, colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 4));
          geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0] + 1, nodeCoarseTuple[1], nodeCoarseTuple[2] + 1, colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 5));
          geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0], nodeCoarseTuple[1] + 1, nodeCoarseTuple[2] + 1, colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 6));
          geoData_.getCoarseTuple2CoarseLID(nodeCoarseTuple[0] + 1, nodeCoarseTuple[1] + 1, nodeCoarseTuple[2] + 1, colIndex_(rowPtr_(nodeIdx * dofsPerNode_ + dof) + 7));
        }
      }
    }
  }
}  // computeGraphDataLinearFunctor::operator()

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_ */
