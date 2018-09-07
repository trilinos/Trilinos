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
#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_


#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_AggregationStructuredAlgorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_IndexManager.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const Teuchos::ParameterList& params, const GraphBase& graph,
                  Aggregates& aggregates, std::vector<unsigned>& aggStat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    RCP<IndexManager> geoData = aggregates.GetIndexManager();
    const bool   coupled      = geoData->isAggregationCoupled();
    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    Array<LO>    ghostedCoarseNodeCoarseLIDs;
    Array<int>   ghostedCoarseNodeCoarsePIDs;
    Array<GO>    ghostedCoarseNodeCoarseGIDs;

    *out << "Extract data for ghosted nodes" << std::endl;
    geoData->getGhostedNodesData(graph.GetDomainMap(), ghostedCoarseNodeCoarseLIDs,
                                 ghostedCoarseNodeCoarsePIDs, ghostedCoarseNodeCoarseGIDs);

    LO rem, rate;
    Array<LO> ghostedIdx(3), coarseIdx(3);
    LO ghostedCoarseNodeCoarseLID, aggId;
    *out << "Loop over fine nodes and assign them to an aggregate and a rank" << std::endl;
    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {
      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, ghostedIdx[0], ghostedIdx[1], ghostedIdx[2]);

      for(int dim = 0; dim < 3; ++dim) {
        coarseIdx[dim] = ghostedIdx[dim] / geoData->getCoarseningRate(dim);
        rem    = ghostedIdx[dim] % geoData->getCoarseningRate(dim);
        if(ghostedIdx[dim] - geoData->getOffset(dim)
           < geoData->getLocalFineNodesInDir(dim) - geoData->getCoarseningEndRate(dim)) {
          rate = geoData->getCoarseningRate(dim);
        } else {
          rate = geoData->getCoarseningEndRate(dim);
        }
        if(rem > (rate / 2)) {++coarseIdx[dim];}
        if(coupled && (geoData->getStartGhostedCoarseNode(dim)*geoData->getCoarseningRate(dim)
                       > geoData->getStartIndex(dim))) {--coarseIdx[dim];}
      }

      geoData->getCoarseNodeGhostedLID(coarseIdx[0], coarseIdx[1], coarseIdx[2],
                                       ghostedCoarseNodeCoarseLID);

      aggId                 = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID];
      vertex2AggId[nodeIdx] = aggId;
      procWinner[nodeIdx]   = ghostedCoarseNodeCoarsePIDs[ghostedCoarseNodeCoarseLID];
      aggStat[nodeIdx]      = AGGREGATED;
      --numNonAggregatedNodes;

    } // Loop over fine points
  } // BuildAggregates()


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  BuildGraph(const GraphBase& graph, RCP<IndexManager>& geoData, RCP<CrsGraph>& myGraph,
             RCP<const Map>& coarseCoordinatesFineMap, RCP<const Map>& coarseCoordinatesMap) const {
    Monitor m(*this, "BuildAggregates");

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    const bool coupled = geoData->isAggregationCoupled();

    // Compute the number of coarse points needed to interpolate quantities to a fine point
    int numInterpolationPoints = 0;
    if(geoData->getInterpolationOrder() == 0) {
      numInterpolationPoints = 1;
    } else if(geoData->getInterpolationOrder() == 1) {
      // Compute 2^numDimensions using bit logic to avoid round-off errors
      numInterpolationPoints = 1 << geoData->getNumDimensions();
    }
    *out << "numInterpolationPoints=" << numInterpolationPoints << std::endl;

    Array<LO> colIndex( geoData->getNumLocalCoarseNodes() + numInterpolationPoints*
                        (geoData->getNumLocalFineNodes() - geoData->getNumLocalCoarseNodes()) );
    Array<size_t> rowPtr(geoData->getNumLocalFineNodes()+1);
    rowPtr[0] = 0;
    ArrayRCP<size_t> nnzOnRow(geoData->getNumLocalFineNodes());

    *out << "Compute prolongatorGraph data" << std::endl;
    if(geoData->getInterpolationOrder() == 0) {
      ComputeGraphDataConstant(graph, geoData, numInterpolationPoints, nnzOnRow, rowPtr, colIndex);
    } else if(geoData->getInterpolationOrder() == 1) {
      ComputeGraphDataLinear(graph, geoData, numInterpolationPoints, nnzOnRow, rowPtr, colIndex);
    }

    // Compute graph's colMap and domainMap
    RCP<Map> colMap, domainMap;
    *out << "Compute domain and column maps of the CrsGraph" << std::endl;
    if(coupled){
      *out << "Extract data for ghosted nodes" << std::endl;
      Array<LO>  ghostedCoarseNodeCoarseLIDs;
      Array<int> ghostedCoarseNodeCoarsePIDs;
      Array<GO>  ghostedCoarseNodeCoarseGIDs;
      geoData->getGhostedNodesData(graph.GetDomainMap(), ghostedCoarseNodeCoarseLIDs,
                                   ghostedCoarseNodeCoarsePIDs, ghostedCoarseNodeCoarseGIDs);

      // In this case we specify the global number of nodes on the coarse mesh
      // as well as the GIDs needed on rank.
      colMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                 geoData->getNumGlobalCoarseNodes(),
                                 ghostedCoarseNodeCoarseGIDs(),
                                 graph.GetDomainMap()->getIndexBase(),
                                 graph.GetDomainMap()->getComm(),
                                 graph.GetDomainMap()->getNode());

      LO coarseNodeIdx = 0;
      Array<GO> coarseNodeCoarseGIDs, coarseNodeFineGIDs;
      geoData->getCoarseNodesData(graph.GetDomainMap(), coarseNodeCoarseGIDs, coarseNodeFineGIDs);
      for(LO nodeIdx = 0; nodeIdx < ghostedCoarseNodeCoarseGIDs.size(); ++nodeIdx) {
        if(ghostedCoarseNodeCoarsePIDs[nodeIdx] == colMap->getComm()->getRank()) {
          coarseNodeCoarseGIDs[coarseNodeIdx] = ghostedCoarseNodeCoarseGIDs[nodeIdx];
          ++coarseNodeIdx;
        }
      }
      domainMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                    geoData->getNumGlobalCoarseNodes(),
                                    coarseNodeCoarseGIDs(),
                                    graph.GetDomainMap()->getIndexBase(),
                                    graph.GetDomainMap()->getComm(),
                                    graph.GetDomainMap()->getNode());
      coarseCoordinatesMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                               geoData->getNumGlobalCoarseNodes(),
                                               coarseNodeCoarseGIDs(),
                                               graph.GetDomainMap()->getIndexBase(),
                                               graph.GetDomainMap()->getComm(),
                                               graph.GetDomainMap()->getNode());
      coarseCoordinatesFineMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                                   geoData->getNumGlobalCoarseNodes(),
                                                   coarseNodeFineGIDs(),
                                                   graph.GetDomainMap()->getIndexBase(),
                                                   graph.GetDomainMap()->getComm(),
                                                   graph.GetDomainMap()->getNode());
    } else {
      // In this case the map will compute the global number of nodes on the coarse mesh
      // since geoData->getNumGlobalCoarseNodes() == Teuchos::OrdinalTraits<GO>::invalid()
      // and it will assign GIDs to the local coarse nodes.
      colMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                 geoData->getNumGlobalCoarseNodes(),
                                 geoData->getNumLocalCoarseNodes(),
                                 graph.GetDomainMap()->getIndexBase(),
                                 graph.GetDomainMap()->getComm(),
                                 graph.GetDomainMap()->getNode());
      domainMap = colMap;

      Array<GO> coarseNodeCoarseGIDs(geoData->getNumLocalCoarseNodes());
      Array<GO> coarseNodeFineGIDs(geoData->getNumLocalCoarseNodes());
      geoData->getCoarseNodesData(graph.GetDomainMap(), coarseNodeCoarseGIDs, coarseNodeFineGIDs);
      coarseCoordinatesMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                               geoData->getNumGlobalCoarseNodes(),
                                               geoData->getNumLocalCoarseNodes(),
                                               graph.GetDomainMap()->getIndexBase(),
                                               graph.GetDomainMap()->getComm(),
                                               graph.GetDomainMap()->getNode());
      coarseCoordinatesFineMap = MapFactory::Build(graph.GetDomainMap()->lib(),
                                                   geoData->getNumGlobalCoarseNodes(),
                                                   coarseNodeFineGIDs(),
                                                   graph.GetDomainMap()->getIndexBase(),
                                                   graph.GetDomainMap()->getComm(),
                                                   graph.GetDomainMap()->getNode());
    }

    myGraph = CrsGraphFactory::Build(graph.GetDomainMap(),
                                     colMap,
                                     nnzOnRow,
                                     Xpetra::DynamicProfile);
    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {
      myGraph->insertLocalIndices(nodeIdx, colIndex(rowPtr[nodeIdx], nnzOnRow[nodeIdx]) );
    }
    myGraph->fillComplete(domainMap, graph.GetDomainMap());

  } // BuildAggregates()


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeGraphDataConstant(const GraphBase& graph, RCP<IndexManager>& geoData,
                           const int numInterpolationPoints, ArrayRCP<size_t>& nnzOnRow,
                           Array<size_t>& rowPtr, Array<LO>& colIndex) const {

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    Array<LO>    ghostedCoarseNodeCoarseLIDs;
    Array<int>   ghostedCoarseNodeCoarsePIDs;
    Array<GO>    ghostedCoarseNodeCoarseGIDs;
    geoData->getGhostedNodesData(graph.GetDomainMap(), ghostedCoarseNodeCoarseLIDs,
                                 ghostedCoarseNodeCoarsePIDs, ghostedCoarseNodeCoarseGIDs);

    LO ghostedCoarseNodeCoarseLID, rem, rate;
    Array<LO> ghostedIdx(3), coarseIdx(3);
    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {
      // For piece-wise constant interpolation we only get one nnz per row
      nnzOnRow[nodeIdx] = Teuchos::as<size_t>(1);
      rowPtr[nodeIdx + 1] = rowPtr[nodeIdx] + 1; // These needs to change for kokkos: rowPtr[nodeIdx + 1] = nodeIdx + 1;

      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, ghostedIdx[0], ghostedIdx[1], ghostedIdx[2]);

      for(int dim = 0; dim < 3; ++dim) {
        coarseIdx[dim] = ghostedIdx[dim] / geoData->getCoarseningRate(dim);
        rem    = ghostedIdx[dim] % geoData->getCoarseningRate(dim);
        if(ghostedIdx[dim] - geoData->getOffset(dim)
           < geoData->getLocalFineNodesInDir(dim) - geoData->getCoarseningEndRate(dim)) {
          rate = geoData->getCoarseningRate(dim);
        } else {
          rate = geoData->getCoarseningEndRate(dim);
        }
        if(rem > (rate / 2)) {++coarseIdx[dim];}
        if( (geoData->getStartGhostedCoarseNode(dim)*geoData->getCoarseningRate(dim)
             > geoData->getStartIndex(dim)) && geoData->isAggregationCoupled() ) {--coarseIdx[dim];}
      }

      geoData->getCoarseNodeGhostedLID(coarseIdx[0], coarseIdx[1], coarseIdx[2],
                                       ghostedCoarseNodeCoarseLID);
      colIndex[rowPtr[nodeIdx]] = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID]; // Here too, substitute nodeIdx for rowPtr[nodeIdx]
    } // Loop over fine points

  } // ComputeGraphDataConstant()


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  ComputeGraphDataLinear(const GraphBase& graph, RCP<IndexManager>& geoData,
                         const int numInterpolationPoints, ArrayRCP<size_t>& nnzOnRow,
                         Array<size_t>& rowPtr, Array<LO>& colIndex) const {

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    const bool coupled       = geoData->isAggregationCoupled();
    const int  numDimensions = geoData->getNumDimensions();
    Array<LO> ghostedIdx(3,0);
    Array<LO> coarseIdx(3,0);
    Array<LO> ijkRem(3,0);
    int rate = 0;

    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {

      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, ghostedIdx[0], ghostedIdx[1], ghostedIdx[2]);
      for(int dim=0; dim < numDimensions; dim++){
        coarseIdx[dim] = ghostedIdx[dim] / geoData->getCoarseningRate(dim);
        ijkRem[dim]    = ghostedIdx[dim] % geoData->getCoarseningRate(dim);
        if(ghostedIdx[dim] - geoData->getOffset(dim)
           < geoData->getLocalFineNodesInDir(dim) - geoData->getCoarseningEndRate(dim)) {
          rate = geoData->getCoarseningRate(dim);
        } else {
          rate = geoData->getCoarseningEndRate(dim);
        }
        if(ijkRem[dim] > (rate / 2)) {++coarseIdx[dim];}
        if(coupled && (geoData->getStartGhostedCoarseNode(dim)*geoData->getCoarseningRate(dim)
                       > geoData->getStartIndex(dim))) {--coarseIdx[dim];}
      }

      // Fill Graph
      // Check if Fine node lies on Coarse Node
      bool allCoarse = true;
      Array<bool> isCoarse(numDimensions);
      for(int dim = 0; dim < numDimensions; ++dim) {
        isCoarse[dim] = false;
        if(ijkRem[dim] == 0)
          isCoarse[dim] = true;

        if(coupled){
          if( ghostedIdx[dim]-geoData->getOffset(dim) == geoData->getLocalFineNodesInDir(dim)-1 &&
              geoData->getMeshEdge(dim*2+1) )
            isCoarse[dim] = true;
        } else {
          if( ghostedIdx[dim]-geoData->getOffset(dim) == geoData->getLocalFineNodesInDir(dim)-1)
            isCoarse[dim] = true;
        }

        if(!isCoarse[dim])
          allCoarse = false;
      }

      if(allCoarse) {
        // Fine node lies on Coarse node, easy case, we only need the LID of the coarse node.
        geoData->getCoarseNodeGhostedLID(coarseIdx[0], coarseIdx[1], coarseIdx[2],
                                         colIndex[rowPtr[nodeIdx]]);
        nnzOnRow[nodeIdx] = Teuchos::as<size_t>(1);
        rowPtr[nodeIdx + 1] = rowPtr[nodeIdx] + 1;
      } else {
        // Harder case, we need the LIDs of all the coarse nodes contributing to the interpolation
        // at the current node.
        nnzOnRow[nodeIdx] = Teuchos::as<size_t>( numInterpolationPoints );
        rowPtr[nodeIdx + 1] = rowPtr[nodeIdx] + Teuchos::as<LO>( numInterpolationPoints );

        for(int dim = 0; dim < numDimensions; ++dim) {
          if(coarseIdx[dim] == geoData->getGhostedNodesInDir(dim) - 1)
            --coarseIdx[dim];
        }
        // Compute Coarse Node LID
        geoData->getCoarseNodeGhostedLID(    coarseIdx[0],   coarseIdx[1],   coarseIdx[2],   colIndex[ rowPtr[nodeIdx]+0]);
        geoData->getCoarseNodeGhostedLID(    coarseIdx[0]+1, coarseIdx[1],   coarseIdx[2],   colIndex[ rowPtr[nodeIdx]+1]);
        if(numDimensions > 1) {
          geoData->getCoarseNodeGhostedLID(  coarseIdx[0],   coarseIdx[1]+1, coarseIdx[2],   colIndex[ rowPtr[nodeIdx]+2]);
          geoData->getCoarseNodeGhostedLID(  coarseIdx[0]+1, coarseIdx[1]+1, coarseIdx[2],   colIndex[ rowPtr[nodeIdx]+3]);
          if(numDimensions > 2) {
            geoData->getCoarseNodeGhostedLID(coarseIdx[0],   coarseIdx[1],   coarseIdx[2]+1, colIndex[ rowPtr[nodeIdx]+4]);
            geoData->getCoarseNodeGhostedLID(coarseIdx[0]+1, coarseIdx[1],   coarseIdx[2]+1, colIndex[ rowPtr[nodeIdx]+5]);
            geoData->getCoarseNodeGhostedLID(coarseIdx[0],   coarseIdx[1]+1, coarseIdx[2]+1, colIndex[ rowPtr[nodeIdx]+6]);
            geoData->getCoarseNodeGhostedLID(coarseIdx[0]+1, coarseIdx[1]+1, coarseIdx[2]+1, colIndex[ rowPtr[nodeIdx]+7]);
          }
        }
      }
    } // Loop over fine points
  } // ComputeGraphDataLinear()

} // end namespace


#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_ */
