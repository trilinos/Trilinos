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
#ifndef MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_StructuredAggregationFactory_decl.hpp"

#include "MueLu_OnePtAggregationAlgorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  StructuredAggregationFactory() : bDefinitionPhase_(true)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
    SET_VALID_ENTRY("aggregation: allow user-specified singletons");
    SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
#undef  SET_VALID_ENTRY

    // general variables needed in AggregationFactory
    validParamList->set<RCP<const FactoryBase> >("Coordinates",             Teuchos::null,
                                                 "Generating factory of problem coordinates");
    validParamList->set<int>                    ("aggregation: number of spatial dimensions", 3,
                                                  "The number of spatial dimensions in the problem");
    validParamList->set<std::string>            ("aggregation: coarsening rate", "{3}",
                                                  "Coarsening rate per spatial dimensions");
    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<std::string >           ("meshLayout",              "Global Lexicographic",
                                                 "Type of mesh ordering");
    validParamList->set<RCP<const FactoryBase> >("meshData",                Teuchos::null,
                                                 "Mesh ordering associated data");

    // special variables necessary for OnePtAggregationAlgorithm
    validParamList->set<std::string>            ("OnePt aggregate map name",         "",
                                                 "Name of input map for single node aggregates. (default='')");
    validParamList->set<std::string>            ("OnePt aggregate map factory",      "",
                                                 "Generating factory of (DOF) map for single node aggregates.");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Coordinates");

    // Request the global number of nodes per dimensions
    if(currentLevel.GetLevelID() == 0) {
      if(currentLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
        currentLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "gNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "gNodesPerDim");
    }

    // Request the local number of nodes per dimensions
    if(currentLevel.GetLevelID() == 0) {
      if(currentLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
        currentLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "lNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "lNodesPerDim");
    }

    const ParameterList& pL = GetParameterList();

    // request special data necessary for OnePtAggregationAlgorithm
    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
    if (mapOnePtName.length() > 0) {
      std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
      if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
        currentLevel.DeclareInput(mapOnePtName, NoFactory::get());
      } else {
        RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
        currentLevel.DeclareInput(mapOnePtName, mapOnePtFact.get());
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    ParameterList pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    // General problem informations are gathered from data stored in the problem matix.
    typedef typename Xpetra::MultiVector<double, LO, GO, NO> xdMV;
    RCP<const xdMV> Coordinates = Get< RCP<const xdMV> >(currentLevel, "Coordinates");
    TEUCHOS_TEST_FOR_EXCEPTION(Coordinates == Teuchos::null, Exceptions::RuntimeError,
                               "Coordinates cannot be accessed from fine level!");
    RCP<const Map> coordMap     = Coordinates->getMap();
    const int myRank            = coordMap->getComm()->getRank();
    const int numRanks          = coordMap->getComm()->getSize();

    // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
    // obtain a nodeMap.
    RCP<GeometricData> geoData = rcp(new GeometricData{});
    std::string meshLayout = pL.get<std::string>("meshLayout");
    const int numDimensions = Coordinates->getNumVectors();
    Array<GO> gFineNodesPerDir(3);
    Array<LO> lFineNodesPerDir(3);
    if(currentLevel.GetLevelID() == 0) {
      // On level 0, data is provided by applications and has no associated factory.
      gFineNodesPerDir = currentLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
      lFineNodesPerDir = currentLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    } else {
      // On level > 0, data is provided directly by generating factories.
      gFineNodesPerDir = Get<Array<GO> >(currentLevel, "gNodesPerDim");
      lFineNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
    }

    // Get the coarsening rate
    std::string coarseningRate = pL.get<std::string>("aggregation: coarsening rate");
    Teuchos::Array<LO> coarseRate;
    try {
      coarseRate = Teuchos::fromStringToArray<LO>(coarseningRate);
    } catch(const Teuchos::InvalidArrayStringRepresentation e) {
      GetOStream(Errors,-1) << " *** \"aggregation: coarsening rate\" must be a string convertible into an array! *** "
                            << std::endl;
      throw e;
    }
    TEUCHOS_TEST_FOR_EXCEPTION((coarseRate.size() > 1) &&
                               (coarseRate.size() < geoData->numDimensions),
                               Exceptions::RuntimeError,
                               "\"aggregation: coarsening rate\" must have at least as many"
                               " components as the number of spatial dimensions in the problem.");

    Array<GO> meshData;
    if(meshLayout == "Local Lexicographic") {
      if(currentLevel.GetLevelID() == 0) {
        // On level 0, data is provided by applications and has no associated factory.
        meshData = currentLevel.Get<Array<GO> >("meshData", NoFactory::get());
        TEUCHOS_TEST_FOR_EXCEPTION(meshData.empty() == true, Exceptions::RuntimeError,
                                   "The meshData array is empty, somehow the input for structured"
                                   " aggregation are not captured correctly.");
      } else {
        // On level > 0, data is provided directly by generating factories.
        meshData = Get<Array<GO> >(currentLevel, "meshData");
      }
    }

    geoData->computeGeometricParameters(meshLayout, numDimensions, myRank, numRanks,
                                        gFineNodesPerDir, lFineNodesPerDir, coarseRate, meshData,
                                        coordMap->getMinGlobalIndex());

    TEUCHOS_TEST_FOR_EXCEPTION(Coordinates->getLocalLength()
                               != static_cast<size_t>(geoData->lNumFineNodes),
                               Exceptions::RuntimeError,
                               "The local number of elements in Coordinates is not equal to the"
                               " number of nodes given by: lNodesPerDim!");
    TEUCHOS_TEST_FOR_EXCEPTION(Coordinates->getGlobalLength()
                               != static_cast<size_t>(geoData->gNumFineNodes),
                               Exceptions::RuntimeError,
                               "The global number of elements in Coordinates is not equal to the"
                               " number of nodes given by: gNodesPerDim!");

    // Build
    RCP<Aggregates> aggregates = rcp(new Aggregates(coordMap));
    aggregates->setObjectLabel("ST");
    aggregates->AggregatesCrossProcessors(true);

    // construct aggStat information
    std::vector<unsigned> aggStat(geoData->lNumFineNodes, READY);

    LO numNonAggregatedNodes = geoData->lNumFineNodes;
    GO numGlobalAggregatedPrev = 0, numGlobalAggsPrev = 0;
    if(geoData->meshLayout == "Global Lexicographic") {
      GlobalLexicographicLayout(coordMap, geoData, aggregates, aggStat, numNonAggregatedNodes);
    } else if(geoData->meshLayout == "Local Lexicographic") {
      LocalLexicographicLayout(coordMap, geoData, aggregates, aggStat, numNonAggregatedNodes);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError,
                               "MueLu::StructuredAggregationFactory::Build: Leftover nodes found! Error!");

    // aggregates->AggregatesCrossProcessors(false);
    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics1) << aggregates->description() << std::endl;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GlobalLexicographicLayout(const RCP<const Map> coordMap, RCP<GeometricData> geoData,
                            RCP<Aggregates> aggregates, std::vector<unsigned>& aggStat,
                            LO& numNonAggregatedNodes) const {

    aggregates->SetNumAggregates(geoData->lNumCoarseNodes);

    // Find the GIDs, LIDs and PIDs of the coarse points on the fine mesh and coarse
    // mesh as this data will be used to fill vertex2AggId and procWinner vectors.
    Array<GO> lCoarseNodeCoarseGIDs(geoData->lNumCoarseNodes),
      lCoarseNodeFineGIDs(geoData->lNumCoarseNodes);
    Array<GO> ghostedCoarseNodeCoarseGIDs(geoData->numGhostedCoarseNodes),
      ghostedCoarseNodeFineGIDs(geoData->numGhostedCoarseNodes);
    Array<LO> ghostedCoarseNodeCoarseIndices(3), ghostedCoarseNodeFineIndices(3), ijk(3);
    LO currentIndex = -1, coarseNodeFineLID = -1, computedCoarseNode = -1;
    for(ijk[2] = 0; ijk[2] < geoData->ghostedCoarseNodesPerDir[2]; ++ijk[2]) {
      for(ijk[1] = 0; ijk[1] < geoData->ghostedCoarseNodesPerDir[1]; ++ijk[1]) {
        for(ijk[0] = 0; ijk[0] < geoData->ghostedCoarseNodesPerDir[0]; ++ijk[0]) {
          currentIndex = (ijk[2]*geoData->ghostedCoarseNodesPerDir[1]
                          *geoData->ghostedCoarseNodesPerDir[0])
            + ijk[1]*geoData->ghostedCoarseNodesPerDir[0]
            + ijk[0];
          ghostedCoarseNodeCoarseIndices[0] = geoData->startGhostedCoarseNode[0] + ijk[0];
          ghostedCoarseNodeCoarseIndices[1] = geoData->startGhostedCoarseNode[1] + ijk[1];
          ghostedCoarseNodeCoarseIndices[2] = geoData->startGhostedCoarseNode[2] + ijk[2];
          GO myCoarseGID = ghostedCoarseNodeCoarseIndices[0]
            + ghostedCoarseNodeCoarseIndices[1]*geoData->gCoarseNodesPerDir[0]
            + (ghostedCoarseNodeCoarseIndices[2]*geoData->gCoarseNodesPerDir[1]
               *geoData->gCoarseNodesPerDir[0]);
          ghostedCoarseNodeCoarseGIDs[currentIndex] = myCoarseGID;
          GO myGID = 0, factor[3] = {};
          factor[2] = geoData->gNumFineNodes10;
          factor[1] = geoData->gFineNodesPerDir[0];
          factor[0] = 1;
          for(int dim = 0; dim < 3; ++dim) {
            if(dim < geoData->numDimensions) {
              if(geoData->startIndices[dim] - geoData->offsets[dim]
                 + ijk[dim]*geoData->coarseRate[dim] < geoData->gFineNodesPerDir[dim] - 1) {
                myGID += (geoData->startIndices[dim] - geoData->offsets[dim]
                          + ijk[dim]*geoData->coarseRate[dim])*factor[dim];
              } else {
                myGID += (geoData->startIndices[dim] - geoData->offsets[dim] + (ijk[dim] - 1)
                          *geoData->coarseRate[dim] + geoData->endRate[dim])*factor[dim];
              }
            }
          }
          if((!geoData->ghostInterface[0] || ijk[0] != 0)
             && (!geoData->ghostInterface[2] || ijk[1] != 0)
             && (!geoData->ghostInterface[4] || ijk[2] != 0)
             && (!geoData->ghostInterface[1] || ijk[0] != geoData->ghostedCoarseNodesPerDir[0] - 1)
             && (!geoData->ghostInterface[3] || ijk[1] != geoData->ghostedCoarseNodesPerDir[1] - 1)
             && (!geoData->ghostInterface[5] || ijk[2] != geoData->ghostedCoarseNodesPerDir[2] - 1)) {
            geoData->getGhostedNodeFineLID(ijk[0], ijk[1], ijk[2], coarseNodeFineLID);
            geoData->getGhostedNodeCoarseLID(ijk[0], ijk[1], ijk[2], computedCoarseNode);

            aggregates->SetIsRoot(coarseNodeFineLID);
            lCoarseNodeCoarseGIDs[computedCoarseNode] = myCoarseGID;
            lCoarseNodeFineGIDs[computedCoarseNode]   = myGID;
          }
          ghostedCoarseNodeFineGIDs[currentIndex] = myGID;
        }
      }
    }

    RCP<const Map> coarseCoordMap = MapFactory::Build (coordMap->lib(),
                                                       geoData->gNumCoarseNodes,
                                                       lCoarseNodeCoarseGIDs(),
                                                       coordMap->getIndexBase(),
                                                       coordMap->getComm());


    Array<int> ghostedCoarseNodeCoarsePIDs(geoData->numGhostedCoarseNodes);
    Array<LO>  ghostedCoarseNodeCoarseLIDs(geoData->numGhostedCoarseNodes);
    coarseCoordMap->getRemoteIndexList(ghostedCoarseNodeCoarseGIDs(),
                                       ghostedCoarseNodeCoarsePIDs(),
                                       ghostedCoarseNodeCoarseLIDs());

    // Now we are ready for the big loop over the fine node that will assign each
    // node on the fine grid to an aggregate and a processor.
    ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()  ->getDataNonConst(0);
    LO iGhosted, jGhosted, kGhosted, iCoarse, jCoarse, kCoarse, iRem, jRem, kRem;
    LO ghostedCoarseNodeCoarseLID, aggId;
    for(LO nodeIdx = 0; nodeIdx < geoData->lNumFineNodes; ++nodeIdx) {
      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, iGhosted, jGhosted, kGhosted);
      iCoarse = iGhosted / geoData->coarseRate[0];
      iRem    = iGhosted % geoData->coarseRate[0];
      if(iRem > (geoData->coarseRate[0] / 2)) { ++iCoarse; }
      jCoarse = jGhosted / geoData->coarseRate[1];
      jRem    = jGhosted % geoData->coarseRate[1];
      if(jRem > (geoData->coarseRate[1] / 2)) { ++jCoarse; }
      kCoarse = kGhosted / geoData->coarseRate[2];
      kRem    = kGhosted % geoData->coarseRate[2];
      if(kRem > (geoData->coarseRate[2] / 2)) { ++kCoarse; }
      geoData->getCoarseNodeGhostedLID(iCoarse, jCoarse, kCoarse, ghostedCoarseNodeCoarseLID);

      aggId                 = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID];
      vertex2AggId[nodeIdx] = aggId;
      procWinner[nodeIdx]   = ghostedCoarseNodeCoarsePIDs[ghostedCoarseNodeCoarseLID];
      aggStat[nodeIdx]      = AGGREGATED;
      --numNonAggregatedNodes;
    }

  } // GlobalLexicographicLayout

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  LocalLexicographicLayout(const RCP<const Map> coordMap, RCP<GeometricData> geoData,
                           RCP<Aggregates> aggregates, std::vector<unsigned>& aggStat,
                           LO& numNonAggregatedNodes) const {

    aggregates->SetNumAggregates(geoData->lNumCoarseNodes);

    // Now the tricky part starts, the coarse nodes / ghosted coarse nodes need to be imported.
    // This requires finding what their GID on the fine mesh is. They need to be ordered
    // lexicographically to allow for fast sweeps through the mesh.

    // We loop over all ghosted coarse nodes by increasing global lexicographic order
    Array<LO>  ghostedCoarseNodeCoarseIndices(3), ghostedCoarseNodeFineIndices(3);
    Array<LO>  lCoarseNodeCoarseIndices(3);
    Array<GO>  lCoarseNodeCoarseGIDs(geoData->lNumCoarseNodes);
    Array<LO>  ghostedCoarseNodeCoarseLIDs(geoData->numGhostedCoarseNodes);
    Array<int> ghostedCoarseNodeCoarsePIDs(geoData->numGhostedCoarseNodes);
    LO currentIndex = -1, countCoarseNodes = 0;
    geoData->computeCoarseLocalLexicographicData();
    for(int k = 0; k < geoData->ghostedCoarseNodesPerDir[2]; ++k) {
      for(int j = 0; j < geoData->ghostedCoarseNodesPerDir[1]; ++j) {
        for(int i = 0; i < geoData->ghostedCoarseNodesPerDir[0]; ++i) {
          currentIndex = k*geoData->ghostedCoarseNodesPerDir[1]*geoData->ghostedCoarseNodesPerDir[0]
            + j*geoData->ghostedCoarseNodesPerDir[0] + i;
          ghostedCoarseNodeCoarseIndices[0] = geoData->startGhostedCoarseNode[0] + i;
          ghostedCoarseNodeFineIndices[0] = ghostedCoarseNodeCoarseIndices[0]*geoData->coarseRate[0];
          if(ghostedCoarseNodeFineIndices[0] > geoData->gFineNodesPerDir[0] - 1) {
            ghostedCoarseNodeFineIndices[0] = geoData->gFineNodesPerDir[0] - 1;
          }
          ghostedCoarseNodeCoarseIndices[1] = geoData->startGhostedCoarseNode[1] + j;
          ghostedCoarseNodeFineIndices[1] = ghostedCoarseNodeCoarseIndices[1]*geoData->coarseRate[1];
          if(ghostedCoarseNodeFineIndices[1] > geoData->gFineNodesPerDir[1] - 1) {
            ghostedCoarseNodeFineIndices[1] = geoData->gFineNodesPerDir[1] - 1;
          }
          ghostedCoarseNodeCoarseIndices[2] = geoData->startGhostedCoarseNode[2] + k;
          ghostedCoarseNodeFineIndices[2] = ghostedCoarseNodeCoarseIndices[2]*geoData->coarseRate[2];
          if(ghostedCoarseNodeFineIndices[2] > geoData->gFineNodesPerDir[2] - 1) {
            ghostedCoarseNodeFineIndices[2] = geoData->gFineNodesPerDir[2] - 1;
          }

          GO myGID = -1, myCoarseGID = -1;
          LO myLID = -1, myPID = -1, myCoarseLID = -1;
          geoData->GetGIDLocalLexicographic(i, j, k, ghostedCoarseNodeFineIndices,
                                            myGID, myPID, myLID);

          int myRankIndex = geoData->rankIndices[myPID];
          for(int dim = 0; dim < 3; ++dim) {
            if(dim < geoData->numDimensions) {
              lCoarseNodeCoarseIndices[dim] = ghostedCoarseNodeCoarseIndices[dim]
                - geoData->coarseMeshData[myRankIndex][3 + 2*dim];
            }
          }
          LO myRankIndexCoarseNodesInDir0 = geoData->coarseMeshData[myRankIndex][4]
            - geoData->coarseMeshData[myRankIndex][3] + 1;
          LO myRankIndexCoarseNodes10 = (geoData->coarseMeshData[myRankIndex][6]
                                         - geoData->coarseMeshData[myRankIndex][5] + 1)
            *myRankIndexCoarseNodesInDir0;
          myCoarseLID = lCoarseNodeCoarseIndices[2]*myRankIndexCoarseNodes10
            + lCoarseNodeCoarseIndices[1]*myRankIndexCoarseNodesInDir0
            + lCoarseNodeCoarseIndices[0];
          myCoarseGID = myCoarseLID + geoData->coarseMeshData[myRankIndex][9];

          ghostedCoarseNodeCoarseLIDs[currentIndex] = myCoarseLID;
          ghostedCoarseNodeCoarsePIDs[currentIndex] = myPID;

          if(myPID == geoData->myRank) {
            lCoarseNodeCoarseGIDs[countCoarseNodes] = myCoarseGID;
            ++countCoarseNodes;
          }
        }
      }
    }

    RCP<const Map> coarseCoordMap = MapFactory::Build (coordMap->lib(),
                                                       geoData->gNumCoarseNodes,
                                                       lCoarseNodeCoarseGIDs(),
                                                       coordMap->getIndexBase(),
                                                       coordMap->getComm());

    // Now we are ready for the big loop over the fine node that will assign each
    // node on the fine grid to an aggregate and a processor.
    ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()  ->getDataNonConst(0);
    LO iGhosted, jGhosted, kGhosted, iCoarse, jCoarse, kCoarse, iRem, jRem, kRem;
    LO ghostedCoarseNodeCoarseLID, aggId;
    for(LO nodeIdx = 0; nodeIdx < geoData->lNumFineNodes; ++nodeIdx) {
      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, iGhosted, jGhosted, kGhosted);
      iCoarse = iGhosted / geoData->coarseRate[0];
      iRem    = iGhosted % geoData->coarseRate[0];
      if(iRem > (geoData->coarseRate[0] / 2)) { ++iCoarse; }
      jCoarse = jGhosted / geoData->coarseRate[1];
      jRem    = jGhosted % geoData->coarseRate[1];
      if(jRem > (geoData->coarseRate[1] / 2)) { ++jCoarse; }
      kCoarse = kGhosted / geoData->coarseRate[2];
      kRem    = kGhosted % geoData->coarseRate[2];
      if(kRem > (geoData->coarseRate[2] / 2)) { ++kCoarse; }
      geoData->getCoarseNodeGhostedLID(iCoarse, jCoarse, kCoarse, ghostedCoarseNodeCoarseLID);

      aggId                 = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID];
      vertex2AggId[nodeIdx] = aggId;
      procWinner[nodeIdx]   = ghostedCoarseNodeCoarsePIDs[ghostedCoarseNodeCoarseLID];
      aggStat[nodeIdx]      = AGGREGATED;
      --numNonAggregatedNodes;
    }

  } // LocalLexicographicLayout

} //namespace MueLu


#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_ */
