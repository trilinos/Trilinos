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
#include "MueLu_LocalLexicographicIndexManager.hpp"
#include "MueLu_GlobalLexicographicIndexManager.hpp"


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
    const GO  minGlobalIndex    = Coordinates->getMap()->getMinGlobalIndex();

    // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
    // obtain a nodeMap.
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
    TEUCHOS_TEST_FOR_EXCEPTION((coarseRate.size() > 1) && (coarseRate.size() < numDimensions),
                               Exceptions::RuntimeError,
                               "\"aggregation: coarsening rate\" must have at least as many"
                               " components as the number of spatial dimensions in the problem.");

    // Now that we have extract info from the level, create the IndexManager
    RCP<MueLu::IndexManager<LO,GO,NO> > geoData;
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
      geoData = rcp(new MueLu::LocalLexicographicIndexManager<LO,GO,NO>(numDimensions, myRank,
                                                                        numRanks, gFineNodesPerDir,
                                                                        lFineNodesPerDir, coarseRate,
                                                                        meshData));
    } else  if(meshLayout == "Global Lexicographic") {
      geoData = rcp(new MueLu::GlobalLexicographicIndexManager<LO,GO,NO>(numDimensions,
                                                                         gFineNodesPerDir,
                                                                         lFineNodesPerDir,
                                                                         coarseRate,
                                                                         minGlobalIndex));
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Coordinates->getLocalLength()
                               != static_cast<size_t>(geoData->getNumLocalFineNodes()),
                               Exceptions::RuntimeError,
                               "The local number of elements in Coordinates is not equal to the"
                               " number of nodes given by: lNodesPerDim!");
    TEUCHOS_TEST_FOR_EXCEPTION(Coordinates->getGlobalLength()
                               != static_cast<size_t>(geoData->getNumGlobalFineNodes()),
                               Exceptions::RuntimeError,
                               "The global number of elements in Coordinates is not equal to the"
                               " number of nodes given by: gNodesPerDim!");

    RCP<const Map> coarseCoordMap;
    Array<LO>  ghostedCoarseNodeCoarseLIDs(geoData->getNumLocalGhostedNodes());
    Array<int> ghostedCoarseNodeCoarsePIDs(geoData->getNumLocalGhostedNodes());
    geoData->getGhostedNodesData(coordMap, coarseCoordMap, ghostedCoarseNodeCoarseLIDs,
                                 ghostedCoarseNodeCoarsePIDs);

    // Create aggregates object and set basic parameters
    RCP<Aggregates> aggregates = rcp(new Aggregates(coordMap));
    aggregates->setObjectLabel("ST");
    aggregates->AggregatesCrossProcessors(true);
    std::vector<unsigned> aggStat(geoData->getNumLocalFineNodes(), READY);
    aggregates->SetNumAggregates(geoData->getNumLocalCoarseNodes());


    // Now we are ready for the big loop over the fine node that will assign each
    // node on the fine grid to an aggregate and a processor.
    LO numNonAggregatedNodes = geoData->getNumLocalFineNodes();
    ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates->GetProcWinner()  ->getDataNonConst(0);
    LO iGhosted, jGhosted, kGhosted, iCoarse, jCoarse, kCoarse, iRem, jRem, kRem;
    LO ghostedCoarseNodeCoarseLID, aggId;
    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {
      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, iGhosted, jGhosted, kGhosted);
      iCoarse = iGhosted / geoData->getCoarseningRate(0);
      iRem    = iGhosted % geoData->getCoarseningRate(0);
      if(iRem > (geoData->getCoarseningRate(0) / 2)) { ++iCoarse; }
      jCoarse = jGhosted / geoData->getCoarseningRate(1);
      jRem    = jGhosted % geoData->getCoarseningRate(1);
      if(jRem > (geoData->getCoarseningRate(1) / 2)) { ++jCoarse; }
      kCoarse = kGhosted / geoData->getCoarseningRate(2);
      kRem    = kGhosted % geoData->getCoarseningRate(2);
      if(kRem > (geoData->getCoarseningRate(2) / 2)) { ++kCoarse; }
      geoData->getCoarseNodeGhostedLID(iCoarse, jCoarse, kCoarse, ghostedCoarseNodeCoarseLID);

      aggId                 = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID];
      vertex2AggId[nodeIdx] = aggId;
      procWinner[nodeIdx]   = ghostedCoarseNodeCoarsePIDs[ghostedCoarseNodeCoarseLID];
      aggStat[nodeIdx]      = AGGREGATED;
      --numNonAggregatedNodes;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError,
                               "MueLu::StructuredAggregationFactory::Build: Leftover nodes found! Error!");

    // aggregates->AggregatesCrossProcessors(false);
    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics1) << aggregates->description() << std::endl;
  }

} //namespace MueLu


#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_ */
