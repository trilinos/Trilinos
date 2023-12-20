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

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_AggregationStructuredAlgorithm.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_UncoupledIndexManager.hpp"
#include "MueLu_LocalLexicographicIndexManager.hpp"
#include "MueLu_GlobalLexicographicIndexManager.hpp"

#include "MueLu_StructuredAggregationFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    StructuredAggregationFactory()
  : bDefinitionPhase_(true) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
  SET_VALID_ENTRY("aggregation: allow user-specified singletons");
  SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
  SET_VALID_ENTRY("aggregation: phase3 avoid singletons");

  // general variables needed in StructuredAggregationFactory
  SET_VALID_ENTRY("aggregation: mesh layout");
  SET_VALID_ENTRY("aggregation: mode");
  SET_VALID_ENTRY("aggregation: output type");
  SET_VALID_ENTRY("aggregation: coarsening rate");
  SET_VALID_ENTRY("aggregation: coarsening order");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("Graph", Teuchos::null,
                                               "Graph of the matrix after amalgamation but without dropping.");
  validParamList->set<RCP<const FactoryBase> >("numDimensions", Teuchos::null,
                                               "Number of spatial dimension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("gNodesPerDim", Teuchos::null,
                                               "Global number of nodes per spatial dimension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null,
                                               "Local number of nodes per spatial dimension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", Teuchos::null,
                                               "Generating factory for variable \'DofsPerNode\', usually the same as the \'Graph\' factory");
  validParamList->set<const bool>("aggregation: single coarse point", false,
                                  "Allows the aggreagtion process to reduce spacial dimensions to a single layer");

  return validParamList;
}  // GetValidParameterList()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");

  ParameterList pL     = GetParameterList();
  std::string coupling = pL.get<std::string>("aggregation: mode");
  const bool coupled   = (coupling == "coupled" ? true : false);
  if (coupled) {
    // Request the global number of nodes per dimensions
    if (currentLevel.GetLevelID() == 0) {
      if (currentLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
        currentLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "gNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "gNodesPerDim");
    }
  }

  // Request the local number of nodes per dimensions
  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("numDimensions", NoFactory::get())) {
      currentLevel.DeclareInput("numDimensions", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("numDimensions", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "numDimensions was not provided by the user on level0!");
    }
    if (currentLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
      currentLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "lNodesPerDim was not provided by the user on level0!");
    }
  } else {
    Input(currentLevel, "numDimensions");
    Input(currentLevel, "lNodesPerDim");
  }
}  // DeclareInput()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_STRUCTUREDAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  *out << "Entering structured aggregation" << std::endl;

  ParameterList pL  = GetParameterList();
  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  // General problem informations are gathered from data stored in the problem matix.
  RCP<const GraphBase> graph = Get<RCP<GraphBase> >(currentLevel, "Graph");
  RCP<const Map> fineMap     = graph->GetDomainMap();
  const int myRank           = fineMap->getComm()->getRank();
  const int numRanks         = fineMap->getComm()->getSize();
  const GO minGlobalIndex    = fineMap->getMinGlobalIndex();
  const LO dofsPerNode       = Get<LO>(currentLevel, "DofsPerNode");

  // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
  // obtain a nodeMap.
  const int interpolationOrder = pL.get<int>("aggregation: coarsening order");
  std::string meshLayout       = pL.get<std::string>("aggregation: mesh layout");
  std::string coupling         = pL.get<std::string>("aggregation: mode");
  const bool coupled           = (coupling == "coupled" ? true : false);
  std::string outputType       = pL.get<std::string>("aggregation: output type");
  const bool outputAggregates  = (outputType == "Aggregates" ? true : false);
  const bool singleCoarsePoint = pL.get<bool>("aggregation: single coarse point");
  int numDimensions;
  Array<GO> gFineNodesPerDir(3);
  Array<LO> lFineNodesPerDir(3);
  if (currentLevel.GetLevelID() == 0) {
    // On level 0, data is provided by applications and has no associated factory.
    numDimensions    = currentLevel.Get<int>("numDimensions", NoFactory::get());
    lFineNodesPerDir = currentLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    if (coupled) {
      gFineNodesPerDir = currentLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
    }
  } else {
    // On level > 0, data is provided directly by generating factories.
    numDimensions    = Get<int>(currentLevel, "numDimensions");
    lFineNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
    if (coupled) {
      gFineNodesPerDir = Get<Array<GO> >(currentLevel, "gNodesPerDim");
    }
  }

  // First make sure that input parameters are set logically based on dimension
  for (int dim = 0; dim < 3; ++dim) {
    if (dim >= numDimensions) {
      gFineNodesPerDir[dim] = 1;
      lFineNodesPerDir[dim] = 1;
    }
  }

  // Get the coarsening rate
  std::string coarseningRate = pL.get<std::string>("aggregation: coarsening rate");
  Teuchos::Array<LO> coarseRate;
  try {
    coarseRate = Teuchos::fromStringToArray<LO>(coarseningRate);
  } catch (const Teuchos::InvalidArrayStringRepresentation& e) {
    GetOStream(Errors, -1) << " *** \"aggregation: coarsening rate\" must be a string convertible into an array! *** "
                           << std::endl;
    throw e;
  }
  TEUCHOS_TEST_FOR_EXCEPTION((coarseRate.size() > 1) && (coarseRate.size() < numDimensions),
                             Exceptions::RuntimeError,
                             "\"aggregation: coarsening rate\" must have at least as many"
                             " components as the number of spatial dimensions in the problem.");

  // Now that we have extracted info from the level, create the IndexManager
  RCP<IndexManager> geoData;
  if (!coupled) {
    geoData = rcp(new MueLu::UncoupledIndexManager<LO, GO, NO>(fineMap->getComm(),
                                                               coupled,
                                                               numDimensions,
                                                               interpolationOrder,
                                                               myRank,
                                                               numRanks,
                                                               gFineNodesPerDir,
                                                               lFineNodesPerDir,
                                                               coarseRate,
                                                               singleCoarsePoint));
  } else if (meshLayout == "Local Lexicographic") {
    Array<GO> meshData;
    if (currentLevel.GetLevelID() == 0) {
      // On level 0, data is provided by applications and has no associated factory.
      meshData = currentLevel.Get<Array<GO> >("aggregation: mesh data", NoFactory::get());
      TEUCHOS_TEST_FOR_EXCEPTION(meshData.empty() == true, Exceptions::RuntimeError,
                                 "The meshData array is empty, somehow the input for structured"
                                 " aggregation are not captured correctly.");
    } else {
      // On level > 0, data is provided directly by generating factories.
      meshData = Get<Array<GO> >(currentLevel, "aggregation: mesh data");
    }
    // Note, LBV Feb 5th 2018:
    // I think that it might make sense to pass ghostInterface rather than interpolationOrder.
    // For that I need to make sure that ghostInterface can be computed with minimal mesh
    // knowledge outside of the IndexManager...
    geoData = rcp(new MueLu::LocalLexicographicIndexManager<LO, GO, NO>(fineMap->getComm(),
                                                                        coupled,
                                                                        numDimensions,
                                                                        interpolationOrder,
                                                                        myRank,
                                                                        numRanks,
                                                                        gFineNodesPerDir,
                                                                        lFineNodesPerDir,
                                                                        coarseRate,
                                                                        meshData));
  } else if (meshLayout == "Global Lexicographic") {
    // Note, LBV Feb 5th 2018:
    // I think that it might make sense to pass ghostInterface rather than interpolationOrder.
    // For that I need to make sure that ghostInterface can be computed with minimal mesh
    // knowledge outside of the IndexManager...
    geoData = rcp(new MueLu::GlobalLexicographicIndexManager<LO, GO, NO>(fineMap->getComm(),
                                                                         coupled,
                                                                         numDimensions,
                                                                         interpolationOrder,
                                                                         gFineNodesPerDir,
                                                                         lFineNodesPerDir,
                                                                         coarseRate,
                                                                         minGlobalIndex));
  }

  *out << "The index manager has now been built" << std::endl;
  *out << "graph num nodes: " << fineMap->getLocalNumElements()
       << ", structured aggregation num nodes: " << geoData->getNumLocalFineNodes() << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getLocalNumElements() != static_cast<size_t>(geoData->getNumLocalFineNodes()),
                             Exceptions::RuntimeError,
                             "The local number of elements in the graph's map is not equal to "
                             "the number of nodes given by: lNodesPerDim!");
  if (coupled) {
    TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getGlobalNumElements() != static_cast<size_t>(geoData->getNumGlobalFineNodes()),
                               Exceptions::RuntimeError,
                               "The global number of elements in the graph's map is not equal to "
                               "the number of nodes given by: gNodesPerDim!");
  }

  *out << "Compute coarse mesh data" << std::endl;
  std::vector<std::vector<GO> > coarseMeshData = geoData->getCoarseMeshData();

  // Now we are ready for the big loop over the fine node that will assign each
  // node on the fine grid to an aggregate and a processor.
  RCP<const FactoryBase> graphFact = GetFactory("Graph");
  RCP<const Map> coarseCoordinatesFineMap, coarseCoordinatesMap;
  RCP<MueLu::AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node> >
      myStructuredAlgorithm = rcp(new AggregationStructuredAlgorithm(graphFact));

  if (interpolationOrder == 0 && outputAggregates) {
    // Create aggregates for prolongation
    *out << "Compute Aggregates" << std::endl;
    RCP<Aggregates> aggregates = rcp(new Aggregates(graph->GetDomainMap()));
    aggregates->setObjectLabel("ST");
    aggregates->SetIndexManager(geoData);
    aggregates->AggregatesCrossProcessors(coupled);
    aggregates->SetNumAggregates(geoData->getNumLocalCoarseNodes());
    std::vector<unsigned> aggStat(geoData->getNumLocalFineNodes(), READY);
    LO numNonAggregatedNodes = geoData->getNumLocalFineNodes();

    myStructuredAlgorithm->BuildAggregates(pL, *graph, *aggregates, aggStat,
                                           numNonAggregatedNodes);

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError,
                               "MueLu::StructuredAggregationFactory::Build: Leftover nodes found! Error!");
    aggregates->ComputeAggregateSizes(true /*forceRecompute*/);
    GetOStream(Statistics1) << aggregates->description() << std::endl;
    Set(currentLevel, "Aggregates", aggregates);

  } else {
    // Create the graph of the prolongator
    *out << "Compute CrsGraph" << std::endl;
    RCP<CrsGraph> myGraph;
    myStructuredAlgorithm->BuildGraph(*graph, geoData, dofsPerNode, myGraph,
                                      coarseCoordinatesFineMap, coarseCoordinatesMap);
    Set(currentLevel, "prolongatorGraph", myGraph);
  }

  if (coupled) {
    Set(currentLevel, "gCoarseNodesPerDim", geoData->getGlobalCoarseNodesPerDir());
  }
  Set(currentLevel, "lCoarseNodesPerDim", geoData->getLocalCoarseNodesPerDir());
  Set(currentLevel, "coarseCoordinatesFineMap", coarseCoordinatesFineMap);
  Set(currentLevel, "coarseCoordinatesMap", coarseCoordinatesMap);
  Set(currentLevel, "structuredInterpolationOrder", interpolationOrder);
  Set(currentLevel, "numDimensions", numDimensions);

}  // Build()
}  // namespace MueLu

#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_ */
