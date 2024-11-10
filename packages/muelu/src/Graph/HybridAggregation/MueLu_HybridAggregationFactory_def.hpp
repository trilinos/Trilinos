// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HYBRIDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_HYBRIDAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_HybridAggregationFactory_decl.hpp"

// Uncoupled Agg
#include "MueLu_InterfaceAggregationAlgorithm.hpp"
#include "MueLu_OnePtAggregationAlgorithm.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm.hpp"

#include "MueLu_AggregationPhase1Algorithm.hpp"
#include "MueLu_AggregationPhase2aAlgorithm.hpp"
#include "MueLu_AggregationPhase2bAlgorithm.hpp"
#include "MueLu_AggregationPhase3Algorithm.hpp"

// Structured Agg
#include "MueLu_AggregationStructuredAlgorithm.hpp"
#include "MueLu_UncoupledIndexManager.hpp"
//#include "MueLu_LocalLexicographicIndexManager.hpp"
//#include "MueLu_GlobalLexicographicIndexManager.hpp"

// Shared
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    HybridAggregationFactory()
  : bDefinitionPhase_(true) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  // From UncoupledAggregationFactory
  SET_VALID_ENTRY("aggregation: max agg size");
  SET_VALID_ENTRY("aggregation: min agg size");
  SET_VALID_ENTRY("aggregation: max selected neighbors");
  SET_VALID_ENTRY("aggregation: ordering");
  validParamList->getEntry("aggregation: ordering").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("natural", "graph", "random"))));
  SET_VALID_ENTRY("aggregation: enable phase 1");
  SET_VALID_ENTRY("aggregation: enable phase 2a");
  SET_VALID_ENTRY("aggregation: enable phase 2b");
  SET_VALID_ENTRY("aggregation: enable phase 3");
  SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
  SET_VALID_ENTRY("aggregation: allow user-specified singletons");
  SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
  SET_VALID_ENTRY("aggregation: match ML phase1");
  SET_VALID_ENTRY("aggregation: match ML phase2a");
  SET_VALID_ENTRY("aggregation: match ML phase2b");
  SET_VALID_ENTRY("aggregation: phase2a agg factor");
  SET_VALID_ENTRY("aggregation: phase3 avoid singletons");

  // From StructuredAggregationFactory
  SET_VALID_ENTRY("aggregation: coarsening rate");
  SET_VALID_ENTRY("aggregation: coarsening order");
  SET_VALID_ENTRY("aggregation: number of spatial dimensions");

  // From HybridAggregationFactory
  SET_VALID_ENTRY("aggregation: use interface aggregation");
#undef SET_VALID_ENTRY

  /* From UncoupledAggregation */
  // general variables needed in AggregationFactory
  validParamList->set<RCP<const FactoryBase> >("Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
  // special variables necessary for OnePtAggregationAlgorithm
  validParamList->set<std::string>("OnePt aggregate map name", "",
                                   "Name of input map for single node aggregates. (default='')");
  validParamList->set<std::string>("OnePt aggregate map factory", "",
                                   "Generating factory of (DOF) map for single node aggregates.");

  // InterfaceAggregation parameters
  validParamList->set<std::string>("Interface aggregate map name", "",
                                   "Name of input map for interface aggregates. (default='')");
  validParamList->set<std::string>("Interface aggregate map factory", "",
                                   "Generating factory of (DOF) map for interface aggregates.");
  validParamList->set<RCP<const FactoryBase> >("interfacesDimensions", Teuchos::null,
                                               "Describes the dimensions of all the interfaces on this rank.");
  validParamList->set<RCP<const FactoryBase> >("nodeOnInterface", Teuchos::null,
                                               "List the LIDs of the nodes on any interface.");

  /* From StructuredAggregation */
  // general variables needed in AggregationFactory
  validParamList->set<RCP<const FactoryBase> >("numDimensions", Teuchos::null,
                                               "Number of spatial dimension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null,
                                               "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");

  // Hybrid Aggregation Params
  validParamList->set<RCP<const FactoryBase> >("aggregationRegionType", Teuchos::null,
                                               "Type of aggregation to use on the region (\"structured\" or \"uncoupled\")");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "Graph");

  ParameterList pL = GetParameterList();

  /* StructuredAggregation */

  // Request the local number of nodes per dimensions
  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("aggregationRegionType", NoFactory::get())) {
      currentLevel.DeclareInput("aggregationRegionType", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("aggregationRegionType", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "Aggregation region type was not provided by the user!");
    }
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
    Input(currentLevel, "aggregationRegionType");
    Input(currentLevel, "numDimensions");
    Input(currentLevel, "lNodesPerDim");
  }

  /* UncoupledAggregation */
  Input(currentLevel, "DofsPerNode");

  // request special data necessary for InterfaceAggregation
  if (pL.get<bool>("aggregation: use interface aggregation") == true) {
    if (currentLevel.GetLevelID() == 0) {
      if (currentLevel.IsAvailable("interfacesDimensions", NoFactory::get())) {
        currentLevel.DeclareInput("interfacesDimensions", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("interfacesDimensions", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "interfacesDimensions was not provided by the user on level0!");
      }
      if (currentLevel.IsAvailable("nodeOnInterface", NoFactory::get())) {
        currentLevel.DeclareInput("nodeOnInterface", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("nodeOnInterface", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "nodeOnInterface was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "interfacesDimensions");
      Input(currentLevel, "nodeOnInterface");
    }
  }

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
}  // DeclareInput()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_HYBRIDAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  *out << "Entering hybrid aggregation" << std::endl;

  ParameterList pL  = GetParameterList();
  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  if (pL.get<int>("aggregation: max agg size") == -1)
    pL.set("aggregation: max agg size", INT_MAX);

  // define aggregation algorithms
  RCP<const FactoryBase> graphFact = GetFactory("Graph");

  // General problem informations are gathered from data stored in the problem matix.
  RCP<const LWGraph> graph = Get<RCP<LWGraph> >(currentLevel, "Graph");
  RCP<const Map> fineMap   = graph->GetDomainMap();
  const int myRank         = fineMap->getComm()->getRank();
  const int numRanks       = fineMap->getComm()->getSize();

  out->setProcRankAndSize(graph->GetImportMap()->getComm()->getRank(),
                          graph->GetImportMap()->getComm()->getSize());

  // Build aggregates
  RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
  aggregates->setObjectLabel("HB");

  // construct aggStat information
  const LO numRows      = graph->GetNodeNumVertices();
  using AggStatHostType = typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType;
  AggStatHostType aggStat(Kokkos::ViewAllocateWithoutInitializing("aggregation status"), numRows);
  Kokkos::deep_copy(aggStat, READY);

  // Get aggregation type for region
  std::string regionType;
  if (currentLevel.GetLevelID() == 0) {
    // On level 0, data is provided by applications and has no associated factory.
    regionType = currentLevel.Get<std::string>("aggregationRegionType", NoFactory::get());
  } else {
    // On level > 0, data is provided directly by generating factories.
    regionType = Get<std::string>(currentLevel, "aggregationRegionType");
  }

  int numDimensions = 0;
  if (currentLevel.GetLevelID() == 0) {
    // On level 0, data is provided by applications and has no associated factory.
    numDimensions = currentLevel.Get<int>("numDimensions", NoFactory::get());
  } else {
    // On level > 0, data is provided directly by generating factories.
    numDimensions = Get<int>(currentLevel, "numDimensions");
  }

  // Get the coarsening rate (potentially used for both structured and uncoupled aggregation if interface)
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

  algos_.clear();
  LO numNonAggregatedNodes = numRows;
  if (regionType == "structured") {
    // Add AggregationStructuredAlgorithm
    algos_.push_back(rcp(new AggregationStructuredAlgorithm(graphFact)));

    // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
    // obtain a nodeMap.
    const int interpolationOrder = pL.get<int>("aggregation: coarsening order");
    Array<LO> lFineNodesPerDir(3);
    if (currentLevel.GetLevelID() == 0) {
      // On level 0, data is provided by applications and has no associated factory.
      lFineNodesPerDir = currentLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    } else {
      // On level > 0, data is provided directly by generating factories.
      lFineNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
    }

    // Set lFineNodesPerDir to 1 for directions beyond numDimensions
    for (int dim = numDimensions; dim < 3; ++dim) {
      lFineNodesPerDir[dim] = 1;
    }

    // Now that we have extracted info from the level, create the IndexManager
    RCP<MueLu::IndexManager<LO, GO, NO> > geoData;
    geoData = rcp(new MueLu::UncoupledIndexManager<LO, GO, NO>(fineMap->getComm(),
                                                               false,
                                                               numDimensions,
                                                               interpolationOrder,
                                                               myRank,
                                                               numRanks,
                                                               Array<GO>(3, -1),
                                                               lFineNodesPerDir,
                                                               coarseRate, false));

    TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getLocalNumElements() != static_cast<size_t>(geoData->getNumLocalFineNodes()),
                               Exceptions::RuntimeError,
                               "The local number of elements in the graph's map is not equal to "
                               "the number of nodes given by: lNodesPerDim!");

    aggregates->SetIndexManager(geoData);
    aggregates->SetNumAggregates(geoData->getNumLocalCoarseNodes());

    Set(currentLevel, "lCoarseNodesPerDim", geoData->getLocalCoarseNodesPerDir());

  }  // end structured aggregation setup

  if (regionType == "uncoupled") {
    // Add unstructred aggregation phases
    algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm(graphFact)));
    if (pL.get<bool>("aggregation: use interface aggregation") == true) algos_.push_back(rcp(new InterfaceAggregationAlgorithm(graphFact)));
    if (pL.get<bool>("aggregation: allow user-specified singletons") == true) algos_.push_back(rcp(new OnePtAggregationAlgorithm(graphFact)));
    if (pL.get<bool>("aggregation: enable phase 1") == true) algos_.push_back(rcp(new AggregationPhase1Algorithm(graphFact)));
    if (pL.get<bool>("aggregation: enable phase 2a") == true) algos_.push_back(rcp(new AggregationPhase2aAlgorithm(graphFact)));
    if (pL.get<bool>("aggregation: enable phase 2b") == true) algos_.push_back(rcp(new AggregationPhase2bAlgorithm(graphFact)));
    if (pL.get<bool>("aggregation: enable phase 3") == true) algos_.push_back(rcp(new AggregationPhase3Algorithm(graphFact)));

    *out << " Build interface aggregates" << std::endl;
    // interface
    if (pL.get<bool>("aggregation: use interface aggregation") == true) {
      BuildInterfaceAggregates(currentLevel, aggregates, aggStat, numNonAggregatedNodes,
                               coarseRate);
    }

    *out << "Treat Dirichlet BC" << std::endl;
    // Dirichlet boundary
    auto dirichletBoundaryMap = graph->GetBoundaryNodeMap();
    for (LO i = 0; i < numRows; i++)
      if (dirichletBoundaryMap[i] == true)
        aggStat[i] = BOUNDARY;

    // OnePt aggregation
    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
    RCP<Map> OnePtMap        = Teuchos::null;
    if (mapOnePtName.length()) {
      std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
      if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
        OnePtMap = currentLevel.Get<RCP<Map> >(mapOnePtName, NoFactory::get());
      } else {
        RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
        OnePtMap                            = currentLevel.Get<RCP<Map> >(mapOnePtName, mapOnePtFact.get());
      }
    }

    LO nDofsPerNode = Get<LO>(currentLevel, "DofsPerNode");
    GO indexBase    = graph->GetDomainMap()->getIndexBase();
    if (OnePtMap != Teuchos::null) {
      for (LO i = 0; i < numRows; i++) {
        // reconstruct global row id (FIXME only works for contiguous maps)
        GO grid = (graph->GetDomainMap()->getGlobalElement(i) - indexBase) * nDofsPerNode + indexBase;
        for (LO kr = 0; kr < nDofsPerNode; kr++)
          if (OnePtMap->isNodeGlobalElement(grid + kr))
            aggStat[i] = ONEPT;
      }
    }

    // Create a fake lCoarseNodesPerDir for CoordinatesTranferFactory
    Array<LO> lCoarseNodesPerDir(3, -1);
    Set(currentLevel, "lCoarseNodesPerDim", lCoarseNodesPerDir);
  }  // end uncoupled aggregation setup

  aggregates->AggregatesCrossProcessors(false);  // No coupled aggregation

  *out << "Run all the algorithms on the local rank" << std::endl;
  for (size_t a = 0; a < algos_.size(); a++) {
    std::string phase = algos_[a]->description();
    SubFactoryMonitor sfm(*this, "Algo \"" + phase + "\"", currentLevel);
    *out << regionType << " | Executing phase " << a << std::endl;

    int oldRank = algos_[a]->SetProcRankVerbose(this->GetProcRankVerbose());
    algos_[a]->BuildAggregatesNonKokkos(pL, *graph, *aggregates, aggStat, numNonAggregatedNodes);
    algos_[a]->SetProcRankVerbose(oldRank);
    *out << regionType << " | Done Executing phase " << a << std::endl;
  }

  *out << "Compute statistics on aggregates" << std::endl;
  aggregates->ComputeAggregateSizes(true /*forceRecompute*/);

  Set(currentLevel, "Aggregates", aggregates);
  Set(currentLevel, "numDimensions", numDimensions);
  Set(currentLevel, "aggregationRegionTypeCoarse", regionType);

  GetOStream(Statistics1) << aggregates->description() << std::endl;
  *out << "HybridAggregation done!" << std::endl;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    BuildInterfaceAggregates(Level& currentLevel, RCP<Aggregates> aggregates,
                             typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat, LO& numNonAggregatedNodes,
                             Array<LO> coarseRate) const {
  FactoryMonitor m(*this, "BuildInterfaceAggregates", currentLevel);

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_HYBRIDAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  // Extract and format input data for algo
  if (coarseRate.size() == 1) {
    coarseRate.resize(3, coarseRate[0]);
  }
  ArrayRCP<LO> vertex2AggId      = aggregates->GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner        = aggregates->GetProcWinner()->getDataNonConst(0);
  Array<LO> interfacesDimensions = Get<Array<LO> >(currentLevel, "interfacesDimensions");
  Array<LO> nodesOnInterfaces    = Get<Array<LO> >(currentLevel, "nodeOnInterface");
  const int numInterfaces        = interfacesDimensions.size() / 3;
  const int myRank               = aggregates->GetMap()->getComm()->getRank();

  // Create coarse level container to gather data on the fly
  Array<LO> coarseInterfacesDimensions(interfacesDimensions.size());
  Array<LO> nodesOnCoarseInterfaces;
  {  // Scoping the temporary variables...
    LO endRate, totalNumCoarseNodes = 0, numCoarseNodes;
    for (int interfaceIdx = 0; interfaceIdx < numInterfaces; ++interfaceIdx) {
      numCoarseNodes = 1;
      for (int dim = 0; dim < 3; ++dim) {
        endRate = (interfacesDimensions[3 * interfaceIdx + dim] - 1) % coarseRate[dim];
        if (interfacesDimensions[3 * interfaceIdx + dim] == 1) {
          coarseInterfacesDimensions[3 * interfaceIdx + dim] = 1;
        } else {
          coarseInterfacesDimensions[3 * interfaceIdx + dim] = (interfacesDimensions[3 * interfaceIdx + dim] - 1) / coarseRate[dim] + 2;
          if (endRate == 0) {
            coarseInterfacesDimensions[3 * interfaceIdx + dim]--;
          }
        }
        numCoarseNodes *= coarseInterfacesDimensions[3 * interfaceIdx + dim];
      }
      totalNumCoarseNodes += numCoarseNodes;
    }
    nodesOnCoarseInterfaces.resize(totalNumCoarseNodes, -1);
  }

  Array<LO> endRate(3);
  LO interfaceOffset = 0, aggregateCount = 0, coarseNodeCount = 0;
  for (int interfaceIdx = 0; interfaceIdx < numInterfaces; ++interfaceIdx) {
    ArrayView<LO> fineNodesPerDim   = interfacesDimensions(3 * interfaceIdx, 3);
    ArrayView<LO> coarseNodesPerDim = coarseInterfacesDimensions(3 * interfaceIdx, 3);
    LO numInterfaceNodes = 1, numCoarseNodes = 1;
    for (int dim = 0; dim < 3; ++dim) {
      numInterfaceNodes *= fineNodesPerDim[dim];
      numCoarseNodes *= coarseNodesPerDim[dim];
      endRate[dim] = (fineNodesPerDim[dim] - 1) % coarseRate[dim];
    }
    ArrayView<LO> interfaceNodes = nodesOnInterfaces(interfaceOffset, numInterfaceNodes);

    interfaceOffset += numInterfaceNodes;

    LO rem, rate, fineNodeIdx;
    Array<LO> nodeIJK(3), coarseIJK(3), rootIJK(3);
    // First find treat coarse nodes as they generate the aggregate IDs
    // and they might be repeated on multiple interfaces (think corners and edges).
    for (LO coarseNodeIdx = 0; coarseNodeIdx < numCoarseNodes; ++coarseNodeIdx) {
      coarseIJK[2] = coarseNodeIdx / (coarseNodesPerDim[0] * coarseNodesPerDim[1]);
      rem          = coarseNodeIdx % (coarseNodesPerDim[0] * coarseNodesPerDim[1]);
      coarseIJK[1] = rem / coarseNodesPerDim[0];
      coarseIJK[0] = rem % coarseNodesPerDim[0];

      for (LO dim = 0; dim < 3; ++dim) {
        if (coarseIJK[dim] == coarseNodesPerDim[dim] - 1) {
          nodeIJK[dim] = fineNodesPerDim[dim] - 1;
        } else {
          nodeIJK[dim] = coarseIJK[dim] * coarseRate[dim];
        }
      }
      fineNodeIdx = (nodeIJK[2] * fineNodesPerDim[1] + nodeIJK[1]) * fineNodesPerDim[0] + nodeIJK[0];

      if (aggStat[interfaceNodes[fineNodeIdx]] == READY) {
        vertex2AggId[interfaceNodes[fineNodeIdx]] = aggregateCount;
        procWinner[interfaceNodes[fineNodeIdx]]   = myRank;
        aggStat[interfaceNodes[fineNodeIdx]]      = AGGREGATED;
        ++aggregateCount;
        --numNonAggregatedNodes;
      }
      nodesOnCoarseInterfaces[coarseNodeCount] = vertex2AggId[interfaceNodes[fineNodeIdx]];
      ++coarseNodeCount;
    }

    // Now loop over all the node on the interface
    // skip the coarse nodes as they are already aggregated
    // and find the appropriate aggregate ID for the fine nodes.
    for (LO nodeIdx = 0; nodeIdx < numInterfaceNodes; ++nodeIdx) {
      // If the node is already aggregated skip it!
      if (aggStat[interfaceNodes[nodeIdx]] == AGGREGATED) {
        continue;
      }

      nodeIJK[2] = nodeIdx / (fineNodesPerDim[0] * fineNodesPerDim[1]);
      rem        = nodeIdx % (fineNodesPerDim[0] * fineNodesPerDim[1]);
      nodeIJK[1] = rem / fineNodesPerDim[0];
      nodeIJK[0] = rem % fineNodesPerDim[0];

      for (int dim = 0; dim < 3; ++dim) {
        coarseIJK[dim] = nodeIJK[dim] / coarseRate[dim];
        rem            = nodeIJK[dim] % coarseRate[dim];
        if (nodeIJK[dim] < fineNodesPerDim[dim] - endRate[dim]) {
          rate = coarseRate[dim];
        } else {
          rate = endRate[dim];
        }
        if (rem > (rate / 2)) {
          ++coarseIJK[dim];
        }
      }

      for (LO dim = 0; dim < 3; ++dim) {
        if (coarseIJK[dim] == coarseNodesPerDim[dim] - 1) {
          nodeIJK[dim] = fineNodesPerDim[dim] - 1;
        } else {
          nodeIJK[dim] = coarseIJK[dim] * coarseRate[dim];
        }
      }
      fineNodeIdx = (nodeIJK[2] * fineNodesPerDim[1] + nodeIJK[1]) * fineNodesPerDim[0] + nodeIJK[0];

      vertex2AggId[interfaceNodes[nodeIdx]] = vertex2AggId[interfaceNodes[fineNodeIdx]];
      procWinner[interfaceNodes[nodeIdx]]   = myRank;
      aggStat[interfaceNodes[nodeIdx]]      = AGGREGATED;
      --numNonAggregatedNodes;
    }  // Loop over interface nodes
  }    // Loop over the interfaces

  // Update aggregates information before subsequent aggregation algorithms
  aggregates->SetNumAggregates(aggregateCount);

  // Set coarse data for next level
  Set(currentLevel, "coarseInterfacesDimensions", coarseInterfacesDimensions);
  Set(currentLevel, "nodeOnCoarseInterface", nodesOnCoarseInterfaces);

}  // BuildInterfaceAggregates()

}  // namespace MueLu

#endif /* MUELU_HYBRIDAGGREGATIONFACTORY_DEF_HPP */
