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
#include "MueLu_IsolatedNodeAggregationAlgorithm.hpp"

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
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_AmalgamationInfo.hpp"


namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  HybridAggregationFactory() : bDefinitionPhase_(true)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    // From UncoupledAggregationFactory
    SET_VALID_ENTRY("aggregation: max agg size");
    SET_VALID_ENTRY("aggregation: min agg size");
    SET_VALID_ENTRY("aggregation: max selected neighbors");
    SET_VALID_ENTRY("aggregation: ordering");
    validParamList->getEntry("aggregation: ordering").setValidator(
      rcp(new validatorType(Teuchos::tuple<std::string>("natural", "graph", "random"), "aggregation: ordering")));
    SET_VALID_ENTRY("aggregation: enable phase 1");
    SET_VALID_ENTRY("aggregation: enable phase 2a");
    SET_VALID_ENTRY("aggregation: enable phase 2b");
    SET_VALID_ENTRY("aggregation: enable phase 3");
    SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
    SET_VALID_ENTRY("aggregation: allow user-specified singletons");
    SET_VALID_ENTRY("aggregation: use interface aggregation");
    SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
    SET_VALID_ENTRY("aggregation: phase3 avoid singletons");

    // From StructuredAggregationFactory
#undef  SET_VALID_ENTRY

    /* From UncoupledAggregation */
    // general variables needed in AggregationFactory
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
    // special variables necessary for OnePtAggregationAlgorithm
    validParamList->set<std::string>            ("OnePt aggregate map name",           "",
                                                 "Name of input map for single node aggregates. (default='')");
    validParamList->set<std::string>            ("OnePt aggregate map factory",        "",
                                                 "Generating factory of (DOF) map for single node aggregates.");

    // InterfaceAggregation parameters
    validParamList->set< std::string >           ("Interface aggregate map name",                "",
                                                  "Name of input map for interface aggregates. (default='')");
    validParamList->set< std::string >           ("Interface aggregate map factory",             "",
                                                  "Generating factory of (DOF) map for interface aggregates.");
    validParamList->set<RCP<const FactoryBase> > ("nodeOnInterface",                  Teuchos::null,
                                                  "Array specifying whether or not a node is on the interface (1 or 0).");

    /* From StructuredAggregation */
    // general variables needed in AggregationFactory
    validParamList->set<std::string >           ("aggregation: mesh layout","Global Lexicographic",
                                                 "Type of mesh ordering");
    validParamList->set<std::string >           ("aggregation: coupling",              "uncoupled",
                                                 "aggregation coupling mode: coupled or uncoupled");
    validParamList->set<int>                    ("aggregation: number of spatial dimensions",    3,
                                                  "The number of spatial dimensions in the problem");
    validParamList->set<int>                    ("aggregation: coarsening order",                0,
                                                  "The interpolation order used to construct grid transfer operators based off these aggregates.");
    validParamList->set<std::string>            ("aggregation: coarsening rate",    "{3}",
                                                  "Coarsening rate per spatial dimensions");
    validParamList->set<RCP<const FactoryBase> >("aggregation: mesh data",  Teuchos::null,
                                                 "Mesh ordering associated data");

    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");


    // Hybrid Aggregation Params
    validParamList->set<RCP<const FactoryBase> >            ("aggregationRegionType",          Teuchos::null,
                                                 "Type of aggregation to use on the region (\"structured\" or \"uncoupled\")");

    return validParamList;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");

    ParameterList pL = GetParameterList();

    /* Hybrid Aggregation */
    if(currentLevel.GetLevelID() == 0){
      if(currentLevel.IsAvailable("aggregationRegionType", NoFactory::get())) {
        currentLevel.DeclareInput("aggregationRegionType", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("aggregationRegionType",NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "Aggregation region type was not provided by the user!");
      }
    } else {
      Input(currentLevel, "aggregationRegionType");
    }


    /* UncoupledAggregation */
    Input(currentLevel, "DofsPerNode");

    // request special data necessary for InterfaceAggregation
    if (pL.get<bool>("aggregation: use interface aggregation")       == true){
      if(currentLevel.GetLevelID() == 0) {
        if(currentLevel.IsAvailable("nodeOnInterface", NoFactory::get())) {
          currentLevel.DeclareInput("nodeOnInterface", NoFactory::get(), this);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("nodeOnInterface", NoFactory::get()),
                                       Exceptions::RuntimeError,
                                       "nodeOnInterface was not provided by the user on level0!");
        }
      } else {
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



    /* StructuredAggregation */
    std::string coupling = pL.get<std::string>("aggregation: coupling");
    const bool coupled = (coupling == "coupled" ? true : false);
    if(coupled) {
      // Request the global number of nodes per dimensions
      if(currentLevel.GetLevelID() == 0) {
        if(currentLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
          currentLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                     Exceptions::RuntimeError,
                                     "gNodesPerDim was not provided by the user on level0!");
        }
      } else {
        Input(currentLevel, "gNodesPerDim");
      }
    }

    // Request the local number of nodes per dimensions
    if(currentLevel.GetLevelID() == 0) {
      if(currentLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
        currentLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(!currentLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "lNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "lNodesPerDim");
    }


  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_HYBRIDAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }
    out->setShowAllFrontMatter(false).setShowProcRank(true);

    *out << "Entering hybrid aggregation" << std::endl;

    ParameterList pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    if (pL.get<int>("aggregation: max agg size") == -1)
      pL.set("aggregation: max agg size", INT_MAX);

    // define aggregation algorithms
    RCP<const FactoryBase> graphFact = GetFactory("Graph");

    // General problem informations are gathered from data stored in the problem matix.
    RCP<const GraphBase> graph  = Get< RCP<GraphBase> >(currentLevel, "Graph");
    RCP<const Map> fineMap      = graph->GetDomainMap();
    const int myRank            = fineMap->getComm()->getRank();
    const int numRanks          = fineMap->getComm()->getSize();
    const GO  minGlobalIndex    = fineMap->getMinGlobalIndex();

    // Build aggregates
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("HB");

    // construct aggStat information
    const LO numRows = graph->GetNodeNumVertices();
    std::vector<unsigned> aggStat(numRows, READY);

    // Get aggregation type for region
    std::string regionType;
    if(currentLevel.GetLevelID() == 0) {
      // On level 0, data is provided by applications and has no associated factory.
      regionType = currentLevel.Get< std::string >("aggregationRegionType", NoFactory::get());
    } else {
      // On level > 0, data is provided directly by generating factories.
      regionType = Get< std::string >(currentLevel, "aggregationRegionType");
    }
    *out<<"p="<< myRank << " | "<<regionType<<" | regionType determined" << std::endl;

    algos_.clear();
    if (regionType == "structured") {
      // Add AggregationStructuredAlgorithm
      algos_.push_back(rcp(new AggregationStructuredAlgorithm(graphFact)));

      // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
      // obtain a nodeMap.
      const int numDimensions      = pL.get<int>("aggregation: number of spatial dimensions");
      const int interpolationOrder = pL.get<int>("aggregation: coarsening order");
      std::string meshLayout       = pL.get<std::string>("aggregation: mesh layout");
      std::string coupling         = pL.get<std::string>("aggregation: coupling");
      const bool coupled = false; //Only support uncoupled
      Array<GO> gFineNodesPerDir(3);
      Array<LO> lFineNodesPerDir(3);
      if(currentLevel.GetLevelID() == 0) {
        // On level 0, data is provided by applications and has no associated factory.
        if(coupled) {
          gFineNodesPerDir = currentLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
        }
        lFineNodesPerDir = currentLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
      } else {
        // On level > 0, data is provided directly by generating factories.
        if(coupled) {
          gFineNodesPerDir = Get<Array<GO> >(currentLevel, "gNodesPerDim");
        }
        lFineNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
      }

      // Set lFineNodesPerDir to 1 for directions beyond numDimensions
      for(int dim = numDimensions; dim < 3; ++dim) {
        lFineNodesPerDir[dim] = 1;
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

      // Now that we have extracted info from the level, create the IndexManager
      RCP<MueLu::IndexManager<LO,GO,NO> > geoData;
      if(!coupled) {
        geoData = rcp(new MueLu::UncoupledIndexManager<LO,GO,NO>(fineMap->getComm(),
                                                                 coupled,
                                                                 numDimensions,
                                                                 interpolationOrder,
                                                                 myRank,
                                                                 numRanks,
                                                                 gFineNodesPerDir,
                                                                 lFineNodesPerDir,
                                                                 coarseRate));
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(coupled, Exceptions::RuntimeError,
                "Coupled aggregation is not yet implemented in hybrid aggregation");
      }

      *out << "The index manager has now been built" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getNodeNumElements()
                                 != static_cast<size_t>(geoData->getNumLocalFineNodes()),
                                 Exceptions::RuntimeError,
                                 "The local number of elements in the graph's map is not equal to "
                                 "the number of nodes given by: lNodesPerDim!");
      if(coupled) {
        TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getGlobalNumElements()
                                   != static_cast<size_t>(geoData->getNumGlobalFineNodes()),
                                   Exceptions::RuntimeError,
                                   "The global number of elements in the graph's map is not equal to "
                                   "the number of nodes given by: gNodesPerDim!");
      }

      *out << "Compute coarse mesh data" << std::endl;
      std::vector<std::vector<GO> > coarseMeshData = geoData->getCoarseMeshData();

      aggregates->SetIndexManager(geoData);
      aggregates->SetNumAggregates(geoData->getNumLocalCoarseNodes());

      Set(currentLevel, "gCoarseNodesPerDim", geoData->getGlobalCoarseNodesPerDir());
      Set(currentLevel, "lCoarseNodesPerDim", geoData->getLocalCoarseNodesPerDir());

    }// end structured aggregation setup

    if (regionType == "uncoupled"){
      // Add unstructred aggregation phases
      algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm(graphFact)));
      if (pL.get<bool>("aggregation: use interface aggregation")       == true)   algos_.push_back(rcp(new InterfaceAggregationAlgorithm         (graphFact)));
      if (pL.get<bool>("aggregation: allow user-specified singletons") == true)   algos_.push_back(rcp(new OnePtAggregationAlgorithm             (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 1" )                 == true)   algos_.push_back(rcp(new AggregationPhase1Algorithm            (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 2a")                 == true)   algos_.push_back(rcp(new AggregationPhase2aAlgorithm           (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 2b")                 == true)   algos_.push_back(rcp(new AggregationPhase2bAlgorithm           (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 3" )                 == true)   algos_.push_back(rcp(new AggregationPhase3Algorithm            (graphFact)));


      // Set map for interface aggregates
      std::string mapInterfaceName = pL.get<std::string>("Interface aggregate map name");
      RCP<Map> InterfaceMap = Teuchos::null;
      // interface
      if (pL.get<bool>("aggregation: use interface aggregation")       == true){
        Teuchos::Array<LO> nodeOnInterface = Get<Array<LO>>(currentLevel,"nodeOnInterface");
        for (LO i = 0; i < numRows; i++) {
          if (nodeOnInterface[i])
            aggStat[i] = INTERFACE;
        }
      }

      // Dirichlet boundary
      ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
      if (dirichletBoundaryMap != Teuchos::null)
        for (LO i = 0; i < numRows; i++)
          if (dirichletBoundaryMap[i] == true)
            aggStat[i] = BOUNDARY;

      // OnePt aggregation
      std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
      RCP<Map> OnePtMap = Teuchos::null;
      if (mapOnePtName.length()) {
        std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
        if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
          OnePtMap = currentLevel.Get<RCP<Map> >(mapOnePtName, NoFactory::get());
        } else {
          RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
          OnePtMap = currentLevel.Get<RCP<Map> >(mapOnePtName, mapOnePtFact.get());
        }
      }
      LO nDofsPerNode = Get<LO>(currentLevel, "DofsPerNode");
      GO indexBase = graph->GetDomainMap()->getIndexBase();
      if (OnePtMap != Teuchos::null) {
        for (LO i = 0; i < numRows; i++) {
          // reconstruct global row id (FIXME only works for contiguous maps)
          GO grid = (graph->GetDomainMap()->getGlobalElement(i)-indexBase) * nDofsPerNode + indexBase;
          for (LO kr = 0; kr < nDofsPerNode; kr++)
            if (OnePtMap->isNodeGlobalElement(grid + kr))
              aggStat[i] = ONEPT;
        }
      }

      // Create a fake lCoarseNodesPerDir for CoordinatesTranferFactory
      Array<LO> lCoarseNodesPerDir(3,-1);
      Set(currentLevel, "lCoarseNodesPerDim", lCoarseNodesPerDir);
    }// end uncoupled aggregation setup

    aggregates->AggregatesCrossProcessors(false); // No coupled aggregation

    LO numNonAggregatedNodes = numRows;
    for (size_t a = 0; a < algos_.size(); a++) {
      std::string phase = algos_[a]->description();
      SubFactoryMonitor sfm(*this, "Algo \"" + phase + "\"", currentLevel);
      *out << "p=" << myRank << " | "<<regionType<<" | Executing phase " << a << std::endl;

      int oldRank = algos_[a]->SetProcRankVerbose(this->GetProcRankVerbose());
      algos_[a]->BuildAggregates(pL, *graph, *aggregates, aggStat, numNonAggregatedNodes);
      algos_[a]->SetProcRankVerbose(oldRank);
      *out << "p=" << myRank << " | "<<regionType<<" | Done Executing phase " << a << std::endl;
    }

    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    Set(currentLevel, "Aggregates",         aggregates);

    Set(currentLevel, "aggregationRegionTypeCoarse", regionType);

    GetOStream(Statistics1) << aggregates->description() << std::endl;
  }

} //namespace MueLu


#endif /* MUELU_HYBRIDAGGREGATIONFACTORY_DEF_HPP */
