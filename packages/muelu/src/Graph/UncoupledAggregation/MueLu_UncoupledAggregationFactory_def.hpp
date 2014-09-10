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
#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UncoupledAggregationFactory_decl.hpp"

#include "MueLu_OnePtAggregationAlgorithm.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm.hpp"
#include "MueLu_MaxLinkAggregationAlgorithm.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm.hpp"
#include "MueLu_EmergencyAggregationAlgorithm.hpp"

#include "MueLu_AggregationPhase1Algorithm.hpp"
#include "MueLu_AggregationPhase2aAlgorithm.hpp"
#include "MueLu_AggregationPhase2bAlgorithm.hpp"
#include "MueLu_AggregationPhase3Algorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::UncoupledAggregationFactory()
  : bDefinitionPhase_(true)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: mode");
    validParamList->getEntry("aggregation: mode").setValidator(
      rcp(new validatorType(Teuchos::tuple<std::string>("new", "old"), "aggregation: mode")));
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
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

    // Aggregation parameters (used in aggregation algorithms)
    // TODO introduce local member function for each aggregation algorithm such that each aggregation algorithm can define its own parameters

    validParamList->set<bool> ("UseOnePtAggregationAlgorithm",             false, "Allow special nodes to be marked for one-to-one transfer to the coarsest level. (default = off)");
    validParamList->set<bool> ("UsePreserveDirichletAggregationAlgorithm", false, "Turn on/off aggregate Dirichlet (isolated nodes) into separate 1pt node aggregates (default = off)");
    validParamList->set<bool> ("UseUncoupledAggregationAlgorithm",          true, "Turn on/off uncoupled aggregation process. Do not turn off: this is "
                               "the main aggregation routine within the uncoupled aggregation process. (default = on)");
    validParamList->set<bool> ("UseMaxLinkAggregationAlgorithm",            true, "Turn on/off MaxLink aggregation algorithm. Adds non-aggregated nodes to "
                               "the next already aggregated neighbour node with the most links. (default = on)");
    validParamList->set<bool> ("UseIsolatedNodeAggregationAlgorithm",       true, "Turn on/off IsolatedNode aggregation algorithm. Ignores isolated "
                               "nodes during aggregation process. (default = on)");
    validParamList->set<bool> ("UseEmergencyAggregationAlgorithm",          true, "Turn on/off Emergency aggregation algorithm. Puts all left over nodes "
                               "into aggregates (including very small aggregates or one-point aggregates). (default = on)");

    validParamList->set< std::string >           ("OnePt aggregate map name",         "", "Name of input map for single node aggregates. (default='')");
    validParamList->set< RCP<const FactoryBase> >("OnePt aggregate map factory",    null, "Generating factory of (DOF) map for single node aggregates.");

    return validParamList;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
    Input(currentLevel, "DofsPerNode");

    const ParameterList& pL = GetParameterList();

    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
    if (mapOnePtName.length() > 0) {
      RCP<const FactoryBase> mapOnePtFact = GetFactory("OnePt aggregate map factory");
      currentLevel.DeclareInput(mapOnePtName, mapOnePtFact.get());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    ParameterList pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    // define aggregation algorithms
    RCP<const FactoryBase> graphFact = GetFactory("Graph");

    // TODO Can we keep different aggregation algorithms over more Build calls?
    algos_.clear();
    algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm(graphFact)));
    if (pL.get<std::string>("aggregation: mode") == "old") {
      if (pL.get<bool>("UseOnePtAggregationAlgorithm")             == true)   algos_.push_back(rcp(new OnePtAggregationAlgorithm             (graphFact)));
      if (pL.get<bool>("UseUncoupledAggregationAlgorithm")         == true)   algos_.push_back(rcp(new AggregationPhase1Algorithm            (graphFact)));
      if (pL.get<bool>("UseMaxLinkAggregationAlgorithm")           == true)   algos_.push_back(rcp(new MaxLinkAggregationAlgorithm           (graphFact)));
      if (pL.get<bool>("UseEmergencyAggregationAlgorithm")         == true)   algos_.push_back(rcp(new EmergencyAggregationAlgorithm         (graphFact)));
                                                                              algos_.push_back(rcp(new IsolatedNodeAggregationAlgorithm      (graphFact)));

    } else {
      if (pL.get<bool>("aggregation: enable phase 1" )             == true)   algos_.push_back(rcp(new AggregationPhase1Algorithm            (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 2a")             == true)   algos_.push_back(rcp(new AggregationPhase2aAlgorithm           (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 2b")             == true)   algos_.push_back(rcp(new AggregationPhase2bAlgorithm           (graphFact)));
      if (pL.get<bool>("aggregation: enable phase 3" )             == true)   algos_.push_back(rcp(new AggregationPhase3Algorithm            (graphFact)));
                                                                              algos_.push_back(rcp(new IsolatedNodeAggregationAlgorithm      (graphFact)));
    }

    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
    RCP<const Map> OnePtMap;
    if (mapOnePtName.length()) {
      RCP<const FactoryBase> mapOnePtFact = GetFactory("OnePt aggregate map factory");
      OnePtMap = currentLevel.Get<RCP<const Map> >(mapOnePtName, mapOnePtFact.get());
    }

    RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");

    // Build
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("UC");

    const LO numRows = graph->GetNodeNumVertices();

    // construct aggStat information
    std::vector<unsigned> aggStat(numRows, READY);

    ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
    if (dirichletBoundaryMap != Teuchos::null)
      for (LO i = 0; i < numRows; i++)
        if (dirichletBoundaryMap[i] == true)
          aggStat[i] = BOUNDARY;

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


    const RCP<const Teuchos::Comm<int> > comm = graph->GetComm();
    GO numGlobalRows = 0;
    if (IsPrint(Statistics1))
      sumAll(comm, as<GO>(numRows), numGlobalRows);

    LO numNonAggregatedNodes = numRows;
    GO numGlobalAggregatedPrev = 0, numGlobalAggsPrev = 0;
    for (size_t a = 0; a < algos_.size(); a++) {
      std::string phase = algos_[a]->description();
      SubFactoryMonitor sfm(*this, "Algo \"" + phase + "\"", currentLevel);

      algos_[a]->BuildAggregates(pL, *graph, *aggregates, aggStat, numNonAggregatedNodes);

      if (IsPrint(Statistics1)) {

        GO numLocalAggregated = numRows - numNonAggregatedNodes, numGlobalAggregated = 0;
        GO numLocalAggs       = aggregates->GetNumAggregates(),  numGlobalAggs = 0;
        sumAll(comm, numLocalAggregated, numGlobalAggregated);
        sumAll(comm, numLocalAggs,       numGlobalAggs);

        double aggPercent = 100*as<double>(numGlobalAggregated)/as<double>(numGlobalRows);
        if (aggPercent > 99.99 && aggPercent < 100.00) {
          // Due to round off (for instance, for 140465733/140466897), we could
          // get 100.00% display even if there are some remaining nodes. This
          // is bad from the users point of view. It is much better to change
          // it to display 99.99%.
          aggPercent = 99.99;
        }
        GetOStream(Statistics1) << "  aggregated : " << (numGlobalAggregated - numGlobalAggregatedPrev) << " (phase), " << std::fixed
                                   << std::setprecision(2) << numGlobalAggregated << "/" << numGlobalRows << " [" << aggPercent << "%] (total)\n"
                                   << "  remaining  : " << numGlobalRows - numGlobalAggregated << "\n"
                                   << "  aggregates : " << numGlobalAggs-numGlobalAggsPrev << " (phase), " << numGlobalAggs << " (total)" << std::endl;
        numGlobalAggregatedPrev = numGlobalAggregated;
        numGlobalAggsPrev       = numGlobalAggs;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError, "MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");

    aggregates->AggregatesCrossProcessors(false);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics0) << aggregates->description() << std::endl;
  }

} //namespace MueLu


#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
