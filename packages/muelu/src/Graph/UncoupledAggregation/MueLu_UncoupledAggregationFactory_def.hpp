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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_UncoupledAggregationFactory_def.hpp
 *
 *  Created on: Sep 17, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UncoupledAggregationFactory_decl.hpp"

#include "MueLu_OnePtAggregationAlgorithm.hpp"
#include "MueLu_SmallAggregationAlgorithm.hpp"
#include "MueLu_UncoupledAggregationAlgorithm.hpp"
#include "MueLu_MaxLinkAggregationAlgorithm.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm.hpp"
#include "MueLu_EmergencyAggregationAlgorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UncoupledAggregationFactory()
: bDefinitionPhase_(true)
  {
  }

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const ParameterList> UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // aggregate parameters (used in aggregation algorithms)
  // TODO introduce local member function for each aggregation algorithm
  //      such that each aggregation algorithm can define its own parameters
  validParamList->set<AggOptions::Ordering>("Ordering", AggOptions::NATURAL, "Ordering strategy (NATURAL|GRAPH|RANDOM)");
  validParamList->set<LocalOrdinal> ("MaxNeighAlreadySelected", 0, "Number of maximum neighbour nodes that are already aggregated already. If a new aggregate has some neighbours that are already aggregated, this node probably can be added to one of these aggregates. We don't need a new one.");
  validParamList->set<LocalOrdinal> ("MinNodesPerAggregate", 2, "Minimum number of nodes for aggregate");

  validParamList->set<bool> ("UseOnePtAggregationAlgorithm", true, "Allow special nodes to be marked for one-to-one transfer to the coarsest level. (default = on)");
  validParamList->set<bool> ("UseSmallAggregatesAggregationAlgorithm", false, "Turn on/off build process for small aggregates in user defined regions. (default = off)");
  validParamList->set<bool> ("UseUncoupledAggregationAlgorithm", true, "Turn on/off uncoupled aggregation process. Do not turn off: this is the main aggregation routine within the uncoupled aggregation process. (default = on)");
  validParamList->set<bool> ("UseMaxLinkAggregationAlgorithm", true, "Turn on/off MaxLink aggregation algorithm. Adds non-aggregated nodes to the next already aggregated neighbour node with the most links. (default = on)");
  validParamList->set<bool> ("UseIsolatedNodeAggregationAlgorithm", true, "Turn on/off IsolatedNode aggregation algorithm. Ignores isolated nodes during aggregation process. (default = on)");
  validParamList->set<bool> ("UseEmergencyAggregationAlgorithm", true, "Turn on/off Emergency aggregation algorithm. Puts all left over nodes into aggregates (including very small aggregates or one-point aggregates). (default = on)");

  // input parameters
  validParamList->set< RCP<const FactoryBase> >("Graph", Teuchos::null, "Generating factory of the graph");
  validParamList->set< RCP<const FactoryBase> >("DofsPerNode", Teuchos::null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

  validParamList->set< std::string >           ("OnePt aggregate map name", "", "Name of input map for single node aggregates. (default='')");
  validParamList->set< RCP<const FactoryBase> >("OnePt aggregate map factory", Teuchos::null, "Generating factory of (DOF) map for single node aggregates.");
  validParamList->set< std::string >           ("SmallAgg aggregate map name", "", "Name of input map for small aggregates. (default='')");
  validParamList->set< RCP<const FactoryBase> >("SmallAgg aggregate map factory", Teuchos::null, "Generating factory of (DOF) map for small aggregates.");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");

  const ParameterList & pL = GetParameterList();
  std::string mapOnePtName                     = pL.get<std::string> ("OnePt aggregate map name");
  Teuchos::RCP<const FactoryBase> mapOnePtFact = GetFactory          ("OnePt aggregate map factory");
  std::string mapSmallAggName                     = pL.get<std::string> ("SmallAgg aggregate map name");
  Teuchos::RCP<const FactoryBase> mapSmallAggFact = GetFactory          ("SmallAgg aggregate map factory");

  if(mapOnePtName.length() > 0)
    currentLevel.DeclareInput(mapOnePtName,mapOnePtFact.get());
  if(mapSmallAggName.length() > 0)
    currentLevel.DeclareInput(mapSmallAggName,mapSmallAggFact.get());

}

/*template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Append(const RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & alg) {
  TEUCHOS_TEST_FOR_EXCEPTION(bDefinitionPhase_==false,Exceptions::RuntimeError,"MueLu::UncoupledAggregationFactory::Build: cannot call Append after Build. Error.");
  algos_.push_back(alg);
}*/

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
{
  FactoryMonitor m(*this, "Build", currentLevel);

  const ParameterList & pL = GetParameterList();
  std::string mapOnePtName    = pL.get<std::string> ("OnePt aggregate map name");
  std::string mapSmallAggName = pL.get<std::string> ("SmallAgg aggregate map name");
  Teuchos::RCP<const FactoryBase> mapOnePtFact    = GetFactory  ("OnePt aggregate map factory");
  Teuchos::RCP<const FactoryBase> mapSmallAggFact = GetFactory ("SmallAgg aggregate map factory");

  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  bool bUseOnePtAggregationAlgorithm        = pL.get<bool> ("UseOnePtAggregationAlgorithm");
  bool bUseSmallAggregationAlgorithm        = pL.get<bool> ("UseSmallAggregatesAggregationAlgorithm");
  bool bUseUncoupledAggregationAglorithm    = pL.get<bool> ("UseUncoupledAggregationAlgorithm");
  bool bUseMaxLinkAggregationAlgorithm      = pL.get<bool> ("UseMaxLinkAggregationAlgorithm");
  bool bUseIsolatedNodeAggregationAglorithm = pL.get<bool> ("UseIsolatedNodeAggregationAlgorithm");
  bool bUseEmergencyAggregationAlgorithm    = pL.get<bool> ("UseEmergencyAggregationAlgorithm");

  // define aggregation algorithms
  Teuchos::RCP<const FactoryBase> graphFact = GetFactory("Graph");
  algos_.clear();  // TODO can we keep different aggregation algorithms over more Build calls?
  if (bUseOnePtAggregationAlgorithm)        algos_.push_back(Teuchos::rcp(new MueLu::OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bUseSmallAggregationAlgorithm)        algos_.push_back(Teuchos::rcp(new MueLu::SmallAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bUseUncoupledAggregationAglorithm)    algos_.push_back(Teuchos::rcp(new MueLu::UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bUseMaxLinkAggregationAlgorithm)      algos_.push_back(Teuchos::rcp(new MueLu::MaxLinkAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bUseIsolatedNodeAggregationAglorithm) algos_.push_back(Teuchos::rcp(new MueLu::IsolatedNodeAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bUseEmergencyAggregationAlgorithm)    algos_.push_back(Teuchos::rcp(new MueLu::EmergencyAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));


  RCP<Aggregates> aggregates;
  {
    // Level Get
    RCP<const GraphBase> graph   = Get< RCP<GraphBase> >(currentLevel, "Graph");
    LocalOrdinal nDofsPerNode = Get<LocalOrdinal>(currentLevel, "DofsPerNode");

    // TODO create a map of Xpetra::Maps for different
    // aggregation information (OnePtAggregegates...)
    Teuchos::RCP<const Map> OnePtMap = Teuchos::null;
    if(mapOnePtName.length() > 0) {
      OnePtMap = currentLevel.Get<Teuchos::RCP<const Map> >(mapOnePtName,mapOnePtFact.get());
    }
    Teuchos::RCP<const Map> SmallAggMap = Teuchos::null;
    if(mapSmallAggName.length() > 0) {
      SmallAggMap = currentLevel.Get<Teuchos::RCP<const Map> >(mapSmallAggName,mapSmallAggFact.get());
    }

    // Build
    aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("UC");


    const LocalOrdinal nRows = graph->GetNodeNumVertices();

    // construct aggStat information
    Teuchos::ArrayRCP<unsigned int> aggStat;
    if(nRows > 0) aggStat = Teuchos::arcp<unsigned int>(nRows);
    ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
    if (dirichletBoundaryMap == Teuchos::null)
      dirichletBoundaryMap = ArrayRCP<bool>(nRows,false);
    for(LocalOrdinal i=0; i<nRows; ++i) {
      if (dirichletBoundaryMap[i] == false)
        aggStat[i] = NodeStats::READY;
      else
        aggStat[i] = NodeStats::BOUNDARY;
      GlobalOrdinal grid = graph->GetDomainMap()->getGlobalElement(i) * nDofsPerNode;
      if(SmallAggMap != Teuchos::null) {
         // reconstruct global row id (FIXME only works for contiguous maps)
        for(LocalOrdinal kr = 0; kr < nDofsPerNode; kr++) {
          if(SmallAggMap->isNodeGlobalElement(grid+kr)) {
            aggStat[i] = MueLu::NodeStats::SMALLAGG;
          }
        }
      }
      if(OnePtMap != Teuchos::null) {
        // reconstruct global row id (FIXME only works for contiguous maps)
        for(LocalOrdinal kr = 0; kr < nDofsPerNode; kr++) {
          if(OnePtMap->isNodeGlobalElement(grid+kr)) {
            aggStat[i] = MueLu::NodeStats::ONEPT;
          }
        }
      }
    }

    // TODO: check return values of functions
    LocalOrdinal nonAggregatedNodes = -1;

    //Teuchos::ParameterList params;
    for(size_t a = 0; a < algos_.size(); a++) {
      nonAggregatedNodes = algos_[a]->BuildAggregates(pL,*graph,*aggregates,aggStat);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(nonAggregatedNodes > 0,Exceptions::RuntimeError,"MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");
  }

  aggregates->AggregatesCrossProcessors(false);
  // Level Set
  Set(currentLevel, "Aggregates", aggregates);

  GetOStream(Statistics0, 0) << aggregates->description() << std::endl;
}

} //namespace MueLu


#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
