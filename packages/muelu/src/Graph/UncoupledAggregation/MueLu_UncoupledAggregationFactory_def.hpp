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
#include "MueLu_EmergencyAggregationAlgorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UncoupledAggregationFactory(RCP<const FactoryBase> graphFact, bool bMaxLinkAggregation, bool bEmergencyAggregation)
: bDefinitionPhase_(true), mapOnePtName_(""), mapOnePtFact_(Teuchos::null), mapSmallAggName_(""), mapSmallAggFact_(Teuchos::null)
  {
  SetFactory("Graph", graphFact);       // for compatibility with old code
  SetFactory("DofsPerNode", graphFact); // for compatibility with old code
  algos_.push_back(Teuchos::rcp(new MueLu::OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  algos_.push_back(Teuchos::rcp(new MueLu::SmallAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  algos_.push_back(Teuchos::rcp(new MueLu::UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bMaxLinkAggregation)   algos_.push_back(Teuchos::rcp(new MueLu::MaxLinkAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
  if (bEmergencyAggregation) algos_.push_back(Teuchos::rcp(new MueLu::EmergencyAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact)));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");

  if(mapOnePtName_.length() > 0)
    currentLevel.DeclareInput(mapOnePtName_,mapOnePtFact_.get());
  if(mapSmallAggName_.length() > 0)
    currentLevel.DeclareInput(mapSmallAggName_,mapSmallAggFact_.get());

}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Append(const RCP<MueLu::AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & alg) {
  TEUCHOS_TEST_FOR_EXCEPTION(bDefinitionPhase_==false,Exceptions::RuntimeError,"MueLu::UncoupledAggregationFactory::Build: cannot call Append after Build. Error.");
  algos_.push_back(alg);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
{
  FactoryMonitor m(*this, "Build", currentLevel);

  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  RCP<Aggregates> aggregates;
  {
    // Level Get
    RCP<const Graph> graph   = Get< RCP<Graph> >(currentLevel, "Graph");
    LocalOrdinal nDofsPerNode = Get<LocalOrdinal>(currentLevel, "DofsPerNode");

    // TODO create a map of Xpetra::Maps for different
    // aggregation information (OnePtAggregegates...)
    Teuchos::RCP<const Map> OnePtMap = Teuchos::null;
    if(mapOnePtName_.length() > 0) {
      OnePtMap = currentLevel.Get<Teuchos::RCP<const Map> >(mapOnePtName_,mapOnePtFact_.get());
    }
    Teuchos::RCP<const Map> SmallAggMap = Teuchos::null;
    if(mapSmallAggName_.length() > 0) {
      SmallAggMap = currentLevel.Get<Teuchos::RCP<const Map> >(mapSmallAggName_,mapSmallAggFact_.get());
    }

    // Build
    aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("UC");


    const LocalOrdinal nRows = graph->GetNodeNumVertices();

    // construct aggStat information
    Teuchos::ArrayRCP<unsigned int> aggStat;
    if(nRows > 0) aggStat = Teuchos::arcp<unsigned int>(nRows);
    for(LocalOrdinal i=0; i<nRows; ++i) {
      aggStat[i] = NodeStats::READY;
      GlobalOrdinal grid = graph->GetDomainMap()->getGlobalElement(i) * nDofsPerNode;
      if(SmallAggMap != Teuchos::null) {
         // reconstruct global row id (FIXME only works for contingoous maps)
         if(SmallAggMap->isNodeGlobalElement(grid)) {
           aggStat[i] = MueLu::NodeStats::SMALLAGG;
         }
       }
      if(OnePtMap != Teuchos::null) {
        // reconstruct global row id (FIXME only works for contingoous maps)
        if(OnePtMap->isNodeGlobalElement(grid)) {
          aggStat[i] = MueLu::NodeStats::ONEPT;
        }
      }
    }

    // TODO: check return values of functions
    LocalOrdinal nonAggregatedNodes = -1;
    for(size_t a = 0; a < algos_.size(); a++) {
      nonAggregatedNodes = algos_[a]->BuildAggregates(*graph,*aggregates,aggStat);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(nonAggregatedNodes > 0,Exceptions::RuntimeError,"MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");
  }

  // Level Set
  Set(currentLevel, "Aggregates", aggregates);

  aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());

}

} //namespace MueLu


#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
