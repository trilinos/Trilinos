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

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  HybridAggregationFactory() : bDefinitionPhase_(true)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
#undef  SET_VALID_ENTRY

    return validParamList;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
  Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDAGGREGATION_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }
    out->setShowAllFrontMatter(false).setShowProcRank(true);

    *out << "Entering hybrid aggregation" << std::endl;

    ParameterList pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    // General problem informations are gathered from data stored in the problem matix.
    RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");
    RCP<const Map> fineMap      = graph->GetDomainMap();
    const int myRank            = fineMap->getComm()->getRank();
    const int numRanks          = fineMap->getComm()->getSize();
    const GO  minGlobalIndex    = fineMap->getMinGlobalIndex();

    // Create aggregates object and set basic parameters
    RCP<Aggregates> aggregates = rcp(new Aggregates(fineMap));
    aggregates->setObjectLabel("ST");
    aggregates->AggregatesCrossProcessors(false);
    // std::vector<unsigned> aggStat(geoData->getNumLocalFineNodes(), READY);
    // aggregates->SetNumAggregates(geoData->getNumLocalCoarseNodes());

    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    Set(currentLevel, "Aggregates",         aggregates);

    GetOStream(Statistics1) << aggregates->description() << std::endl;
  }

} //namespace MueLu


#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DEF_HPP_ */
