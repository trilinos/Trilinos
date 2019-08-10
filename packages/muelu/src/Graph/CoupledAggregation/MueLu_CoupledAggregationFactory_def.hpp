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
#ifndef MUELU_COUPLEDAGGREGATIONFACTORY_DEF_HPP
#define MUELU_COUPLEDAGGREGATIONFACTORY_DEF_HPP

#include "MueLu_CoupledAggregationFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::CoupledAggregationFactory()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(algo2_.GetMinNodesPerAggregate() != algo1_.GetMinNodesPerAggregate(), Exceptions::RuntimeError, "");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "Graph");
    //Input(currentLevel, "UnAmalgamationInfo"); // TODO, only provided by CoalesceDropFactory2
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Aggregates> aggregates;
    {
      //TODO check for reuse of aggregates here

      // Level Get
      RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");

      // Build
      aggregates = rcp(new Aggregates(*graph));
      aggregates->setObjectLabel("UC");

      algo1_.CoarsenUncoupled(*graph, *aggregates);
      algo2_.AggregateLeftovers(*graph, *aggregates);

    }

    aggregates->AggregatesCrossProcessors(true);
    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    // Level Set
    Set(currentLevel, "Aggregates", aggregates);

    if (IsPrint(Statistics0)) {
      aggregates->describe(GetOStream(Statistics0), getVerbLevel());
    }

  }

} //namespace MueLu

#endif // MUELU_COUPLEDAGGREGATIONFACTORY_DEF_HPP
