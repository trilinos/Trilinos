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
 * MueLu_ExperimentalAggregationFactory_def.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: wiesner
 */

#ifndef MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_

#include <bitset>

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ExperimentalAggregationFactory_decl.hpp"
#include "MueLu_CheapAggregationAlgorithm.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ExperimentalAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ExperimentalAggregationFactory(RCP<const FactoryBase> graphFact)
    : graphFact_(graphFact)
  {
    algo1_ = Teuchos::rcp(new MueLu::CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(graphFact));
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ExperimentalAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("Graph", graphFact_.get(), this); // we should request data...

    if (currentLevel.GetLevelID() == 0) currentLevel.DeclareInput("coarseAggStat", MueLu::NoFactory::get(), this);
    else                                currentLevel.DeclareInput("coarseAggStat", this, this);

  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ExperimentalAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Aggregation (Experimental)", currentLevel);

    RCP<Aggregates> aggregates;
    {
      // Level Get
      RCP<const Graph> graph = currentLevel.Get< RCP<Graph> >("Graph", graphFact_.get());

      // Build
      aggregates = rcp(new Aggregates(*graph));
      aggregates->setObjectLabel("UC");


      const LocalOrdinal nRows = graph->GetNodeNumVertices();

      Teuchos::ArrayRCP<unsigned int> aggStat;
      if(currentLevel.GetLevelID() == 0 && currentLevel.IsAvailable("coarseAggStat",MueLu::NoFactory::get())) {
        aggStat = currentLevel.Get<Teuchos::ArrayRCP<unsigned int> >("coarseAggStat",MueLu::NoFactory::get());
      } else if (currentLevel.IsAvailable("coarseAggStat", this)) {
        aggStat = currentLevel.Get<Teuchos::ArrayRCP<unsigned int> >("coarseAggStat",this);
      } else {
        if(nRows > 0) aggStat = Teuchos::arcp<unsigned int>(nRows);
        for(LocalOrdinal i=0; i<nRows; ++i) {
          aggStat[i] = (unsigned int) 0;
        }
      }

      Teuchos::ArrayRCP<unsigned int> coarse_aggStat = Teuchos::arcp<unsigned int>(nRows);
      for(LocalOrdinal i=0; i<nRows; ++i) {
        coarse_aggStat[i] = (unsigned int) 0;
      }

      // TODO: check return values of functions
      algo1_->PhaseOnePt(*graph,*aggregates,aggStat, coarse_aggStat);
      algo1_->Phase1(*graph,*aggregates,aggStat, coarse_aggStat);
      algo1_->Phase2_maxlink(*graph,*aggregates,aggStat, coarse_aggStat);
      algo1_->Phase3(*graph,*aggregates,aggStat, coarse_aggStat);

      LocalOrdinal numAggs = aggregates->GetNumAggregates();
      coarse_aggStat.resize(Teuchos::as<int>(numAggs));
      currentLevel.Set("coarseAggStat", coarse_aggStat, this);
    }

    // Level Set
    currentLevel.Set("Aggregates", aggregates, this);

    aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());

  }

} //namespace MueLu

#endif /* MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_ */
