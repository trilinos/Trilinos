/*
 * MueLu_ExperimentalAggregationFactory_def.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: wiesner
 */

#ifndef MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_

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
      Teuchos::ArrayRCP<NodeState> aggStat;

      if(currentLevel.GetLevelID() == 0 && currentLevel.IsAvailable("coarseAggStat",MueLu::NoFactory::get())) {
        aggStat = currentLevel.Get<Teuchos::ArrayRCP<NodeState> >("coarseAggStat",MueLu::NoFactory::get());
      } else if (currentLevel.IsAvailable("coarseAggStat", this)) {
        std::cout << "coarseAggStat: found on level" << std::endl;
        aggStat = currentLevel.Get<Teuchos::ArrayRCP<NodeState> >("coarseAggStat",this);
      } else {
        std::cout << "use default coarseAggStat" << std::endl;
        if(nRows > 0) aggStat = Teuchos::arcp<NodeState>(nRows);
        for(LocalOrdinal i=0; i<nRows; ++i) {
          aggStat[i] = READY;
        }
      }

      // we cannot have more aggregates than nodes on the current proc
      Teuchos::ArrayRCP<NodeState> coarse_aggStat = Teuchos::arcp<NodeState>(nRows);

      LocalOrdinal ret = -1;
      ret = algo1_->Phase1a(*graph,*aggregates,aggStat, coarse_aggStat);
      ret = algo1_->Phase2_maxlink(*graph,*aggregates,aggStat);
      ret = algo1_->Phase3(*graph,*aggregates,aggStat);


      LocalOrdinal numAggs = aggregates->GetNumAggregates();
      coarse_aggStat.resize(Teuchos::as<int>(numAggs));

      currentLevel.Set("coarseAggStat", coarse_aggStat, this);

      //ret = algo1_->Phase4(*graph,*aggregates,aggStat);

    }

    // Level Set
    currentLevel.Set("Aggregates", aggregates, this);

    aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());

  }

} //namespace MueLu

#endif /* MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_ */
