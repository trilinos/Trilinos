/*
 * MueLu_ExperimentalAggregationFactory_def.hpp
 *
 *  Created on: Jul 26, 2012
 *      Author: wiesner
 */

#ifndef MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_


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
    algo1_ = Teuchos::rcp(new CheapAggregationAlgorithm(graphFact));
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ExperimentalAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("Graph", graphFact_.get(), this); // we should request data...
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

      algo1_->CoarsenUncoupled(*graph, *aggregates);

    }

    // Level Set
    currentLevel.Set("Aggregates", aggregates, this);

    /*if (IsPrint(Statistics0)) {
      aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());
    }*/
    aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());

  }

} //namespace MueLu

#endif /* MUELU_EXPERIMENTALAGGREGATIONFACTORY_DEF_HPP_ */
