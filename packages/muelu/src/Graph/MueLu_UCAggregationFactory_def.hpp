#ifndef MUELU_UCAGGREGATIONFACTORY_DEF_HPP
#define MUELU_UCAGGREGATIONFACTORY_DEF_HPP

#include "MueLu_UCAggregationFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UCAggregationFactory(RCP<const FactoryBase> graphFact)
    : graphFact_(graphFact)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(algo2_.GetMinNodesPerAggregate() != algo1_.GetMinNodesPerAggregate(), Exceptions::RuntimeError, "");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    //if(currentLevel.IsAvailable("Aggregates",this)) return; //TODO: Why??????

    currentLevel.DeclareInput("Graph", graphFact_.get(), this); // we should request data...
    //currentLevel.DeclareInput("UnAmalgamationInfo", graphFact_.get(), this); // TODO, only provided by CoalesceDropFactory2
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Aggregation", currentLevel);

    RCP<Aggregates> aggregates;
    {
      //TODO check for reuse of aggregates here
      //FIXME should there be some way to specify the name of the graph in the needs table, i.e., could
      //FIXME there ever be more than one graph?
      //FIXME TAW: The graph is always labeled with "Graph". There can be more than one graph of course
      //FIXME TAW: We can distinguish them by their factory!
      
      // Level Get
      RCP<const Graph> graph = currentLevel.Get< RCP<Graph> >("Graph", graphFact_.get());
      
      //if(currentLevel.IsAvailable("UnAmalgamationInfo", graphFact_.get())) {
      //  RCP<const AmalgamationInfo> graph = currentLevel.Get< RCP<AmalgamationInfo> >("UnAmalgamationInfo", graphFact_.get());
      //}

      // Build
      aggregates = rcp(new Aggregates(*graph)); 
      aggregates->setObjectLabel("UC");
      
      algo1_.CoarsenUncoupled(*graph, *aggregates);
      algo2_.AggregateLeftovers(*graph, *aggregates);

    }

    // Level Set
    currentLevel.Set("Aggregates", aggregates, this);

    if (IsPrint(Statistics0)) {
      aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());
    }

  }

} //namespace MueLu

#endif // MUELU_UCAGGREGATIONFACTORY_DEF_HPP
