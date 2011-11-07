#ifndef MUELU_UCAGGREGATIONFACTORY_DEF_HPP
#define MUELU_UCAGGREGATIONFACTORY_DEF_HPP

#include "MueLu_UCAggregationFactory_decl.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UCAggregationFactory(RCP<FactoryBase> graphFact)
    : graphFact_(graphFact)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(algo2_.GetMinNodesPerAggregate() != algo1_.GetMinNodesPerAggregate(), Exceptions::RuntimeError, "");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    //if(currentLevel.IsAvailable("Aggregates",this)) return; //TODO: Why??????

    currentLevel.DeclareInput("Graph", graphFact_.get()); // we should request data...
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void UCAggregationFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const
  {
    Monitor m(*this, "Aggregation");

    //TODO check for reuse of aggregates here
    //FIXME should there be some way to specify the name of the graph in the needs table, i.e., could
    //FIXME there ever be more than one graph?
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("UCAggregationFactory::Build_" + Teuchos::toString(currentLevel.GetLevelID())));
    timer->start(true);

    // Level Get
    RCP<const Graph> graph = currentLevel.Get< RCP<Graph> >("Graph", graphFact_.get());

    // Build
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph)); 
    aggregates->setObjectLabel("UC");

    algo1_.CoarsenUncoupled(*graph, *aggregates);
    algo2_.AggregateLeftovers(*graph, *aggregates);

    // Level Set
    currentLevel.Set("Aggregates", aggregates, this);

    //
    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(graph->GetComm()));

    if (IsPrint(Statistics0)) {
      aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());
    }

  }

} //namespace MueLu

#endif // MUELU_UCAGGREGATIONFACTORY_DEF_HPP
