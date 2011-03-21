// Type definitions for templated classes (generally graph-related) that do not require a scalar.
#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>                          Graph;
#endif

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               Aggregates;
#endif

#ifdef MUELU_AGGREGATIONFACTORY_SHORT
typedef MueLu::AggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       AggregationFactory;
#endif
