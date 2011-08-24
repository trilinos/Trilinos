// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Graph;
#endif

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Aggregates;
#endif

#ifdef MUELU_LOCALAGGREGATIONFACTORY_SHORT
typedef MueLu::LocalAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LocalAggregationFactory;
#endif

#ifdef MUELU_UCAGGREGATIONFACTORY_SHORT
typedef MueLu::UCAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationFactory;
#endif

#ifdef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::UCAggregationCommHelper<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationCommHelper;
#endif

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level Level;
#endif

#ifdef MUELU_DEFAULTFACTORYHANDLER_SHORT
typedef MueLu::DefaultFactoryHandlerBase DefaultFactoryHandlerBase;
#endif

#ifdef MUELU_TWOLEVELFACTORY_SHORT
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
#endif

#ifdef MUELU_SINGLELEVELFACTORY_SHORT
typedef MueLu::SingleLevelFactoryBase SingleLevelFactoryBase;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory PRFactory;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory PFactory;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory RFactory;
#endif
