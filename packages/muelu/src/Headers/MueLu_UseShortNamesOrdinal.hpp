// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node> Aggregates;
#endif
#ifdef MUELU_AGGREGATES_KOKKOS_SHORT
typedef MueLu::Aggregates_kokkos<LocalOrdinal,GlobalOrdinal,Node> Aggregates_kokkos;
#endif
#ifdef MUELU_AGGREGATIONPHASE1ALGORITHM_SHORT
typedef MueLu::AggregationPhase1Algorithm<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase1Algorithm;
#endif
#ifdef MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_SHORT
typedef MueLu::AggregationPhase1Algorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase1Algorithm_kokkos;
#endif
#ifdef MUELU_AGGREGATIONPHASE2AALGORITHM_SHORT
typedef MueLu::AggregationPhase2aAlgorithm<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase2aAlgorithm;
#endif
#ifdef MUELU_AGGREGATIONPHASE2AALGORITHM_KOKKOS_SHORT
typedef MueLu::AggregationPhase2aAlgorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase2aAlgorithm_kokkos;
#endif
#ifdef MUELU_AGGREGATIONPHASE2BALGORITHM_SHORT
typedef MueLu::AggregationPhase2bAlgorithm<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase2bAlgorithm;
#endif
#ifdef MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_SHORT
typedef MueLu::AggregationPhase2bAlgorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase2bAlgorithm_kokkos;
#endif
#ifdef MUELU_AGGREGATIONPHASE3ALGORITHM_SHORT
typedef MueLu::AggregationPhase3Algorithm<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase3Algorithm;
#endif
#ifdef MUELU_AGGREGATIONPHASE3ALGORITHM_KOKKOS_SHORT
typedef MueLu::AggregationPhase3Algorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> AggregationPhase3Algorithm_kokkos;
#endif
#ifdef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_SHORT
typedef MueLu::AggregationStructuredAlgorithm<LocalOrdinal,GlobalOrdinal,Node> AggregationStructuredAlgorithm;
#endif
#ifdef MUELU_AMALGAMATIONINFO_SHORT
typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node> AmalgamationInfo;
#endif
#ifdef MUELU_COUPLEDAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::CoupledAggregationCommHelper<LocalOrdinal,GlobalOrdinal,Node> CoupledAggregationCommHelper;
#endif
#ifdef MUELU_COUPLEDAGGREGATIONFACTORY_SHORT
typedef MueLu::CoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node> CoupledAggregationFactory;
#endif
#ifdef MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_SHORT
typedef MueLu::GlobalLexicographicIndexManager<LocalOrdinal,GlobalOrdinal,Node> GlobalLexicographicIndexManager;
#endif
#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node> Graph;
#endif
#ifdef MUELU_GRAPHBASE_SHORT
typedef MueLu::GraphBase<LocalOrdinal,GlobalOrdinal,Node> GraphBase;
#endif
#ifdef MUELU_HYBRIDAGGREGATIONFACTORY_SHORT
typedef MueLu::HybridAggregationFactory<LocalOrdinal,GlobalOrdinal,Node> HybridAggregationFactory;
#endif
#ifdef MUELU_INDEXMANAGER_SHORT
typedef MueLu::IndexManager<LocalOrdinal,GlobalOrdinal,Node> IndexManager;
#endif
#ifdef MUELU_INTERFACEAGGREGATIONALGORITHM_SHORT
typedef MueLu::InterfaceAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> InterfaceAggregationAlgorithm;
#endif
#ifdef MUELU_ISOLATEDNODEAGGREGATIONALGORITHM_SHORT
typedef MueLu::IsolatedNodeAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> IsolatedNodeAggregationAlgorithm;
#endif
#ifdef MUELU_ISOLATEDNODEAGGREGATIONALGORITHM_KOKKOS_SHORT
typedef MueLu::IsolatedNodeAggregationAlgorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> IsolatedNodeAggregationAlgorithm_kokkos;
#endif
#ifdef MUELU_ISORROPIAINTERFACE_SHORT
typedef MueLu::IsorropiaInterface<LocalOrdinal,GlobalOrdinal,Node> IsorropiaInterface;
#endif
#ifdef MUELU_LWGRAPH_SHORT
typedef MueLu::LWGraph<LocalOrdinal,GlobalOrdinal,Node> LWGraph;
#endif
#ifdef MUELU_LWGRAPH_KOKKOS_SHORT
typedef MueLu::LWGraph_kokkos<LocalOrdinal,GlobalOrdinal,Node> LWGraph_kokkos;
#endif
#ifdef MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
typedef MueLu::LeftoverAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> LeftoverAggregationAlgorithm;
#endif
#ifdef MUELU_LOCALAGGREGATIONALGORITHM_SHORT
typedef MueLu::LocalAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> LocalAggregationAlgorithm;
#endif
#ifdef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_SHORT
typedef MueLu::LocalLexicographicIndexManager<LocalOrdinal,GlobalOrdinal,Node> LocalLexicographicIndexManager;
#endif
#ifdef MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
typedef MueLu::OnePtAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> OnePtAggregationAlgorithm;
#endif
#ifdef MUELU_ONEPTAGGREGATIONALGORITHM_KOKKOS_SHORT
typedef MueLu::OnePtAggregationAlgorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> OnePtAggregationAlgorithm_kokkos;
#endif
#ifdef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_SHORT
typedef MueLu::PreserveDirichletAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node> PreserveDirichletAggregationAlgorithm;
#endif
#ifdef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_KOKKOS_SHORT
typedef MueLu::PreserveDirichletAggregationAlgorithm_kokkos<LocalOrdinal,GlobalOrdinal,Node> PreserveDirichletAggregationAlgorithm_kokkos;
#endif
#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<LocalOrdinal,GlobalOrdinal,Node> PRFactory;
#endif
#ifdef MUELU_REBALANCEMAPFACTORY_SHORT
typedef MueLu::RebalanceMapFactory<LocalOrdinal,GlobalOrdinal,Node> RebalanceMapFactory;
#endif
#ifdef MUELU_REPARTITIONINTERFACE_SHORT
typedef MueLu::RepartitionInterface<LocalOrdinal,GlobalOrdinal,Node> RepartitionInterface;
#endif
#ifdef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
typedef MueLu::UncoupledAggregationFactory<LocalOrdinal,GlobalOrdinal,Node> UncoupledAggregationFactory;
#endif
#ifdef MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_SHORT
typedef MueLu::UncoupledAggregationFactory_kokkos<LocalOrdinal,GlobalOrdinal,Node> UncoupledAggregationFactory_kokkos;
#endif
#ifdef MUELU_UNCOUPLEDINDEXMANAGER_SHORT
typedef MueLu::UncoupledIndexManager<LocalOrdinal,GlobalOrdinal,Node> UncoupledIndexManager;
#endif
#ifdef MUELU_USERAGGREGATIONFACTORY_SHORT
typedef MueLu::UserAggregationFactory<LocalOrdinal,GlobalOrdinal,Node> UserAggregationFactory;
#endif
#ifdef MUELU_FACTORY_SHORT
typedef MueLu::Factory Factory;
#endif

#ifdef MUELU_FACTORYBASE_SHORT
typedef MueLu::FactoryBase FactoryBase;
#endif

#ifdef MUELU_FACTORYMANAGERBASE_SHORT
typedef MueLu::FactoryManagerBase FactoryManagerBase;
#endif

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level Level;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory PFactory;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory RFactory;
#endif

#ifdef MUELU_SINGLELEVELFACTORYBASE_SHORT
typedef MueLu::SingleLevelFactoryBase SingleLevelFactoryBase;
#endif

#ifdef MUELU_TWOLEVELFACTORYBASE_SHORT
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
#endif

#ifdef MUELU_VARIABLECONTAINER_SHORT
typedef MueLu::VariableContainer VariableContainer;
#endif

#ifdef MUELU_SMOOTHERFACTORYBASE_SHORT
typedef MueLu::SmootherFactoryBase SmootherFactoryBase;
#endif

#ifdef MUELU_AMESOSSMOOTHER_SHORT
typedef MueLu::AmesosSmoother<Node> AmesosSmoother;
#endif
#ifdef MUELU_IFPACKSMOOTHER_SHORT
typedef MueLu::IfpackSmoother<Node> IfpackSmoother;
#endif
