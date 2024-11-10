// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_AGGREGATES_SHORT
using Aggregates [[maybe_unused]] = MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONPHASE1ALGORITHM_SHORT
using AggregationPhase1Algorithm [[maybe_unused]] = MueLu::AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONPHASE2AALGORITHM_SHORT
using AggregationPhase2aAlgorithm [[maybe_unused]] = MueLu::AggregationPhase2aAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONPHASE2BALGORITHM_SHORT
using AggregationPhase2bAlgorithm [[maybe_unused]] = MueLu::AggregationPhase2bAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONPHASE3ALGORITHM_SHORT
using AggregationPhase3Algorithm [[maybe_unused]] = MueLu::AggregationPhase3Algorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_SHORT
using AggregationStructuredAlgorithm [[maybe_unused]] = MueLu::AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AMALGAMATIONINFO_SHORT
using AmalgamationInfo [[maybe_unused]] = MueLu::AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_SHORT
using GlobalLexicographicIndexManager [[maybe_unused]] = MueLu::GlobalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_HYBRIDAGGREGATIONFACTORY_SHORT
using HybridAggregationFactory [[maybe_unused]] = MueLu::HybridAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INDEXMANAGER_SHORT
using IndexManager [[maybe_unused]] = MueLu::IndexManager<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INDEXMANAGER_KOKKOS_SHORT
using IndexManager_kokkos [[maybe_unused]] = MueLu::IndexManager_kokkos<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INTERFACEAGGREGATIONALGORITHM_SHORT
using InterfaceAggregationAlgorithm [[maybe_unused]] = MueLu::InterfaceAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INTERFACEMAPPINGTRANSFERFACTORY_SHORT
using InterfaceMappingTransferFactory [[maybe_unused]] = MueLu::InterfaceMappingTransferFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ISORROPIAINTERFACE_SHORT
using IsorropiaInterface [[maybe_unused]] = MueLu::IsorropiaInterface<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LWGRAPH_SHORT
using LWGraph [[maybe_unused]] = MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LWGRAPH_KOKKOS_SHORT
using LWGraph_kokkos [[maybe_unused]] = MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_SHORT
using LocalLexicographicIndexManager [[maybe_unused]] = MueLu::LocalLexicographicIndexManager<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LOCALORDINALTRANSFERFACTORY_SHORT
using LocalOrdinalTransferFactory [[maybe_unused]] = MueLu::LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ONEPTAGGREGATIONALGORITHM_SHORT
using OnePtAggregationAlgorithm [[maybe_unused]] = MueLu::OnePtAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_SHORT
using PreserveDirichletAggregationAlgorithm [[maybe_unused]] = MueLu::PreserveDirichletAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PRFACTORY_SHORT
using PRFactory [[maybe_unused]] = MueLu::PRFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCEMAPFACTORY_SHORT
using RebalanceMapFactory [[maybe_unused]] = MueLu::RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REPARTITIONINTERFACE_SHORT
using RepartitionInterface [[maybe_unused]] = MueLu::RepartitionInterface<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_STRUCTUREDAGGREGATIONFACTORY_KOKKOS_SHORT
using StructuredAggregationFactory_kokkos [[maybe_unused]] = MueLu::StructuredAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UNCOUPLEDAGGREGATIONFACTORY_SHORT
using UncoupledAggregationFactory [[maybe_unused]] = MueLu::UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UNCOUPLEDINDEXMANAGER_SHORT
using UncoupledIndexManager [[maybe_unused]] = MueLu::UncoupledIndexManager<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_USERAGGREGATIONFACTORY_SHORT
using UserAggregationFactory [[maybe_unused]] = MueLu::UserAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACTORY_SHORT
using Factory [[maybe_unused]] = MueLu::Factory;
#endif
#ifdef MUELU_FACTORYBASE_SHORT
using FactoryBase [[maybe_unused]] = MueLu::FactoryBase;
#endif
#ifdef MUELU_FACTORYMANAGERBASE_SHORT
using FactoryManagerBase [[maybe_unused]] = MueLu::FactoryManagerBase;
#endif
#ifdef MUELU_LEVEL_SHORT
using Level [[maybe_unused]] = MueLu::Level;
#endif
#ifdef MUELU_PFACTORY_SHORT
using PFactory [[maybe_unused]] = MueLu::PFactory;
#endif
#ifdef MUELU_RFACTORY_SHORT
using RFactory [[maybe_unused]] = MueLu::RFactory;
#endif
#ifdef MUELU_SINGLELEVELFACTORYBASE_SHORT
using SingleLevelFactoryBase [[maybe_unused]] = MueLu::SingleLevelFactoryBase;
#endif
#ifdef MUELU_TWOLEVELFACTORYBASE_SHORT
using TwoLevelFactoryBase [[maybe_unused]] = MueLu::TwoLevelFactoryBase;
#endif
#ifdef MUELU_VARIABLECONTAINER_SHORT
using VariableContainer [[maybe_unused]] = MueLu::VariableContainer;
#endif
#ifdef MUELU_SMOOTHERFACTORYBASE_SHORT
using SmootherFactoryBase [[maybe_unused]] = MueLu::SmootherFactoryBase;
#endif
#ifdef MUELU_AMESOSSMOOTHER_SHORT
typedef MueLu::AmesosSmoother<Node> AmesosSmoother;
#endif
#ifdef MUELU_IFPACKSMOOTHER_SHORT
typedef MueLu::IfpackSmoother<Node> IfpackSmoother;
#endif
