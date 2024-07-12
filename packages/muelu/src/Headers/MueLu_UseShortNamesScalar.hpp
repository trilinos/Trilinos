// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node of the current context.

#include <Xpetra_UseShortNamesScalar.hpp>

#ifdef MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_SHORT
using AdaptiveSaMLParameterListInterpreter [[maybe_unused]] = MueLu::AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
using AggregationExportFactory [[maybe_unused]] = MueLu::AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AGGREGATEQUALITYESTIMATEFACTORY_SHORT
using AggregateQualityEstimateFactory [[maybe_unused]] = MueLu::AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AMALGAMATIONFACTORY_SHORT
using AmalgamationFactory [[maybe_unused]] = MueLu::AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AMESOS2SMOOTHER_SHORT
using Amesos2Smoother [[maybe_unused]] = MueLu::Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_AMGXOPERATOR_SHORT
using AMGXOperator [[maybe_unused]] = MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ALGEBRAICPERMUTATIONSTRATEGY_SHORT
using AlgebraicPermutationStrategy [[maybe_unused]] = MueLu::AlgebraicPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BELOSSMOOTHER_SHORT
using BelosSmoother [[maybe_unused]] = MueLu::BelosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLACKBOXPFACTORY_SHORT
using BlackBoxPFactory [[maybe_unused]] = MueLu::BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDCOARSEMAPFACTORY_SHORT
using BlockedCoarseMapFactory [[maybe_unused]] = MueLu::BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDCOORDINATESTRANSFERFACTORY_SHORT
using BlockedCoordinatesTransferFactory [[maybe_unused]] = MueLu::BlockedCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDDIRECTSOLVER_SHORT
using BlockedDirectSolver [[maybe_unused]] = MueLu::BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_SHORT
using BlockedGaussSeidelSmoother [[maybe_unused]] = MueLu::BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDJACOBISMOOTHER_SHORT
using BlockedJacobiSmoother [[maybe_unused]] = MueLu::BlockedJacobiSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDPFACTORY_SHORT
using BlockedPFactory [[maybe_unused]] = MueLu::BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BLOCKEDRAPFACTORY_SHORT
using BlockedRAPFactory [[maybe_unused]] = MueLu::BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BRICKAGGREGATIONFACTORY_SHORT
using BrickAggregationFactory [[maybe_unused]] = MueLu::BrickAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_BRAESSSARAZINSMOOTHER_SHORT
using BraessSarazinSmoother [[maybe_unused]] = MueLu::BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CGSOLVER_SHORT
using CGSolver [[maybe_unused]] = MueLu::CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CLASSICALMAPFACTORY_SHORT
using ClassicalMapFactory [[maybe_unused]] = MueLu::ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CLASSICALPFACTORY_SHORT
using ClassicalPFactory [[maybe_unused]] = MueLu::ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CLONEREPARTITIONINTERFACE_SHORT
using CloneRepartitionInterface [[maybe_unused]] = MueLu::CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COALESCEDROPFACTORY_SHORT
using CoalesceDropFactory [[maybe_unused]] = MueLu::CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COALESCEDROPFACTORY_KOKKOS_SHORT
using CoalesceDropFactory_kokkos [[maybe_unused]] = MueLu::CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COARSEMAPFACTORY_SHORT
using CoarseMapFactory [[maybe_unused]] = MueLu::CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COARSENINGVISUALIZATIONFACTORY_SHORT
using CoarseningVisualizationFactory [[maybe_unused]] = MueLu::CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COMBINEPFACTORY_SHORT
using CombinePFactory [[maybe_unused]] = MueLu::CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CONSTRAINT_SHORT
using Constraint [[maybe_unused]] = MueLu::Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_CONSTRAINTFACTORY_SHORT
using ConstraintFactory [[maybe_unused]] = MueLu::ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COORDINATESTRANSFERFACTORY_SHORT
using CoordinatesTransferFactory [[maybe_unused]] = MueLu::CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_COUPLEDRBMFACTORY_SHORT
using CoupledRBMFactory [[maybe_unused]] = MueLu::CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_DEMOFACTORY_SHORT
using DemoFactory [[maybe_unused]] = MueLu::DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_DIRECTSOLVER_SHORT
using DirectSolver [[maybe_unused]] = MueLu::DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
using DropNegativeEntriesFactory [[maybe_unused]] = MueLu::DropNegativeEntriesFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_EMINPFACTORY_SHORT
using EminPFactory [[maybe_unused]] = MueLu::EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACADEBGS2X2_SHORT
using FacadeBGS2x2 [[maybe_unused]] = MueLu::FacadeBGS2x2<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACADECLASSBASE_SHORT
using FacadeClassBase [[maybe_unused]] = MueLu::FacadeClassBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACADECLASSFACTORY_SHORT
using FacadeClassFactory [[maybe_unused]] = MueLu::FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACADESIMPLE_SHORT
using FacadeSimple [[maybe_unused]] = MueLu::FacadeSimple<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACTORYFACTORY_SHORT
using FactoryFactory [[maybe_unused]] = MueLu::FactoryFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FACTORYMANAGER_SHORT
using FactoryManager [[maybe_unused]] = MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FAKESMOOTHERPROTOTYPE_SHORT
using FakeSmootherPrototype [[maybe_unused]] = MueLu::FakeSmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FILTEREDAFACTORY_SHORT
using FilteredAFactory [[maybe_unused]] = MueLu::FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_FINELEVELINPUTDATAFACTORY_SHORT
using FineLevelInputDataFactory [[maybe_unused]] = MueLu::FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GENERALGEOMETRICPFACTORY_SHORT
using GeneralGeometricPFactory [[maybe_unused]] = MueLu::GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GENERICRFACTORY_SHORT
using GenericRFactory [[maybe_unused]] = MueLu::GenericRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GEOMETRICINTERPOLATIONPFACTORY_SHORT
using GeometricInterpolationPFactory [[maybe_unused]] = MueLu::GeometricInterpolationPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_SHORT
using GeometricInterpolationPFactory_kokkos [[maybe_unused]] = MueLu::GeometricInterpolationPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_GMRESSOLVER_SHORT
using GMRESSolver [[maybe_unused]] = MueLu::GMRESSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_HIERARCHY_SHORT
using Hierarchy [[maybe_unused]] = MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_HIERARCHYMANAGER_SHORT
using HierarchyManager [[maybe_unused]] = MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_HIERARCHYFACTORY_SHORT
using HierarchyFactory [[maybe_unused]] = MueLu::HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_HIERARCHYUTILS_SHORT
using HierarchyUtils [[maybe_unused]] = MueLu::HierarchyUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INTERFACEAGGREGATIONFACTORY_SHORT
using InterfaceAggregationFactory [[maybe_unused]] = MueLu::InterfaceAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_IFPACK2SMOOTHER_SHORT
using Ifpack2Smoother [[maybe_unused]] = MueLu::Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_SHORT
using IndefBlockedDiagonalSmoother [[maybe_unused]] = MueLu::IndefBlockedDiagonalSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INITIALBLOCKNUMBERFACTORY_SHORT
using InitialBlockNumberFactory [[maybe_unused]] = MueLu::InitialBlockNumberFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INTREPIDPCOARSENFACTORY_SHORT
using IntrepidPCoarsenFactory [[maybe_unused]] = MueLu::IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_INVERSEAPPROXIMATIONFACTORY_SHORT
using InverseApproximationFactory [[maybe_unused]] = MueLu::InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LINEDETECTIONFACTORY_SHORT
using LineDetectionFactory [[maybe_unused]] = MueLu::LineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LOCALPERMUTATIONSTRATEGY_SHORT
using LocalPermutationStrategy [[maybe_unused]] = MueLu::LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_LOWPRECISIONFACTORY_SHORT
using LowPrecisionFactory [[maybe_unused]] = MueLu::LowPrecisionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MAPTRANSFERFACTORY_SHORT
using MapTransferFactory [[maybe_unused]] = MueLu::MapTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MATRIXANALYSISFACTORY_SHORT
using MatrixAnalysisFactory [[maybe_unused]] = MueLu::MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MERGEDBLOCKEDMATRIXFACTORY_SHORT
using MergedBlockedMatrixFactory [[maybe_unused]] = MueLu::MergedBlockedMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MERGEDSMOOTHER_SHORT
using MergedSmoother [[maybe_unused]] = MueLu::MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MLPARAMETERLISTINTERPRETER_SHORT
using MLParameterListInterpreter [[maybe_unused]] = MueLu::MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MULTIVECTORTRANSFERFACTORY_SHORT
using MultiVectorTransferFactory [[maybe_unused]] = MueLu::MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_NOTAYAGGREGATIONFACTORY_SHORT
using NotayAggregationFactory [[maybe_unused]] = MueLu::NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_NULLSPACEFACTORY_SHORT
using NullspaceFactory [[maybe_unused]] = MueLu::NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_NULLSPACEFACTORY_KOKKOS_SHORT
using NullspaceFactory_kokkos [[maybe_unused]] = MueLu::NullspaceFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_NULLSPACEPRESMOOTHFACTORY_SHORT
using NullspacePresmoothFactory [[maybe_unused]] = MueLu::NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PARAMETERLISTINTERPRETER_SHORT
using ParameterListInterpreter [[maybe_unused]] = MueLu::ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PATTERNFACTORY_SHORT
using PatternFactory [[maybe_unused]] = MueLu::PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PERFUTILS_SHORT
using PerfUtils [[maybe_unused]] = MueLu::PerfUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PERFMODELS_SHORT
using PerfModels [[maybe_unused]] = MueLu::PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PERMUTATIONFACTORY_SHORT
using PermutationFactory [[maybe_unused]] = MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PERMUTINGSMOOTHER_SHORT
using PermutingSmoother [[maybe_unused]] = MueLu::PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PGPFACTORY_SHORT
using PgPFactory [[maybe_unused]] = MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
using PreDropFunctionBaseClass [[maybe_unused]] = MueLu::PreDropFunctionBaseClass<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
using PreDropFunctionConstVal [[maybe_unused]] = MueLu::PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_PROJECTORSMOOTHER_SHORT
using ProjectorSmoother [[maybe_unused]] = MueLu::ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_RAPFACTORY_SHORT
using RAPFactory [[maybe_unused]] = MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_RAPSHIFTFACTORY_SHORT
using RAPShiftFactory [[maybe_unused]] = MueLu::RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCEACFACTORY_SHORT
using RebalanceAcFactory [[maybe_unused]] = MueLu::RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCEBLOCKACFACTORY_SHORT
using RebalanceBlockAcFactory [[maybe_unused]] = MueLu::RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_SHORT
using RebalanceBlockInterpolationFactory [[maybe_unused]] = MueLu::RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_SHORT
using RebalanceBlockRestrictionFactory [[maybe_unused]] = MueLu::RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REBALANCETRANSFERFACTORY_SHORT
using RebalanceTransferFactory [[maybe_unused]] = MueLu::RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REFMAXWELLSMOOTHER_SHORT
using RefMaxwellSmoother [[maybe_unused]] = MueLu::RefMaxwellSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REGIONRFACTORY_SHORT
using RegionRFactory [[maybe_unused]] = MueLu::RegionRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REGIONRFACTORY_KOKKOS_SHORT
using RegionRFactory_kokkos [[maybe_unused]] = MueLu::RegionRFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REITZINGERPFACTORY_SHORT
using ReitzingerPFactory [[maybe_unused]] = MueLu::ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REORDERBLOCKAFACTORY_SHORT
using ReorderBlockAFactory [[maybe_unused]] = MueLu::ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REPARTITIONFACTORY_SHORT
using RepartitionFactory [[maybe_unused]] = MueLu::RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REPARTITIONBLOCKDIAGONALFACTORY_SHORT
using RepartitionBlockDiagonalFactory [[maybe_unused]] = MueLu::RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REPARTITIONHEURISTICFACTORY_SHORT
using RepartitionHeuristicFactory [[maybe_unused]] = MueLu::RepartitionHeuristicFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REPLICATEPFACTORY_SHORT
using ReplicatePFactory [[maybe_unused]] = MueLu::ReplicatePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_RIGIDBODYMODEFACTORY_SHORT
using RigidBodyModeFactory [[maybe_unused]] = MueLu::RigidBodyModeFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SAPFACTORY_SHORT
using SaPFactory [[maybe_unused]] = MueLu::SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SAPFACTORY_KOKKOS_SHORT
using SaPFactory_kokkos [[maybe_unused]] = MueLu::SaPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SCALEDNULLSPACEFACTORY_SHORT
using ScaledNullspaceFactory [[maybe_unused]] = MueLu::ScaledNullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
using SchurComplementFactory [[maybe_unused]] = MueLu::SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SEGREGATEDAFACTORY_SHORT
using SegregatedAFactory [[maybe_unused]] = MueLu::SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SHIFTEDLAPLACIAN_SHORT
using ShiftedLaplacian [[maybe_unused]] = MueLu::ShiftedLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SHIFTEDLAPLACIANOPERATOR_SHORT
using ShiftedLaplacianOperator [[maybe_unused]] = MueLu::ShiftedLaplacianOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SIMPLESMOOTHER_SHORT
using SimpleSmoother [[maybe_unused]] = MueLu::SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SMOOTHER_SHORT
using Smoother [[maybe_unused]] = MueLu::Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SMOOTHERBASE_SHORT
using SmootherBase [[maybe_unused]] = MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SMOOTHERFACTORY_SHORT
using SmootherFactory [[maybe_unused]] = MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
using SmootherPrototype [[maybe_unused]] = MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SMOOVECCOALESCEDROPFACTORY_SHORT
using SmooVecCoalesceDropFactory [[maybe_unused]] = MueLu::SmooVecCoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SOLVERBASE_SHORT
using SolverBase [[maybe_unused]] = MueLu::SolverBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_STEEPESTDESCENTSOLVER_SHORT
using SteepestDescentSolver [[maybe_unused]] = MueLu::SteepestDescentSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_STRATIMIKOSSMOOTHER_SHORT
using StratimikosSmoother [[maybe_unused]] = MueLu::StratimikosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
using StructuredAggregationFactory [[maybe_unused]] = MueLu::StructuredAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_STRUCTUREDLINEDETECTIONFACTORY_SHORT
using StructuredLineDetectionFactory [[maybe_unused]] = MueLu::StructuredLineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SUBBLOCKAFACTORY_SHORT
using SubBlockAFactory [[maybe_unused]] = MueLu::SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TEKOSMOOTHER_SHORT
using TekoSmoother [[maybe_unused]] = MueLu::TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TENTATIVEPFACTORY_SHORT
using TentativePFactory [[maybe_unused]] = MueLu::TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
using TentativePFactory_kokkos [[maybe_unused]] = MueLu::TentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MATRIXFREETENTATIVEP_SHORT
using MatrixFreeTentativeP [[maybe_unused]] = MueLu::MatrixFreeTentativeP<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MATRIXFREETENTATIVEPFACTORY_SHORT
using MatrixFreeTentativePFactory [[maybe_unused]] = MueLu::MatrixFreeTentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_THRESHOLDAFILTERFACTORY_SHORT
using ThresholdAFilterFactory [[maybe_unused]] = MueLu::ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TOGGLECOORDINATESTRANSFERFACTORY_SHORT
using ToggleCoordinatesTransferFactory [[maybe_unused]] = MueLu::ToggleCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TOGGLEPFACTORY_SHORT
using TogglePFactory [[maybe_unused]] = MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TOPRAPFACTORY_SHORT
using TopRAPFactory [[maybe_unused]] = MueLu::TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TOPSMOOTHERFACTORY_SHORT
using TopSmootherFactory [[maybe_unused]] = MueLu::TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TPETRAOPERATOR_SHORT
using TpetraOperator [[maybe_unused]] = MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TRANSPFACTORY_SHORT
using TransPFactory [[maybe_unused]] = MueLu::TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_RFROMP_OR_TRANSP_SHORT
using RfromP_Or_TransP [[maybe_unused]] = MueLu::RfromP_Or_TransP<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TRILINOSSMOOTHER_SHORT
using TrilinosSmoother [[maybe_unused]] = MueLu::TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UNSMOOSHFACTORY_SHORT
using UnsmooshFactory [[maybe_unused]] = MueLu::UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_USERPFACTORY_SHORT
using UserPFactory [[maybe_unused]] = MueLu::UserPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UTILITIES_SHORT
using Utilities [[maybe_unused]] = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UTILITIESBASE_SHORT
using UtilitiesBase [[maybe_unused]] = MueLu::UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT
using VariableDofLaplacianFactory [[maybe_unused]] = MueLu::VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SEMICOARSENPFACTORY_SHORT
using SemiCoarsenPFactory [[maybe_unused]] = MueLu::SemiCoarsenPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_SEMICOARSENPFACTORY_KOKKOS_SHORT
using SemiCoarsenPFactory_kokkos [[maybe_unused]] = MueLu::SemiCoarsenPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_UZAWASMOOTHER_SHORT
using UzawaSmoother [[maybe_unused]] = MueLu::UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_VISUALIZATIONHELPERS_SHORT
using VisualizationHelpers [[maybe_unused]] = MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ZEROSUBBLOCKAFACTORY_SHORT
using ZeroSubBlockAFactory [[maybe_unused]] = MueLu::ZeroSubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ZOLTANINTERFACE_SHORT
using ZoltanInterface [[maybe_unused]] = MueLu::ZoltanInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_ZOLTAN2INTERFACE_SHORT
using Zoltan2Interface [[maybe_unused]] = MueLu::Zoltan2Interface<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_NODEPARTITIONINTERFACE_SHORT
using NodePartitionInterface [[maybe_unused]] = MueLu::NodePartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_XPETRAOPERATOR_SHORT
using XpetraOperator [[maybe_unused]] = MueLu::XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_REFMAXWELL_SHORT
using RefMaxwell [[maybe_unused]] = MueLu::RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MAXWELL1_SHORT
using Maxwell1 [[maybe_unused]] = MueLu::Maxwell1<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MULTIPHYS_SHORT
using MultiPhys [[maybe_unused]] = MueLu::MultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_MAXWELL_UTILS_SHORT
using Maxwell_Utils [[maybe_unused]] = MueLu::Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#ifdef MUELU_TWOLEVELMATLABFACTORY_SHORT
typedef MueLu::TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TwoLevelMatlabFactory;
#endif
#ifdef MUELU_SINGLELEVELMATLABFACTORY_SHORT
typedef MueLu::SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> SingleLevelMatlabFactory;
#endif
#ifdef MUELU_MATLABSMOOTHER_SHORT
typedef MueLu::MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatlabSmoother;
#endif
