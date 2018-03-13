// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node of the current context.

#include <Xpetra_UseShortNamesScalar.hpp>

#ifdef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
typedef MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> AggregationExportFactory;
#endif
#ifdef MUELU_AMALGAMATIONFACTORY_SHORT
typedef MueLu::AmalgamationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> AmalgamationFactory;
#endif
#ifdef MUELU_AMESOS2SMOOTHER_SHORT
typedef MueLu::Amesos2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> Amesos2Smoother;
#endif
#ifdef MUELU_AMGXOPERATOR_SHORT
typedef MueLu::AMGXOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> AMGXOperator;
#endif
#ifdef MUELU_ALGEBRAICPERMUTATIONSTRATEGY_SHORT
typedef MueLu::AlgebraicPermutationStrategy<Scalar,LocalOrdinal,GlobalOrdinal,Node> AlgebraicPermutationStrategy;
#endif
#ifdef MUELU_BLACKBOXPFACTORY_SHORT
typedef MueLu::BlackBoxPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlackBoxPFactory;
#endif
#ifdef MUELU_BLOCKEDCOARSEMAPFACTORY_SHORT
typedef MueLu::BlockedCoarseMapFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedCoarseMapFactory;
#endif
#ifdef MUELU_BLOCKEDDIRECTSOLVER_SHORT
typedef MueLu::BlockedDirectSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedDirectSolver;
#endif
#ifdef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_SHORT
typedef MueLu::BlockedGaussSeidelSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedGaussSeidelSmoother;
#endif
#ifdef MUELU_BLOCKEDJACOBISMOOTHER_SHORT
typedef MueLu::BlockedJacobiSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedJacobiSmoother;
#endif
#ifdef MUELU_BLOCKEDPFACTORY_SHORT
typedef MueLu::BlockedPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedPFactory;
#endif
#ifdef MUELU_BLOCKEDRAPFACTORY_SHORT
typedef MueLu::BlockedRAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedRAPFactory;
#endif
#ifdef MUELU_BRICKAGGREGATIONFACTORY_SHORT
typedef MueLu::BrickAggregationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BrickAggregationFactory;
#endif
#ifdef MUELU_BRAESSSARAZINSMOOTHER_SHORT
typedef MueLu::BraessSarazinSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> BraessSarazinSmoother;
#endif
#ifdef MUELU_CGSOLVER_SHORT
typedef MueLu::CGSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> CGSolver;
#endif
#ifdef MUELU_CLONEREPARTITIONINTERFACE_SHORT
typedef MueLu::CloneRepartitionInterface<Scalar,LocalOrdinal,GlobalOrdinal,Node> CloneRepartitionInterface;
#endif
#ifdef MUELU_COALESCEDROPFACTORY_SHORT
typedef MueLu::CoalesceDropFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoalesceDropFactory;
#endif
#ifdef MUELU_COALESCEDROPFACTORY_KOKKOS_SHORT
typedef MueLu::CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoalesceDropFactory_kokkos;
#endif
#ifdef MUELU_COARSEMAPFACTORY_SHORT
typedef MueLu::CoarseMapFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoarseMapFactory;
#endif
#ifdef MUELU_COARSEMAPFACTORY_KOKKOS_SHORT
typedef MueLu::CoarseMapFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoarseMapFactory_kokkos;
#endif
#ifdef MUELU_COARSENINGVISUALIZATIONFACTORY_SHORT
typedef MueLu::CoarseningVisualizationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoarseningVisualizationFactory;
#endif
#ifdef MUELU_CONSTRAINT_SHORT
typedef MueLu::Constraint<Scalar,LocalOrdinal,GlobalOrdinal,Node> Constraint;
#endif
#ifdef MUELU_CONSTRAINTFACTORY_SHORT
typedef MueLu::ConstraintFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> ConstraintFactory;
#endif
#ifdef MUELU_COORDINATESTRANSFERFACTORY_SHORT
typedef MueLu::CoordinatesTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoordinatesTransferFactory;
#endif
#ifdef MUELU_COORDINATESTRANSFERFACTORY_KOKKOS_SHORT
typedef MueLu::CoordinatesTransferFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoordinatesTransferFactory_kokkos;
#endif
#ifdef MUELU_COUPLEDRBMFACTORY_SHORT
typedef MueLu::CoupledRBMFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> CoupledRBMFactory;
#endif
#ifdef MUELU_DEMOFACTORY_SHORT
typedef MueLu::DemoFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> DemoFactory;
#endif
#ifdef MUELU_DIRECTSOLVER_SHORT
typedef MueLu::DirectSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> DirectSolver;
#endif
#ifdef MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
typedef MueLu::DropNegativeEntriesFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> DropNegativeEntriesFactory;
#endif
#ifdef MUELU_EMINPFACTORY_SHORT
typedef MueLu::EminPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> EminPFactory;
#endif
#ifdef MUELU_FACADECLASSFACTORY_SHORT
typedef MueLu::FacadeClassFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> FacadeClassFactory;
#endif
#ifdef MUELU_FACTORYMANAGER_SHORT
typedef MueLu::FactoryManager<Scalar,LocalOrdinal,GlobalOrdinal,Node> FactoryManager;
#endif
#ifdef MUELU_FAKESMOOTHERPROTOTYPE_SHORT
typedef MueLu::FakeSmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node> FakeSmootherPrototype;
#endif
#ifdef MUELU_FILTEREDAFACTORY_SHORT
typedef MueLu::FilteredAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> FilteredAFactory;
#endif
#ifdef MUELU_FINELEVELINPUTDATAFACTORY_SHORT
typedef MueLu::FineLevelInputDataFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> FineLevelInputDataFactory;
#endif
#ifdef MUELU_GENERALGEOMETRICPFACTORY_SHORT
typedef MueLu::GeneralGeometricPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> GeneralGeometricPFactory;
#endif
#ifdef MUELU_GENERICRFACTORY_SHORT
typedef MueLu::GenericRFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> GenericRFactory;
#endif
#ifdef MUELU_GMRESSOLVER_SHORT
typedef MueLu::GMRESSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> GMRESSolver;
#endif
#ifdef MUELU_HIERARCHY_SHORT
typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> Hierarchy;
#endif
#ifdef MUELU_HIERARCHYMANAGER_SHORT
typedef MueLu::HierarchyManager<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyManager;
#endif
#ifdef MUELU_HIERARCHYFACTORY_SHORT
typedef MueLu::HierarchyFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyFactory;
#endif
#ifdef MUELU_HIERARCHYUTILS_SHORT
typedef MueLu::HierarchyUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyUtils;
#endif
#ifdef MUELU_IFPACK2SMOOTHER_SHORT
typedef MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> Ifpack2Smoother;
#endif
#ifdef MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_SHORT
typedef MueLu::IndefBlockedDiagonalSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> IndefBlockedDiagonalSmoother;
#endif
#ifdef MUELU_INTREPIDPCOARSENFACTORY_SHORT
typedef MueLu::IntrepidPCoarsenFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> IntrepidPCoarsenFactory;
#endif
#ifdef MUELU_LINEDETECTIONFACTORY_SHORT
typedef MueLu::LineDetectionFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> LineDetectionFactory;
#endif
#ifdef MUELU_LOCALPERMUTATIONSTRATEGY_SHORT
typedef MueLu::LocalPermutationStrategy<Scalar,LocalOrdinal,GlobalOrdinal,Node> LocalPermutationStrategy;
#endif
#ifdef MUELU_MAPTRANSFERFACTORY_SHORT
typedef MueLu::MapTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> MapTransferFactory;
#endif
#ifdef MUELU_MATRIXANALYSISFACTORY_SHORT
typedef MueLu::MatrixAnalysisFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> MatrixAnalysisFactory;
#endif
#ifdef MUELU_MERGEDBLOCKEDMATRIXFACTORY_SHORT
typedef MueLu::MergedBlockedMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> MergedBlockedMatrixFactory;
#endif
#ifdef MUELU_MERGEDSMOOTHER_SHORT
typedef MueLu::MergedSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> MergedSmoother;
#endif
#ifdef MUELU_MULTIVECTORTRANSFERFACTORY_SHORT
typedef MueLu::MultiVectorTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVectorTransferFactory;
#endif
#ifdef MUELU_NULLSPACEFACTORY_SHORT
typedef MueLu::NullspaceFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> NullspaceFactory;
#endif
#ifdef MUELU_NULLSPACEFACTORY_KOKKOS_SHORT
typedef MueLu::NullspaceFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> NullspaceFactory_kokkos;
#endif
#ifdef MUELU_NULLSPACEPRESMOOTHFACTORY_SHORT
typedef MueLu::NullspacePresmoothFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> NullspacePresmoothFactory;
#endif
#ifdef MUELU_PATTERNFACTORY_SHORT
typedef MueLu::PatternFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> PatternFactory;
#endif
#ifdef MUELU_PERFUTILS_SHORT
typedef MueLu::PerfUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node> PerfUtils;
#endif
#ifdef MUELU_PERMUTATIONFACTORY_SHORT
typedef MueLu::PermutationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> PermutationFactory;
#endif
#ifdef MUELU_PERMUTINGSMOOTHER_SHORT
typedef MueLu::PermutingSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> PermutingSmoother;
#endif
#ifdef MUELU_PGPFACTORY_SHORT
typedef MueLu::PgPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> PgPFactory;
#endif
#ifdef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
typedef MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node> PreDropFunctionBaseClass;
#endif
#ifdef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
typedef MueLu::PreDropFunctionConstVal<Scalar,LocalOrdinal,GlobalOrdinal,Node> PreDropFunctionConstVal;
#endif
#ifdef MUELU_PROJECTORSMOOTHER_SHORT
typedef MueLu::ProjectorSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> ProjectorSmoother;
#endif
#ifdef MUELU_RAPFACTORY_SHORT
typedef MueLu::RAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RAPFactory;
#endif
#ifdef MUELU_RAPSHIFTFACTORY_SHORT
typedef MueLu::RAPShiftFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RAPShiftFactory;
#endif
#ifdef MUELU_REBALANCEACFACTORY_SHORT
typedef MueLu::RebalanceAcFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RebalanceAcFactory;
#endif
#ifdef MUELU_REBALANCEBLOCKACFACTORY_SHORT
typedef MueLu::RebalanceBlockAcFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RebalanceBlockAcFactory;
#endif
#ifdef MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_SHORT
typedef MueLu::RebalanceBlockInterpolationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RebalanceBlockInterpolationFactory;
#endif
#ifdef MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_SHORT
typedef MueLu::RebalanceBlockRestrictionFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RebalanceBlockRestrictionFactory;
#endif
#ifdef MUELU_REBALANCETRANSFERFACTORY_SHORT
typedef MueLu::RebalanceTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RebalanceTransferFactory;
#endif
#ifdef MUELU_REORDERBLOCKAFACTORY_SHORT
typedef MueLu::ReorderBlockAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> ReorderBlockAFactory;
#endif
#ifdef MUELU_REPARTITIONFACTORY_SHORT
typedef MueLu::RepartitionFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RepartitionFactory;
#endif
#ifdef MUELU_REPARTITIONHEURISTICFACTORY_SHORT
typedef MueLu::RepartitionHeuristicFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RepartitionHeuristicFactory;
#endif
#ifdef MUELU_RIGIDBODYMODEFACTORY_SHORT
typedef MueLu::RigidBodyModeFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> RigidBodyModeFactory;
#endif
#ifdef MUELU_SAPFACTORY_SHORT
typedef MueLu::SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SaPFactory;
#endif
#ifdef MUELU_SAPFACTORY_KOKKOS_SHORT
typedef MueLu::SaPFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> SaPFactory_kokkos;
#endif
#ifdef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
typedef MueLu::SchurComplementFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SchurComplementFactory;
#endif
#ifdef MUELU_SEGREGATEDAFACTORY_SHORT
typedef MueLu::SegregatedAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SegregatedAFactory;
#endif
#ifdef MUELU_SHIFTEDLAPLACIAN_SHORT
typedef MueLu::ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node> ShiftedLaplacian;
#endif
#ifdef MUELU_SHIFTEDLAPLACIANOPERATOR_SHORT
typedef MueLu::ShiftedLaplacianOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> ShiftedLaplacianOperator;
#endif
#ifdef MUELU_SIMPLESMOOTHER_SHORT
typedef MueLu::SimpleSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> SimpleSmoother;
#endif
#ifdef MUELU_SMOOTHER_SHORT
typedef MueLu::Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> Smoother;
#endif
#ifdef MUELU_SMOOTHERBASE_SHORT
typedef MueLu::SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherBase;
#endif
#ifdef MUELU_SMOOTHERFACTORY_SHORT
typedef MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherFactory;
#endif
#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
typedef MueLu::SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherPrototype;
#endif
#ifdef MUELU_SOLVERBASE_SHORT
typedef MueLu::SolverBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> SolverBase;
#endif
#ifdef MUELU_STEEPESTDESCENTSOLVER_SHORT
typedef MueLu::SteepestDescentSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> SteepestDescentSolver;
#endif
#ifdef MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
typedef MueLu::StructuredAggregationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> StructuredAggregationFactory;
#endif
#ifdef MUELU_STRUCTUREDLINEDETECTIONFACTORY_SHORT
typedef MueLu::StructuredLineDetectionFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> StructuredLineDetectionFactory;
#endif
#ifdef MUELU_SUBBLOCKAFACTORY_SHORT
typedef MueLu::SubBlockAFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SubBlockAFactory;
#endif
#ifdef MUELU_TEKOSMOOTHER_SHORT
typedef MueLu::TekoSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> TekoSmoother;
#endif
#ifdef MUELU_TENTATIVEPFACTORY_SHORT
typedef MueLu::TentativePFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TentativePFactory;
#endif
#ifdef MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
typedef MueLu::TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> TentativePFactory_kokkos;
#endif
#ifdef MUELU_THRESHOLDAFILTERFACTORY_SHORT
typedef MueLu::ThresholdAFilterFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> ThresholdAFilterFactory;
#endif
#ifdef MUELU_TOGGLECOORDINATESTRANSFERFACTORY_SHORT
typedef MueLu::ToggleCoordinatesTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> ToggleCoordinatesTransferFactory;
#endif
#ifdef MUELU_TOGGLEPFACTORY_SHORT
typedef MueLu::TogglePFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TogglePFactory;
#endif
#ifdef MUELU_TOPRAPFACTORY_SHORT
typedef MueLu::TopRAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TopRAPFactory;
#endif
#ifdef MUELU_TOPSMOOTHERFACTORY_SHORT
typedef MueLu::TopSmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TopSmootherFactory;
#endif
#ifdef MUELU_TPETRAOPERATOR_SHORT
typedef MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraOperator;
#endif
#ifdef MUELU_TRANSPFACTORY_SHORT
typedef MueLu::TransPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TransPFactory;
#endif
#ifdef MUELU_TRILINOSSMOOTHER_SHORT
typedef MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> TrilinosSmoother;
#endif
#ifdef MUELU_UNSMOOSHFACTORY_SHORT
typedef MueLu::UnsmooshFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> UnsmooshFactory;
#endif
#ifdef MUELU_USERPFACTORY_SHORT
typedef MueLu::UserPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> UserPFactory;
#endif
#ifdef MUELU_UTILITIES_SHORT
typedef MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node> Utilities;
#endif
#ifdef MUELU_UTILITIESBASE_SHORT
typedef MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> UtilitiesBase;
#endif
#ifdef MUELU_UTILITIES_KOKKOS_SHORT
typedef MueLu::Utilities_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node> Utilities_kokkos;
#endif
#ifdef MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT
typedef MueLu::VariableDofLaplacianFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> VariableDofLaplacianFactory;
#endif
#ifdef MUELU_SEMICOARSENPFACTORY_SHORT
typedef MueLu::SemiCoarsenPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SemiCoarsenPFactory;
#endif
#ifdef MUELU_UZAWASMOOTHER_SHORT
typedef MueLu::UzawaSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> UzawaSmoother;
#endif
#ifdef MUELU_VISUALIZATIONHELPERS_SHORT
typedef MueLu::VisualizationHelpers<Scalar,LocalOrdinal,GlobalOrdinal,Node> VisualizationHelpers;
#endif
#ifdef MUELU_ZOLTANINTERFACE_SHORT
typedef MueLu::ZoltanInterface<Scalar,LocalOrdinal,GlobalOrdinal,Node> ZoltanInterface;
#endif
#ifdef MUELU_ZOLTAN2INTERFACE_SHORT
typedef MueLu::Zoltan2Interface<Scalar,LocalOrdinal,GlobalOrdinal,Node> Zoltan2Interface;
#endif
#ifdef MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_SHORT
typedef MueLu::AdaptiveSaMLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> AdaptiveSaMLParameterListInterpreter;
#endif
#ifdef MUELU_FACTORYFACTORY_SHORT
typedef MueLu::FactoryFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> FactoryFactory;
#endif
#ifdef MUELU_MLPARAMETERLISTINTERPRETER_SHORT
typedef MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> MLParameterListInterpreter;
#endif
#ifdef MUELU_PARAMETERLISTINTERPRETER_SHORT
typedef MueLu::ParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> ParameterListInterpreter;
#endif
#ifdef MUELU_TWOLEVELMATLABFACTORY_SHORT
typedef MueLu::TwoLevelMatlabFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> TwoLevelMatlabFactory;
#endif
#ifdef MUELU_SINGLELEVELMATLABFACTORY_SHORT
typedef MueLu::SingleLevelMatlabFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SingleLevelMatlabFactory;
#endif
#ifdef MUELU_MATLABSMOOTHER_SHORT
typedef MueLu::MatlabSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node> MatlabSmoother;
#endif
