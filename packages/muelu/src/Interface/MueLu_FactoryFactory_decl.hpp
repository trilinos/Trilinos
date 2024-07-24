// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FACTORYFACTORY_DECL_HPP
#define MUELU_FACTORYFACTORY_DECL_HPP

#include <string>
#include <vector>

#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryFactory_fwd.hpp"

#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Hierarchy_fwd.hpp"

#include "MueLu_Monitor.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_AggregateQualityEstimateFactory_fwd.hpp"
#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_BlackBoxPFactory_fwd.hpp"
#include "MueLu_BlockedCoarseMapFactory_fwd.hpp"
#include "MueLu_BlockedCoordinatesTransferFactory_fwd.hpp"
#include "MueLu_BlockedDirectSolver_fwd.hpp"
#include "MueLu_BlockedGaussSeidelSmoother_fwd.hpp"
#include "MueLu_BlockedJacobiSmoother_fwd.hpp"
#include "MueLu_BlockedPFactory_fwd.hpp"
#include "MueLu_BlockedRAPFactory_fwd.hpp"
#include "MueLu_BraessSarazinSmoother_fwd.hpp"
#include "MueLu_BrickAggregationFactory_fwd.hpp"
#include "MueLu_ClassicalMapFactory_fwd.hpp"
#include "MueLu_ClassicalPFactory_fwd.hpp"
#include "MueLu_CloneRepartitionInterface_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_SmooVecCoalesceDropFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_CoarseningVisualizationFactory_fwd.hpp"
#include "MueLu_ConstraintFactory_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_DropNegativeEntriesFactory_fwd.hpp"
#include "MueLu_EminPFactory_fwd.hpp"
#include "MueLu_FilteredAFactory_fwd.hpp"
#include "MueLu_FineLevelInputDataFactory_fwd.hpp"
#include "MueLu_GeneralGeometricPFactory_fwd.hpp"
#include "MueLu_ReplicatePFactory_fwd.hpp"
#include "MueLu_CombinePFactory_fwd.hpp"
#include "MueLu_GenericRFactory_fwd.hpp"
#include "MueLu_GeometricInterpolationPFactory_fwd.hpp"
#include "MueLu_InterfaceAggregationFactory_fwd.hpp"
#include "MueLu_InterfaceMappingTransferFactory_fwd.hpp"
#include "MueLu_InitialBlockNumberFactory_fwd.hpp"
#include "MueLu_IndefBlockedDiagonalSmoother_fwd.hpp"
#include "MueLu_InverseApproximationFactory_fwd.hpp"
#include "MueLu_IsorropiaInterface_fwd.hpp"
#include "MueLu_LineDetectionFactory_fwd.hpp"
#include "MueLu_LocalOrdinalTransferFactory_fwd.hpp"
#include "MueLu_RepartitionInterface_fwd.hpp"
#include "MueLu_RepartitionBlockDiagonalFactory_fwd.hpp"
#include "MueLu_MapTransferFactory_fwd.hpp"
#include "MueLu_MatrixAnalysisFactory_fwd.hpp"
#include "MueLu_MultiVectorTransferFactory_fwd.hpp"
#include "MueLu_NotayAggregationFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_NullspacePresmoothFactory_fwd.hpp"
#include "MueLu_PatternFactory_fwd.hpp"
#include "MueLu_PgPFactory_fwd.hpp"
#include "MueLu_RebalanceBlockInterpolationFactory_fwd.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory_fwd.hpp"
#include "MueLu_RebalanceBlockAcFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"
#include "MueLu_RegionRFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_RepartitionHeuristicFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_RAPShiftFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_ReorderBlockAFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_ScaledNullspaceFactory_fwd.hpp"
#include "MueLu_SegregatedAFactory_fwd.hpp"
#include "MueLu_SemiCoarsenPFactory_fwd.hpp"
#include "MueLu_SchurComplementFactory_fwd.hpp"
#include "MueLu_SimpleSmoother_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_StructuredAggregationFactory_fwd.hpp"
#include "MueLu_StructuredLineDetectionFactory_fwd.hpp"
#include "MueLu_SubBlockAFactory_fwd.hpp"
#ifdef HAVE_MUELU_TEKO
#include "MueLu_TekoSmoother_fwd.hpp"
#endif
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_ToggleCoordinatesTransferFactory_fwd.hpp"
#include "MueLu_TogglePFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_RfromP_Or_TransP_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_HybridAggregationFactory_fwd.hpp"
#include "MueLu_UnsmooshFactory_fwd.hpp"
#include "MueLu_UserAggregationFactory_fwd.hpp"
#include "MueLu_UserPFactory_fwd.hpp"
#include "MueLu_UzawaSmoother_fwd.hpp"
#include "MueLu_VariableDofLaplacianFactory_fwd.hpp"
#include "MueLu_ZeroSubBlockAFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_Zoltan2Interface_fwd.hpp"
#include "MueLu_NodePartitionInterface_fwd.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"
#include "MueLu_GeometricInterpolationPFactory_kokkos_fwd.hpp"
#ifdef HAVE_MUELU_DEPRECATED_CODE
#include "MueLu_NullspaceFactory_kokkos_fwd.hpp"
#include "MueLu_SaPFactory_kokkos_fwd.hpp"
#endif
#include "MueLu_SemiCoarsenPFactory_kokkos_fwd.hpp"
#include "MueLu_StructuredAggregationFactory_kokkos_fwd.hpp"
#include "MueLu_TentativePFactory_kokkos_fwd.hpp"
#include "MueLu_MatrixFreeTentativePFactory_fwd.hpp"
#include "MueLu_RegionRFactory_kokkos_fwd.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory_fwd.hpp"
#endif

namespace MueLu {

/*! class FactoryFactory

@brief Factory that can generate other factories from


*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FactoryFactory : public BaseClass {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap;  // TODO: remove
  typedef std::map<std::string, RCP<FactoryManagerBase> > FactoryManagerMap;

 public:
  /// \brief: Interpret Factory parameter list and build new factory
  ///
  /// \param param [in]: ParameterEntry being either the parameter list containing the "factory" parameter declaring the factory type (e.g., "TrilinosSmoother") or being a plain Parameter containing the factory type as value
  /// \param factoryMapIn [in]: FactoryMap containing a map between factory name (e.g., "smootherFact1") and corresponding factory of all previously defined factories
  /// \param factoryManagersIn [in]: FactoryManagerMap containing a map between group names and Factory manager objects. Needed for factories with sub-factory managers.
  ///
  /// Parameter List Parsing:
  /// ---------
  ///     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
  ///
  /// or:
  ///
  ///     <ParameterList name="smootherFact1">
  ///       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
  ///       ...
  ///     </ParameterList>
  ///
  virtual RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry& param, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  //
  //
  //

  // FOLLOWING FUNCTIONS SHOULD LIVE WITH THE CORRESPONDING CLASS

  //
  //
  //

  template <class T>  // T must implement the Factory interface
  RCP<T> Build(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  template <class T>  // T must implement the Factory interface
  RCP<T> Build2(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  template <class T>  // T must implement the Factory interface
  RCP<T> BuildRAPFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  template <class T>  // T must implement the Factory interface
  RCP<T> BuildTogglePFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  RCP<ToggleCoordinatesTransferFactory> BuildToggleCoordinatesTransferFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  //! TrilinosSmoother
  // Parameter List Parsing:
  //     <ParameterList name="smootherFact1">
  //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
  //       <Parameter name="verbose" type="string" value="Warnings"/>
  //       <Parameter name="type" type="string" value="RELAXATION"/>
  //       <ParameterList name="ParameterList">
  //       ...
  //       </ParameterList>
  //     </ParameterList>
  RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

#ifdef HAVE_MUELU_MATLAB
  //! MatlabSmoother
  // Parameter List Parsing:
  //     <ParameterList name="smootherFact1">
  //       <Parameter name="factory" type="string" value="MatlabSmoother"/>
  //       <Parameter name="Setup Function" type="string" value="mySmootherSetup.m"/>
  //       <Parameter name="Solve Function" type="string" value="mySmootherSolve.m"/>
  //       <!--A is implicitly included in this list and nothing else is needed to get diagonal-->
  //       <Parameter name="Needs" type="string" value=""/>
  //       <!--A,x,b are also assumed inputs to the solver: only one additional arg then (diag)-->
  //       <Parameter name="Number of Solver Args"               type="int" value="1"/>
  //     </ParameterList>
  RCP<FactoryBase> BuildMatlabSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;
#endif

  RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList& paramList, const FactoryMap& /* factoryMapIn */, const FactoryManagerMap& /* factoryManagersIn */) const;

  template <class T>  // T must implement the Factory interface
  RCP<FactoryBase> BuildBlockedSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

#ifdef HAVE_MUELU_TEKO
  RCP<FactoryBase> BuildTekoSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;
#endif

  RCP<FactoryBase> BuildBlockedDirectSolver(const Teuchos::ParameterList& paramList, const FactoryMap& /* factoryMapIn */, const FactoryManagerMap& /* factoryManagersIn */) const;

  // RCP<FactoryBase> BuildBlockedPFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
  //   RCP<BlockedPFactory> pfac = rcp(new BlockedPFactory());

  template <class T>  // T must implement the Factory interface
  RCP<T> BuildBlockedFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

  template <class T>  // T must implement the Factory interface
  RCP<T> BuildBlockedCoordFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const;

};  // class
}  // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif  // MUELU_FACTORYFACTORY_DECL_HPP

// TODO: handle factory parameters
// TODO: parameter validator
// TODO: static
// TODO: default parameters should not be duplicated here and on the Factory (ex: default for overlap (=0) is defined both here and on TrilinosSmoother constructors)
