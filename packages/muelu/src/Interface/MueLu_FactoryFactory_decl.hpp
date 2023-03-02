// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

#include "MueLu_AggregateQualityEstimateFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_BlackBoxPFactory.hpp"
#include "MueLu_BlockedCoarseMapFactory.hpp"
#include "MueLu_BlockedCoordinatesTransferFactory.hpp"
#include "MueLu_BlockedDirectSolver.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_BlockedJacobiSmoother.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"
#include "MueLu_BrickAggregationFactory.hpp"
#include "MueLu_ClassicalMapFactory.hpp"
#include "MueLu_ClassicalPFactory.hpp"
#include "MueLu_CloneRepartitionInterface.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_SmooVecCoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CoarseningVisualizationFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_DropNegativeEntriesFactory.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_FineLevelInputDataFactory.hpp"
#include "MueLu_GeneralGeometricPFactory.hpp"
#include "MueLu_ReplicatePFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_GeometricInterpolationPFactory.hpp"
#include "MueLu_InterfaceAggregationFactory.hpp"
#include "MueLu_InterfaceMappingTransferFactory.hpp"
#include "MueLu_InitialBlockNumberFactory.hpp"
#include "MueLu_IndefBlockedDiagonalSmoother.hpp"
#include "MueLu_InverseApproximationFactory.hpp"
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_LineDetectionFactory.hpp"
#include "MueLu_LocalOrdinalTransferFactory.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_RepartitionBlockDiagonalFactory.hpp"
#include "MueLu_MapTransferFactory.hpp"
#include "MueLu_MatrixAnalysisFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_NotayAggregationFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_NullspacePresmoothFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_RebalanceBlockInterpolationFactory.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RegionRFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RAPShiftFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_ReorderBlockAFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_ScaledNullspaceFactory.hpp"
#include "MueLu_SegregatedAFactory.hpp"
#include "MueLu_SemiCoarsenPFactory.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_SimpleSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_StructuredAggregationFactory.hpp"
#include "MueLu_StructuredLineDetectionFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#ifdef HAVE_MUELU_TEKO
#include "MueLu_TekoSmoother.hpp"
#endif
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_ToggleCoordinatesTransferFactory.hpp"
#include "MueLu_TogglePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RfromP_Or_TransP.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_HybridAggregationFactory.hpp"
#include "MueLu_UnsmooshFactory.hpp"
#include "MueLu_UserAggregationFactory.hpp"
#include "MueLu_UserPFactory.hpp"
#include "MueLu_UzawaSmoother.hpp"
#include "MueLu_VariableDofLaplacianFactory.hpp"
#include "MueLu_ZeroSubBlockAFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_NodePartitionInterface.hpp"


#include "MueLu_AmalgamationFactory_kokkos.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
#include "MueLu_GeometricInterpolationPFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_SaPFactory_kokkos.hpp"
#include "MueLu_SemiCoarsenPFactory_kokkos.hpp"
#include "MueLu_StructuredAggregationFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_MatrixFreeTentativePFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"
#include "MueLu_RegionRFactory_kokkos.hpp"

#ifdef HAVE_MUELU_MATLAB
// This is distasteful, but (sadly) neccesary due to peculiarities in MueLu's build system.
#include "../matlab/src/MueLu_SingleLevelMatlabFactory_decl.hpp"
#include "../matlab/src/MueLu_SingleLevelMatlabFactory_def.hpp"
#include "../matlab/src/MueLu_TwoLevelMatlabFactory_decl.hpp"
#include "../matlab/src/MueLu_TwoLevelMatlabFactory_def.hpp"
#include "../matlab/src/MueLu_MatlabSmoother_decl.hpp"
#include "../matlab/src/MueLu_MatlabSmoother_def.hpp"
#endif

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#endif

namespace MueLu {

  /*! class FactoryFactory

  @brief Factory that can generate other factories from


  */
  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class FactoryFactory : public BaseClass {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase>  > FactoryMap; // TODO: remove
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
    virtual RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry& param, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      // Find factory
      std::string            factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList()) {
        factoryName = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      }

      // TODO: see how Teko handles this (=> register factories).
      if (factoryName == "AggregateQualityEstimateFactory")       return Build2<AggregateQualityEstimateFactory>       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "AggregationExportFactory")              return Build2<AggregationExportFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "AmalgamationFactory")                   return Build2<AmalgamationFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedCoarseMapFactory")               return Build2<BlockedCoarseMapFactory>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedRAPFactory")                     return BuildRAPFactory<BlockedRAPFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BrickAggregationFactory")               return Build2<BrickAggregationFactory>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ClassicalMapFactory")                   return Build2<ClassicalMapFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ClassicalPFactory")                     return Build2<ClassicalPFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CloneRepartitionInterface")             return Build2<CloneRepartitionInterface>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoarseMapFactory")                      return Build2<CoarseMapFactory>                      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoarseningVisualizationFactory")        return Build2<CoarseningVisualizationFactory>        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoalesceDropFactory")                   return Build2<CoalesceDropFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SmooVecCoalesceDropFactory")           return Build2<SmooVecCoalesceDropFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ConstraintFactory")                     return Build2<ConstraintFactory>                     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoupledAggregationFactory")             return BuildCoupledAggregationFactory                (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoordinatesTransferFactory")            return Build2<CoordinatesTransferFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "DirectSolver")                          return BuildDirectSolver                             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "DropNegativeEntriesFactory")            return Build2<DropNegativeEntriesFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "EminPFactory")                          return Build2<EminPFactory>                          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "FilteredAFactory")                      return Build2<FilteredAFactory>                      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "FineLevelInputDataFactory")             return Build2<FineLevelInputDataFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "GeneralGeometricPFactory")              return Build2<GeneralGeometricPFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ReplicatePFactory")                     return Build2<ReplicatePFactory>                     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "GenericRFactory")                       return Build2<GenericRFactory>                       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "GeometricInterpolationPFactory")        return Build2<GeometricInterpolationPFactory>        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "HybridAggregationFactory")              return Build2<HybridAggregationFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "InterfaceAggregationFactory")           return Build2<InterfaceAggregationFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "InterfaceMappingTransferFactory")       return Build2<InterfaceMappingTransferFactory>       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "InverseApproximationFactory")           return Build2<InverseApproximationFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "InitialBlockNumberFactory")              return Build2<InitialBlockNumberFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "LineDetectionFactory")                  return Build2<LineDetectionFactory>                  (paramList, factoryMapIn, factoryManagersIn);
      // LocalOrdinalTransferFactory is a utility factory that can be used for multiple things, so there is no default
      //      if (factoryName == "LocalOrdinalTransferFactory")           return Build2<LocalOrdinalTransferFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MapTransferFactory")                    return Build2<MapTransferFactory>                    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MatrixAnalysisFactory")                 return Build2<MatrixAnalysisFactory>                 (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MultiVectorTransferFactory")            return Build2<MultiVectorTransferFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NoFactory")                             return MueLu::NoFactory::getRCP();
      if (factoryName == "NoSmoother")                            return rcp(new SmootherFactory(Teuchos::null));
      if (factoryName == "NotayAggregationFactory")               return Build2<NotayAggregationFactory>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NullspaceFactory")                      return Build2<NullspaceFactory>                      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NullspacePresmoothFactory")             return Build2<NullspacePresmoothFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "PatternFactory")                        return Build2<PatternFactory>                        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "PgPFactory")                            return Build2<PgPFactory>                            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SaPFactory")                            return Build2<SaPFactory>                            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RAPFactory")                            return BuildRAPFactory<RAPFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RAPShiftFactory")                       return BuildRAPFactory<RAPShiftFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceAcFactory")                    return Build2<RebalanceAcFactory>                    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceTransferFactory")              return Build2<RebalanceTransferFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RegionRFactory")                        return Build2<RegionRFactory>                        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RegionRFactory_kokkos")                 return Build2<RegionRFactory_kokkos>                 (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ReorderBlockAFactory")                  return Build2<ReorderBlockAFactory>                  (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RepartitionInterface")                  return Build2<RepartitionInterface>                  (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ScaledNullspaceFactory")                return Build2<ScaledNullspaceFactory>                (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SegregatedAFactory")                    return Build2<SegregatedAFactory>                    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SemiCoarsenPFactory")                   return Build2<SemiCoarsenPFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "StructuredAggregationFactory")          return Build2<StructuredAggregationFactory>          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "StructuredLineDetectionFactory")        return Build2<StructuredLineDetectionFactory>        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SubBlockAFactory")                      return Build2<SubBlockAFactory>                      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TentativePFactory")                     return Build2<TentativePFactory>                     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ToggleCoordinatesTransferFactory")      return BuildToggleCoordinatesTransferFactory         (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TogglePFactory")                        return BuildTogglePFactory<TogglePFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TransPFactory")                         return Build2<TransPFactory>                         (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RfromP_Or_TransP")                      return Build2<RfromP_Or_TransP>                      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TrilinosSmoother")                      return BuildTrilinosSmoother                         (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UncoupledAggregationFactory")           return Build2<UncoupledAggregationFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UnsmooshFactory")                       return Build2<UnsmooshFactory>                       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UserAggregationFactory")                return Build2<UserAggregationFactory>                (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UserPFactory")                          return Build2<UserPFactory>                          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "VariableDofLaplacianFactory")           return Build2<VariableDofLaplacianFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ZeroSubBlockAFactory")                  return Build2<ZeroSubBlockAFactory>                  (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "AmalgamationFactory_kokkos")            return Build2<AmalgamationFactory_kokkos>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoalesceDropFactory_kokkos")            return Build2<CoalesceDropFactory_kokkos>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoarseMapFactory_kokkos")               return Build2<CoarseMapFactory_kokkos>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoordinatesTransferFactory_kokkos")     return Build2<CoordinatesTransferFactory_kokkos>     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "GeometricInterpolationPFactory_kokkos") return Build2<GeometricInterpolationPFactory_kokkos> (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NullspaceFactory_kokkos")               return Build2<NullspaceFactory_kokkos>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SaPFactory_kokkos")                     return Build2<SaPFactory_kokkos>                     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SemiCoarsenPFactory_kokkos")            return Build2<SemiCoarsenPFactory_kokkos>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "StructuredAggregationFactory_kokkos")   return Build2<StructuredAggregationFactory_kokkos>   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TentativePFactory_kokkos")              return Build2<TentativePFactory_kokkos>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MatrixFreeTentativePFactory_kokkos")    return Build2<MatrixFreeTentativePFactory_kokkos>    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UncoupledAggregationFactory_kokkos")    return Build2<UncoupledAggregationFactory_kokkos>    (paramList, factoryMapIn, factoryManagersIn);

      if (factoryName == "ZoltanInterface") {
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
        return Build2<ZoltanInterface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a ZoltanInterface object: Zoltan is disabled: HAVE_MUELU_ZOLTAN && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN && HAVE_MPI
      }
      if (factoryName == "Zoltan2Interface") {
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
        return Build2<Zoltan2Interface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a Zoltan2Interface object: Zoltan2 is disabled: HAVE_MUELU_ZOLTAN2 && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN2 && HAVE_MPI
      }
      if (factoryName == "IsorropiaInterface") {
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
        return Build2<IsorropiaInterface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a IsorropiaInterface object: Isorropia is disabled: HAVE_MUELU_ISORROPIA && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN2 && HAVE_MPI
      }

      if (factoryName == "NodePartitionInterface") {
#if defined(HAVE_MPI)
        return Build2<NodePartitionInterface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a NodePartitionInterface object: HAVE_MPI == false.");
#endif // HAVE_MPI
      }

      if (factoryName == "RepartitionFactory") {
#ifdef HAVE_MPI
        return Build2<RepartitionFactory>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a RepartitionFactory object: HAVE_MPI == false.");
#endif // HAVE_MPI
      }
      if (factoryName == "RepartitionHeuristicFactory") {
#ifdef HAVE_MPI
        return Build2<RepartitionHeuristicFactory>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a RepartitionHeuristicFactory object: HAVE_MPI == false.");
#endif // HAVE_MPI
      }
      // Blocked factories
      if (factoryName == "BlockedCoordinatesTransferFactory")     return BuildBlockedCoordFactory<BlockedCoordinatesTransferFactory>     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedDirectSolver")             return BuildBlockedDirectSolver(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedGaussSeidelSmoother")      return BuildBlockedSmoother<BlockedGaussSeidelSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedJacobiSmoother")           return BuildBlockedSmoother<BlockedJacobiSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedPFactory")                 return BuildBlockedFactory<BlockedPFactory>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BraessSarazinSmoother")           return BuildBlockedSmoother<BraessSarazinSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "IndefiniteBlockDiagonalSmoother") return BuildBlockedSmoother<IndefBlockedDiagonalSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SimpleSmoother")                  return BuildBlockedSmoother<SimpleSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SchurComplementFactory")          return Build2<SchurComplementFactory> (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceBlockRestrictionFactory")return BuildBlockedFactory<RebalanceBlockRestrictionFactory>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceBlockAcFactory")         return BuildBlockedFactory<RebalanceBlockAcFactory>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceBlockInterpolationFactory") return BuildBlockedFactory<RebalanceBlockInterpolationFactory>(paramList, factoryMapIn, factoryManagersIn);
#ifdef HAVE_MPI
      if (factoryName == "RepartitionBlockDiagonalFactory")    return Build2<RepartitionBlockDiagonalFactory>    (paramList, factoryMapIn, factoryManagersIn);
#endif
#ifdef HAVE_MUELU_TEKO
      if (factoryName == "TekoSmoother")                    return BuildTekoSmoother(paramList, factoryMapIn, factoryManagersIn);
#endif
      if (factoryName == "UzawaSmoother")                   return BuildBlockedSmoother<UzawaSmoother>(paramList, factoryMapIn, factoryManagersIn);

      // Matlab factories
#ifdef HAVE_MUELU_MATLAB
      if (factoryName == "TwoLevelMatlabFactory")           return Build2<TwoLevelMatlabFactory>        (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SingleLevelMatlabFactory")        return Build2<SingleLevelMatlabFactory>     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MatlabSmoother")                  return BuildMatlabSmoother                  (paramList, factoryMapIn, factoryManagersIn);
#endif

#ifdef HAVE_MUELU_INTREPID2
      if (factoryName == "IntrepidPCoarsenFactory")           return Build2<IntrepidPCoarsenFactory>        (paramList, factoryMapIn, factoryManagersIn);
#endif

      // Use a user defined factories (in <Factories> node)
      if (factoryMapIn.find(factoryName) != factoryMapIn.end()) {
        TEUCHOS_TEST_FOR_EXCEPTION((param.isList() && (++paramList.begin() != paramList.end())), Exceptions::RuntimeError,
                                   "MueLu::FactoryFactory: Error during the parsing of: " << std::endl << paramList << std::endl
                                   << "'" << factoryName << "' is not a factory name but an existing instance of a factory." << std::endl
                                   << "Extra parameters cannot be specified after the creation of the object." << std::endl << std::endl
                                   << "Correct syntaxes includes:" << std::endl
                                   << " <Parameter name=\"...\" type=\"string\" value=\"" << factoryName << "\"/>" << std::endl
                                   << "or" << std::endl
                                   << " <ParameterList name=\"...\"><Parameter name=\"factory\" type=\"string\" value=\"" << factoryName << "\"/></ParameterList>" << std::endl
                                   );

        return factoryMapIn.find(factoryName)->second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory: unknown factory name : " << factoryName);

      TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
    }

    //
    //
    //

    // FOLLOWING FUNCTIONS SHOULD LIVE WITH THE CORRESPONDING CLASS

    //
    //
    //

#define arraysize(ar)  (sizeof(ar) / sizeof(ar[0]))

    template <class T> // T must implement the Factory interface
    RCP<T> Build(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory = rcp(new T());

      const char* strarray[] = {"A", "P", "R", "Graph", "UnAmalgamationInfo", "Aggregates", "Nullspace", "TransferFactory", "DofsPerNode"};
      std::vector<std::string> v(strarray, strarray + arraysize(strarray));
      for (size_t i = 0; i < v.size(); ++i)
        if (paramList.isParameter(v[i]))
          factory->SetFactory(v[i], BuildFactory(paramList.getEntry(v[i]), factoryMapIn, factoryManagersIn));

      return factory;
    }

    template <class T> // T must implement the Factory interface
    RCP<T> Build2(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory = rcp(new T());

      ParameterList paramListWithFactories;

      // Read the RCP<Factory> parameters of the class T
      RCP<const ParameterList> validParamList = factory->GetValidParameterList(); // TODO check for Teuchos::null (no parameter list validation)
      TEUCHOS_TEST_FOR_EXCEPTION(validParamList == Teuchos::null, Exceptions::RuntimeError, "FactoryFactory::Build2: default parameter list is null. Please fix this.");
      for (ParameterList::ConstIterator param = validParamList->begin(); param != validParamList->end(); ++param) {
        const std::string& pName = validParamList->name(param);

        if (!paramList.isParameter(pName)) {
          // Ignore unknown parameters
          continue;
        }

        if (validParamList->isType< RCP<const FactoryBase> >(pName)) {
          // Generate or get factory described by param
          RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry(pName), factoryMapIn, factoryManagersIn);
          paramListWithFactories.set(pName, generatingFact);
        } else if (validParamList->isType<RCP<const ParameterList> >(pName)) {
          if (pName == "ParameterList") {
            // NOTE: we cannot use
            //     subList = sublist(rcpFromRef(paramList), pName)
            // here as that would result in sublist also being a reference to a temporary object.
            // The resulting dereferencing in the corresponding factory would then segfault
            RCP<const ParameterList> subList = Teuchos::sublist(rcp(new ParameterList(paramList)), pName);
            paramListWithFactories.set(pName, subList);
          }
        } else {
          paramListWithFactories.setEntry(pName, paramList.getEntry(pName));
        }
      }

      // Configure the factory
      factory->SetParameterList(paramListWithFactories);

      return factory;
    }

    template <class T> // T must implement the Factory interface
    RCP<T> BuildRAPFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory;
      if (paramList.isSublist("TransferFactories") == false) {
        factory = Build2<T>(paramList, factoryMapIn, factoryManagersIn);

      } else {
        RCP<Teuchos::ParameterList>       paramListNonConst = rcp(new Teuchos::ParameterList(paramList));
        RCP<const Teuchos::ParameterList> transferFactories = rcp(new Teuchos::ParameterList(*sublist(paramListNonConst, "TransferFactories")));

        paramListNonConst->remove("TransferFactories");

        factory = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

        for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
          RCP<const FactoryBase> p = BuildFactory(transferFactories->entry(param), factoryMapIn, factoryManagersIn);
          factory->AddTransferFactory(p);
        }
      }

      return factory;
    }

    template <class T> // T must implement the Factory interface
    RCP<T> BuildTogglePFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory;
      if (paramList.isSublist("TransferFactories") == false) {
          //TODO put in an error message: the TogglePFactory needs a TransferFactories sublist!
        factory = Build2<T>(paramList, factoryMapIn, factoryManagersIn);

      } else {
        RCP<Teuchos::ParameterList>       paramListNonConst = rcp(new Teuchos::ParameterList(paramList));
        RCP<const Teuchos::ParameterList> transferFactories = rcp(new Teuchos::ParameterList(*sublist(paramListNonConst, "TransferFactories")));

        paramListNonConst->remove("TransferFactories");

        // build TogglePFactory
        factory = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

        // count how many prolongation factories and how many coarse null space factories have been declared.
        // the numbers must match!
        int numProlongatorFactories = 0;
        int numPtentFactories = 0;
        int numCoarseNspFactories   = 0;
        for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
          size_t foundNsp = transferFactories->name(param).find("Nullspace");
          if (foundNsp != std::string::npos && foundNsp == 0 && transferFactories->name(param).length()==10) {
            numCoarseNspFactories++;
            continue;
          }
          size_t foundPtent   = transferFactories->name(param).find("Ptent");
          if (foundPtent != std::string::npos && foundPtent == 0 && transferFactories->name(param).length()==6) {
            numPtentFactories++;
            continue;
          }
          size_t foundP   = transferFactories->name(param).find("P");
          if (foundP != std::string::npos && foundP == 0 && transferFactories->name(param).length()==2) {
            numProlongatorFactories++;
            continue;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(numProlongatorFactories!=numCoarseNspFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: The user has to provide the same number of prolongator and coarse nullspace factories!");
        TEUCHOS_TEST_FOR_EXCEPTION(numPtentFactories!=numCoarseNspFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: The user has to provide the same number of ptent and coarse nullspace factories!");
        TEUCHOS_TEST_FOR_EXCEPTION(numProlongatorFactories < 2, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: The TogglePFactory needs at least two different prolongation operators. The factories have to be provided using the names P%i and Nullspace %i, where %i denotes a number between 1 and 9.");

        // create empty vectors with data
        std::vector<Teuchos::ParameterEntry> prolongatorFactoryNames(numProlongatorFactories);
        std::vector<Teuchos::ParameterEntry> coarseNspFactoryNames(numProlongatorFactories);
        std::vector<Teuchos::ParameterEntry> ptentFactoryNames(numProlongatorFactories);

        for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
          size_t foundNsp = transferFactories->name(param).find("Nullspace");
          if (foundNsp != std::string::npos && foundNsp == 0 && transferFactories->name(param).length()==10) {
            int number = atoi(&(transferFactories->name(param).at(9)));
                TEUCHOS_TEST_FOR_EXCEPTION(number < 1 || number > numProlongatorFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: Please use the format Nullspace%i with %i an integer between 1 and the maximum number of prolongation operators in TogglePFactory!");
                coarseNspFactoryNames[number-1] = transferFactories->entry(param);
                continue;
          }
          size_t foundPtent   = transferFactories->name(param).find("Ptent");
          if (foundPtent != std::string::npos && foundPtent == 0 && transferFactories->name(param).length()==6) {
            int number = atoi(&(transferFactories->name(param).at(5)));
                TEUCHOS_TEST_FOR_EXCEPTION(number < 1 || number > numPtentFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: Please use the format Ptent%i with %i an integer between 1 and the maximum number of prolongation operators in TogglePFactory!");
                ptentFactoryNames[number-1] = transferFactories->entry(param);
                continue;
          }
          size_t foundP   = transferFactories->name(param).find("P");
          if (foundP != std::string::npos && foundP == 0 && transferFactories->name(param).length()==2) {
            int number = atoi(&(transferFactories->name(param).at(1)));
                TEUCHOS_TEST_FOR_EXCEPTION(number < 1 || number > numProlongatorFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleP: Please use the format P%i with %i an integer between 1 and the maximum number of prolongation operators in TogglePFactory!");
                prolongatorFactoryNames[number-1] = transferFactories->entry(param);
                continue;
          }
        }

        // register all prolongation factories in TogglePFactory
        for (std::vector<Teuchos::ParameterEntry>::const_iterator it = prolongatorFactoryNames.begin(); it != prolongatorFactoryNames.end(); ++it) {
          RCP<const FactoryBase> p = BuildFactory(*it, factoryMapIn, factoryManagersIn);
          factory->AddProlongatorFactory(p);
        }

        // register all tentative prolongation factories in TogglePFactory
        for (std::vector<Teuchos::ParameterEntry>::const_iterator it = ptentFactoryNames.begin(); it != ptentFactoryNames.end(); ++it) {
          RCP<const FactoryBase> p = BuildFactory(*it, factoryMapIn, factoryManagersIn);
          factory->AddPtentFactory(p);
        }

        // register all coarse nullspace factories in TogglePFactory
        for (std::vector<Teuchos::ParameterEntry>::const_iterator it = coarseNspFactoryNames.begin(); it != coarseNspFactoryNames.end(); ++it) {
          RCP<const FactoryBase> p = BuildFactory(*it, factoryMapIn, factoryManagersIn);
          factory->AddCoarseNullspaceFactory(p);
        }
      }
      return factory;
    }

    RCP<ToggleCoordinatesTransferFactory> BuildToggleCoordinatesTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<ToggleCoordinatesTransferFactory> factory;
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("TransferFactories") == false, Exceptions::RuntimeError, "FactoryFactory::BuildToggleCoordinatesTransferFactory: the ToggleCoordinatesTransferFactory needs a sublist 'TransferFactories' containing information about the subfactories for coordinate transfer!");

      RCP<Teuchos::ParameterList>       paramListNonConst = rcp(new Teuchos::ParameterList(paramList));
      RCP<const Teuchos::ParameterList> transferFactories = rcp(new Teuchos::ParameterList(*sublist(paramListNonConst, "TransferFactories")));
      paramListNonConst->remove("TransferFactories");

      // build CoordinatesTransferFactory
      factory = Build2<ToggleCoordinatesTransferFactory>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // count how many coordinate transfer factories have been declared.
      // the numbers must match!
      int numCoordTransferFactories = 0;
      for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
        size_t foundCoordinates = transferFactories->name(param).find("Coordinates");
        if (foundCoordinates != std::string::npos && foundCoordinates == 0 && transferFactories->name(param).length()==12) {
          numCoordTransferFactories++;
          continue;
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(numCoordTransferFactories != 2, Exceptions::RuntimeError, "FactoryFactory::BuildToggleCoordinatesTransfer: The ToggleCoordinatesTransferFactory needs two (different) coordinate transfer factories. The factories have to be provided using the names Coordinates%i, where %i denotes a number between 1 and 9.");

      // create empty vectors with data
      std::vector<Teuchos::ParameterEntry> coarseCoordsFactoryNames(numCoordTransferFactories);

      for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
        size_t foundCoords = transferFactories->name(param).find("Coordinates");
        if (foundCoords != std::string::npos && foundCoords == 0 && transferFactories->name(param).length()==12) {
          int number = atoi(&(transferFactories->name(param).at(11)));
              TEUCHOS_TEST_FOR_EXCEPTION(number < 1 || number > numCoordTransferFactories, Exceptions::RuntimeError, "FactoryFactory::BuildToggleCoordinatesTransfer: Please use the format Coordinates%i with %i an integer between 1 and the maximum number of coordinate transfer factories in ToggleCoordinatesTransferFactory!");
              coarseCoordsFactoryNames[number-1] = transferFactories->entry(param);
              continue;
        }
      }

      // register all coarse nullspace factories in TogglePFactory
      for (std::vector<Teuchos::ParameterEntry>::const_iterator it = coarseCoordsFactoryNames.begin(); it != coarseCoordsFactoryNames.end(); ++it) {
        RCP<const FactoryBase> p = BuildFactory(*it, factoryMapIn, factoryManagersIn);
        factory->AddCoordTransferFactory(p);
      }

      return factory;
    }

    //! CoupledAggregationFactory
    RCP<FactoryBase> BuildCoupledAggregationFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<CoupledAggregationFactory> factory = Build<CoupledAggregationFactory>(paramList, factoryMapIn, factoryManagersIn);

      if (paramList.isParameter("aggregation: ordering"))
        factory->SetOrdering(paramList.get<std::string>("aggregation: ordering"));

      if (paramList.isParameter("aggregation: max selected neighbors"))
        factory->SetMaxNeighAlreadySelected(paramList.get<int>("aggregation: max selected neighbors"));

      if (paramList.isParameter("Phase3AggCreation"))
        factory->SetPhase3AggCreation(paramList.get<double>("Phase3AggCreation"));

      if(paramList.isParameter("aggregation: min agg size"))
        factory->SetMinNodesPerAggregate(paramList.get<int>("aggregation: min agg size"));

      return factory;
    }

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
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      // Is it true? TEUCHOS_TEST_FOR_EXCEPTION(!paramList.isParameter("type"), Exceptions::RuntimeError, "TrilinosSmoother: parameter 'type' is mandatory");
      // type="" is default in TrilinosSmoother, but what happen then?

      std::string type="";            if(paramList.isParameter("type"))          type    = paramList.get<std::string>("type");
      int         overlap=0;          if(paramList.isParameter("overlap"))       overlap = paramList.get<int>        ("overlap");
      // std::string verbose;         if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params;  if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      // parameters from SmootherFactory
      //bool bKeepSmootherData = false; if(paramList.isParameter("keep smoother data")) bKeepSmootherData = paramList.get<bool>("keep smoother data");

      // Read in factory information for smoothers (if available...)
      // NOTE: only a selected number of factories can be used with the Trilinos smoother
      //       smoothers usually work with the global data available (which is A and the transfers P and R)

      Teuchos::RCP<TrilinosSmoother> trilSmoo = Teuchos::rcp(new TrilinosSmoother(type, params, overlap));

      if (paramList.isParameter("LineDetection_Layers")) {
        RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry("LineDetection_Layers"), factoryMapIn, factoryManagersIn);
        trilSmoo->SetFactory("LineDetection_Layers", generatingFact);
      }
      if (paramList.isParameter("LineDetection_VertLineIds")) {
        RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry("LineDetection_Layers"), factoryMapIn, factoryManagersIn);
        trilSmoo->SetFactory("LineDetection_Layers", generatingFact);
      }
      if (paramList.isParameter("CoarseNumZLayers")) {
        RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry("CoarseNumZLayers"), factoryMapIn, factoryManagersIn);
        trilSmoo->SetFactory("CoarseNumZLayers", generatingFact);
      }

      RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(Teuchos::null));
      Teuchos::ParameterList smooFactParams;
      //smooFactParams.setEntry("keep smoother data", paramList.getEntry("keep smoother data"));
      smooFact->SetParameterList(smooFactParams);
      smooFact->SetSmootherPrototypes(trilSmoo);
      return smooFact;
    }

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
    RCP<FactoryBase> BuildMatlabSmoother(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new MatlabSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "MatlabSmoother", Exceptions::RuntimeError, "");

      // Read in factory information for smoothers (if available...)
      // NOTE: only a selected number of factories can be used with the Trilinos smoother
      //       smoothers usually work with the global data available (which is A and the transfers P and R)

      Teuchos::RCP<MatlabSmoother> matSmoo = Teuchos::rcp(new MatlabSmoother(paramList));

      return rcp(new SmootherFactory(matSmoo));
    }
#endif

    RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList& paramList, const FactoryMap& /* factoryMapIn */, const FactoryManagerMap& /* factoryManagersIn */) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new DirectSolver()), Teuchos::null));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params)), Teuchos::null));
    }

    template <class T> // T must implement the Factory interface
    RCP<FactoryBase> BuildBlockedSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));

      // internal vector of factory managers
      std::vector<RCP<FactoryManager> > facManagers;

      // loop over all "block%i" sublists in parameter list
      int blockid = 1;
      bool blockExists = true;
      while (blockExists == true) {
        std::stringstream ss;
        ss << "block" << blockid;

        if(paramList.isSublist(ss.str()) == true) {
          // we either have a parameter group or we have a list of factories in here
          RCP<const ParameterList> b = rcp(new ParameterList(*sublist(paramListNonConst, ss.str())));

          RCP<FactoryManager> M = Teuchos::null;

          if (b->isParameter("group")) {
            // use a factory manager
            std::string facManagerName = b->get< std::string >("group");
            TEUCHOS_TEST_FOR_EXCEPTION(factoryManagersIn.count(facManagerName) != 1, Exceptions::RuntimeError, "Factory manager has not been found. Please check the spelling of the factory managers in your xml file.");
            RCP<FactoryManagerBase> Mb = factoryManagersIn.find(facManagerName)->second;
            M = Teuchos::rcp_dynamic_cast<FactoryManager>(Mb);
            TEUCHOS_TEST_FOR_EXCEPTION(M==Teuchos::null, Exceptions::RuntimeError, "Failed to cast FactoryManagerBase object to FactoryManager.");
          } else {
            // read in the list of factories
            M = rcp(new FactoryManager());
            for (ParameterList::ConstIterator param = b->begin(); param != b->end(); ++param) {
              RCP<const FactoryBase> p = BuildFactory(b->entry(param), factoryMapIn, factoryManagersIn);
              M->SetFactory(b->name(param),p);
            }
          }

          // add factory manager to internal vector of factory managers
          M->SetIgnoreUserData(true);
          facManagers.push_back(M);
          paramListNonConst->remove(ss.str());
          blockid++;
        } else {
          blockExists = false;
          break;
        }

      }

      // create a new blocked smoother
      RCP<T> bs = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // important: set block factory for A here!
      // TAW: 7/6/2016: We should not need to set/hardcode the blocked operator here.
      //                The user might want to overwrite this in the xml file, so just
      //                use what is declared as "A"
      //bs->SetFactory("A", MueLu::NoFactory::getRCP());

      for (int i = 0; i<Teuchos::as<int>(facManagers.size()); i++) {
        bs->AddFactoryManager(facManagers[i],i);
      }

      return rcp(new SmootherFactory(bs));
    }

#ifdef HAVE_MUELU_TEKO
    RCP<FactoryBase> BuildTekoSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));
      RCP<ParameterList> tekoParams = rcp(new ParameterList(paramListNonConst->sublist("Inverse Factory Library")));
      paramListNonConst->remove("Inverse Factory Library");

      // create a new blocked smoother
      RCP<TekoSmoother> bs = Build2<TekoSmoother>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // important: set block factory for A here!
      // TAW: 7/6/2016: We should not need to set/hardcode the blocked operator here.
      //                The user might want to overwrite this in the xml file, so just
      //                use what is declared as "A"
      //bs->SetFactory("A", MueLu::NoFactory::getRCP());

      // Set Teko parameters ("Inverse Factory Library")
      bs->SetTekoParameters(tekoParams);

      return rcp(new SmootherFactory(bs));
    }
#endif

    RCP<FactoryBase> BuildBlockedDirectSolver(const Teuchos::ParameterList& /* paramList */, const FactoryMap& /* factoryMapIn */, const FactoryManagerMap& /* factoryManagersIn */) const {
      //if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new BlockedDirectSolver())));

      /*TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params))));*/
    }

    //RCP<FactoryBase> BuildBlockedPFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
    //  RCP<BlockedPFactory> pfac = rcp(new BlockedPFactory());

    template <class T> // T must implement the Factory interface
    RCP<T> BuildBlockedFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> pfac = Teuchos::null;

      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));

      // internal vector of factory managers
      std::vector<RCP<FactoryManager> > facManagers;

      // loop over all "block%i" sublists in parameter list
      int blockid = 1;
      bool blockExists = true;
      while (blockExists == true) {
        std::stringstream ss;
        ss << "block" << blockid;

        if(paramList.isSublist(ss.str()) == true) {
          // we either have a parameter group or we have a list of factories in here
          RCP<const ParameterList> b = rcp(new ParameterList(*sublist(paramListNonConst, ss.str())));

          RCP<FactoryManager> M = Teuchos::null;

          if (b->isParameter("group")) {
            // use a factory manager
            std::string facManagerName = b->get< std::string >("group");
            TEUCHOS_TEST_FOR_EXCEPTION(factoryManagersIn.count(facManagerName) != 1, Exceptions::RuntimeError, "Factory manager has not been found. Please check the spelling of the factory managers in your xml file.");
            RCP<FactoryManagerBase> Mb = factoryManagersIn.find(facManagerName)->second;
            M = Teuchos::rcp_dynamic_cast<FactoryManager>(Mb);
            TEUCHOS_TEST_FOR_EXCEPTION(M==Teuchos::null, Exceptions::RuntimeError, "Failed to cast FactoryManagerBase object to FactoryManager.");
          } else {
            // read in the list of factories
            M = rcp(new FactoryManager());
            for (ParameterList::ConstIterator param = b->begin(); param != b->end(); ++param) {
              RCP<const FactoryBase> p = BuildFactory(b->entry(param), factoryMapIn, factoryManagersIn);
              M->SetFactory(b->name(param),p);
            }
          }

          // add factory manager to internal vector of factory managers
          M->SetIgnoreUserData(true);
          facManagers.push_back(M);
          paramListNonConst->remove(ss.str());
          blockid++;
        } else {
          blockExists = false;
          break;
        }

      }

      // build BlockedPFactory (without sub block information)
      pfac = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // add FactoryManager objects
      for(size_t i = 0; i<facManagers.size(); i++) {
        pfac->AddFactoryManager(facManagers[i]); // add factory manager
      }

      return pfac;
    }


    template <class T> // T must implement the Factory interface
    RCP<T> BuildBlockedCoordFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> pfac = Teuchos::null;

      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));

      // internal vector of factory managers
      std::vector<RCP<const FactoryBase> > facBase;

      // loop over all "block%i" sublists in parameter list
      int blockid = 1;
      bool blockExists = true;
      while (blockExists == true) {
        std::stringstream ss;
        ss << "block" << blockid;

        if(paramList.isSublist(ss.str()) == true) {
          // we either have a parameter group or we have a list of factories in here
          RCP<const ParameterList> b = rcp(new ParameterList(*sublist(paramListNonConst, ss.str())));

            // read in the list of factories
            for (ParameterList::ConstIterator param = b->begin(); param != b->end(); ++param) {
              RCP<const FactoryBase> p = BuildFactory(b->entry(param), factoryMapIn, factoryManagersIn);
              facBase.push_back(p);
            }

          // add factory manager to internal vector of factory managers
          paramListNonConst->remove(ss.str());
          blockid++;
        } else {
          blockExists = false;
          break;
        }

      }

      // build BlockedPFactory (without sub block information)
      pfac = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // add FactoryManager objects
      for(size_t i = 0; i<facBase.size(); i++) {
        pfac->AddFactory(facBase[i]); // add factory manager
      }

      return pfac;
    }

  }; // class
} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

  // TODO: handle factory parameters
  // TODO: parameter validator
  // TODO: static
  // TODO: default parameters should not be duplicated here and on the Factory (ex: default for overlap (=0) is defined both here and on TrilinosSmoother constructors)
