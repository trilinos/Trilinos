
INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

SET(HAVE_STOKHOS_DAKOTA ${${PACKAGE_NAME}_ENABLE_TriKota})

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add path to Dakota headers in installed directory
IF (${PACKAGE_NAME}_ENABLE_TriKota)
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_BINARY_DIR}/../../TriKota/Dakota/install/include
  )
ENDIF()


SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Stokhos.hpp
  Stokhos_ConfigDefs.h
  Stokhos_Version.hpp
  Stokhos_gaussq.h
  Stokhos_DynamicArrayTraits.hpp
  Stokhos_OneDOrthogPolyBasis.hpp
  Stokhos_RecurrenceBasis.hpp
  Stokhos_RecurrenceBasisImp.hpp
  Stokhos_HermiteBasis.hpp
  Stokhos_HermiteBasisImp.hpp
  Stokhos_JacobiBasis.hpp
  Stokhos_JacobiBasisImp.hpp
  Stokhos_DiscretizedStieltjesBasis.hpp
  Stokhos_DiscretizedStieltjesBasisImp.hpp
  Stokhos_RysBasis.hpp
  Stokhos_LegendreBasis.hpp
  Stokhos_LegendreBasisImp.hpp
  Stokhos_ClenshawCurtisLegendreBasis.hpp
  Stokhos_ClenshawCurtisLegendreBasisImp.hpp
  Stokhos_StieltjesPCEBasis.hpp
  Stokhos_StieltjesPCEBasisImp.hpp
  Stokhos_StieltjesBasis.hpp
  Stokhos_StieltjesBasisImp.hpp
  Stokhos_Lanczos.hpp
  Stokhos_LanczosPCEBasis.hpp
  Stokhos_LanczosPCEBasisImp.hpp
  Stokhos_LanczosProjPCEBasis.hpp
  Stokhos_LanczosProjPCEBasisImp.hpp
  Stokhos_PecosOneDOrthogPolyBasis.hpp
  Stokhos_PecosOneDOrthogPolyBasisImp.hpp
  Stokhos_OrthogPolyBasis.hpp
  Stokhos_ProductBasis.hpp
  Stokhos_DerivBasis.hpp
  Stokhos_ProductBasisUtils.hpp
  Stokhos_CompletePolynomialBasis.hpp
  Stokhos_CompletePolynomialBasisImp.hpp
  Stokhos_GramSchmidtBasis.hpp
  Stokhos_GramSchmidtBasisImp.hpp
  Stokhos_ReducedPCEBasis.hpp
  Stokhos_SDMUtils.hpp
  Stokhos_OrthogonalizationFactory.hpp
  Stokhos_GSReducedPCEBasisBase.hpp
  Stokhos_GSReducedPCEBasisBaseImp.hpp
  Stokhos_MonomialProjGramSchmidtPCEBasis.hpp
  Stokhos_MonomialProjGramSchmidtPCEBasisImp.hpp
  Stokhos_MonomialProjGramSchmidtPCEBasis2.hpp
  Stokhos_MonomialProjGramSchmidtPCEBasis2Imp.hpp
  Stokhos_MonomialGramSchmidtPCEBasis.hpp
  Stokhos_MonomialGramSchmidtPCEBasisImp.hpp
  Stokhos_ProductLanczosPCEBasis.hpp
  Stokhos_ProductLanczosPCEBasisImp.hpp
  Stokhos_ProductLanczosGramSchmidtPCEBasis.hpp
  Stokhos_ProductLanczosGramSchmidtPCEBasisImp.hpp
  Stokhos_ReducedBasisFactory.hpp
  Stokhos_ReducedBasisFactoryImp.hpp
  Stokhos_ReducedQuadratureFactory.hpp
  Stokhos_ReducedQuadratureFactoryImp.hpp
  Stokhos_BasisFactory.hpp
  Stokhos_BasisFactoryImp.hpp
  Stokhos_SparseArray.hpp
  Stokhos_Sparse3Tensor.hpp
  Stokhos_Sparse3TensorImp.hpp
  Stokhos_Dense3Tensor.hpp
  Stokhos_Dense3TensorImp.hpp
  Stokhos_Sparse3TensorUtilities.hpp
  Stokhos_Quadrature.hpp
  Stokhos_TensorProductQuadrature.hpp
  Stokhos_TensorProductQuadratureImp.hpp
  Stokhos_SparseGridQuadrature.hpp
  Stokhos_SparseGridQuadratureImp.hpp
  Stokhos_AnisoSparseGridQuadrature.hpp
  Stokhos_AnisoSparseGridQuadratureImp.hpp
  Stokhos_UserDefinedQuadrature.hpp
  Stokhos_UserDefinedQuadratureImp.hpp
  Stokhos_QuadratureFactory.hpp
  Stokhos_QuadratureFactoryImp.hpp
  Stokhos_StandardStorage.hpp
  Stokhos_StaticStandardStorage.hpp
  Stokhos_StaticFixedStandardStorage.hpp
  Stokhos_OrthogPolyApprox.hpp
  Stokhos_OrthogPolyApproxImp.hpp
  Stokhos_OrthogPolyExpansion.hpp
  Stokhos_DivisionExpansionStrategy.hpp
  Stokhos_MeanBasedDivisionExpansionStrategy.hpp
  Stokhos_DenseDirectDivisionExpansionStrategy.hpp
  Stokhos_SPDDenseDirectDivisionExpansionStrategy.hpp
  Stokhos_GMRESDivisionExpansionStrategy.hpp
  Stokhos_CGDivisionExpansionStrategy.hpp
  Stokhos_OrthogPolyExpansionBase.hpp
  Stokhos_OrthogPolyExpansionBaseImp.hpp
  Stokhos_ConstantOrthogPolyExpansion.hpp
  Stokhos_ConstantOrthogPolyExpansionImp.hpp
  Stokhos_AlgebraicOrthogPolyExpansion.hpp
  Stokhos_AlgebraicOrthogPolyExpansionImp.hpp
  Stokhos_QuadOrthogPolyExpansion.hpp
  Stokhos_QuadOrthogPolyExpansionImp.hpp
  Stokhos_DerivOrthogPolyExpansion.hpp
  Stokhos_DerivOrthogPolyExpansionImp.hpp
  Stokhos_ForUQTKOrthogPolyExpansion.hpp
  Stokhos_ForUQTKOrthogPolyExpansionImp.hpp
  Stokhos_ExpansionFactory.hpp
  Stokhos_ExpansionFactoryImp.hpp
  Stokhos_StieltjesGramSchmidtBuilder.hpp
  Stokhos_StieltjesGramSchmidtBuilderImp.hpp
  Stokhos_EpetraSparse3Tensor.hpp
  Stokhos_ParallelData.hpp
  Stokhos_SGOperator.hpp
  Stokhos_MatrixFreeOperator.hpp
  Stokhos_KLMatrixFreeOperator.hpp
  Stokhos_KLReducedMatrixFreeOperator.hpp
  Stokhos_FullyAssembledOperator.hpp
  Stokhos_SGOperatorFactory.hpp
  Stokhos_SGPreconditioner.hpp
  Stokhos_MeanBasedPreconditioner.hpp
  Stokhos_ApproxGaussSeidelPreconditioner.hpp
  Stokhos_ApproxJacobiPreconditioner.hpp
  Stokhos_ApproxSchurComplementPreconditioner.hpp
  Stokhos_KroneckerProductPreconditioner.hpp
  Stokhos_FullyAssembledPreconditioner.hpp
  Stokhos_SGPreconditionerFactory.hpp
  Stokhos_AbstractPreconditionerFactory.hpp
  Stokhos_IfpackPreconditionerFactory.hpp
  Stokhos_MLPreconditionerFactory.hpp
  Stokhos_PreconditionerFactory.hpp
  Stokhos_MLPrecOp.hpp
  Stokhos_SGModelEvaluator.hpp
  Stokhos_SGQuadModelEvaluator.hpp
  Stokhos_SGInverseModelEvaluator.hpp
  Stokhos_ResponseStatisticModelEvaluator.hpp
  Stokhos_ProductContainer.hpp
  Stokhos_ProductContainerImp.hpp
  Stokhos_VectorOrthogPoly.hpp
  Stokhos_VectorOrthogPolyImp.hpp
  Stokhos_Epetra.hpp
  Stokhos_VectorOrthogPolyTraitsEpetra.hpp
  Stokhos_ProductEpetraVector.hpp
  Stokhos_ProductEpetraMultiVector.hpp
  Stokhos_ProductEpetraOperator.hpp
  Stokhos_ProductEpetraMultiVectorOperator.hpp
  Stokhos_EpetraVectorOrthogPoly.hpp
  Stokhos_EpetraMultiVectorOrthogPoly.hpp
  Stokhos_EpetraOperatorOrthogPoly.hpp
  Stokhos_EpetraMultiVectorOperatorOrthogPoly.hpp
  Stokhos_PCECovarianceOp.hpp
  Stokhos_PCEAnasaziKL.hpp
  Stokhos_KL_OneDEigenPair.hpp
  Stokhos_KL_ProductEigenPair.hpp
  Stokhos_KL_OneDExponentialCovarianceFunction.hpp
  Stokhos_KL_OneDExponentialCovarianceFunctionImp.hpp
  Stokhos_KL_ExponentialRandomField.hpp
  Stokhos_KL_ExponentialRandomFieldImp.hpp
  Stokhos_EpetraMultiVectorOperator.hpp
  Stokhos_BlockDiagonalOperator.hpp
  Stokhos_MPPreconditioner.hpp
  Stokhos_MPBlockDiagonalPreconditioner.hpp
  Stokhos_MPMeanBasedPreconditioner.hpp
  Stokhos_MPPreconditionerFactory.hpp
  Stokhos_MPModelEvaluator.hpp
  Stokhos_SGQuadMPModelEvaluator.hpp
  Stokhos_MPModelEvaluatorAdapter.hpp
  Stokhos_MPInverseModelEvaluator.hpp
  Stokhos_InterlacedOperator.hpp
  Stokhos_SGModelEvaluator_Interlaced.hpp
  Stokhos_SGModelEvaluator_Adaptive.hpp
  Stokhos_BasisInteractionGraph.hpp
  Stokhos_AdaptivityUtils.hpp
  Stokhos_AdaptivityManager.hpp
  Stokhos_Operator.hpp
  Stokhos_DiagPreconditioner.hpp
  Stokhos_JacobiPreconditioner.hpp
  Stokhos_GSPreconditioner.hpp
  Stokhos_SchurPreconditioner.hpp
  Stokhos_SchurPreconditionerImp.hpp
  Stokhos_InversePreconditioner.hpp
  Stokhos_BlockPreconditioner.hpp
  Stokhos_BlockPreconditionerImp.hpp
  )

SET(SOURCES ${SOURCES}
  Stokhos_Version.cpp
  Stokhos_SparseGridQuadrature.cpp
  Stokhos_AnisoSparseGridQuadrature.cpp
  Stokhos_EpetraSparse3Tensor.cpp
  Stokhos_ParallelData.cpp
  Stokhos_MatrixFreeOperator.cpp
  Stokhos_KLMatrixFreeOperator.cpp
  Stokhos_KLReducedMatrixFreeOperator.cpp
  Stokhos_FullyAssembledOperator.cpp
  Stokhos_SGOperatorFactory.cpp
  Stokhos_MeanBasedPreconditioner.cpp
  Stokhos_ApproxGaussSeidelPreconditioner.cpp
  Stokhos_ApproxJacobiPreconditioner.cpp
  Stokhos_ApproxSchurComplementPreconditioner.cpp
  Stokhos_KroneckerProductPreconditioner.cpp
  Stokhos_FullyAssembledPreconditioner.cpp
  Stokhos_SGPreconditionerFactory.cpp
  Stokhos_IfpackPreconditionerFactory.cpp
  Stokhos_MLPreconditionerFactory.cpp
  Stokhos_PreconditionerFactory.cpp
  Stokhos_MLPrecOp.cpp
  Stokhos_ProductEpetraVector.cpp
  Stokhos_ProductEpetraMultiVector.cpp
  Stokhos_ProductEpetraOperator.cpp
  Stokhos_ProductEpetraMultiVectorOperator.cpp
  Stokhos_EpetraVectorOrthogPoly.cpp
  Stokhos_EpetraMultiVectorOrthogPoly.cpp
  Stokhos_EpetraOperatorOrthogPoly.cpp
  Stokhos_EpetraMultiVectorOperatorOrthogPoly.cpp
  Stokhos_SGModelEvaluator.cpp
  Stokhos_SGQuadModelEvaluator.cpp
  Stokhos_SGInverseModelEvaluator.cpp
  Stokhos_ResponseStatisticModelEvaluator.cpp
  Stokhos_PCECovarianceOp.cpp
  Stokhos_PCEAnasaziKL.cpp
  Stokhos_EpetraMultiVectorOperator.cpp
  Stokhos_BlockDiagonalOperator.cpp
  Stokhos_MPBlockDiagonalPreconditioner.cpp
  Stokhos_MPMeanBasedPreconditioner.cpp
  Stokhos_MPPreconditionerFactory.cpp
  Stokhos_MPModelEvaluator.cpp
  Stokhos_SGQuadMPModelEvaluator.cpp
  Stokhos_MPModelEvaluatorAdapter.cpp
  Stokhos_MPInverseModelEvaluator.cpp
  Stokhos_InterlacedOperator.cpp
  Stokhos_SGModelEvaluator_Interlaced.cpp
  Stokhos_SGModelEvaluator_Adaptive.cpp
  Stokhos_BasisInteractionGraph.cpp
  Stokhos_AdaptivityUtils.cpp
  Stokhos_AdaptivityManager.cpp
)

# Sacado
IF (Stokhos_ENABLE_Sacado)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado)
  SET(HEADERS ${HEADERS}
    sacado/Stokhos_Sacado.hpp
    sacado/Stokhos_Sacado_MathFunctions.hpp
    sacado/Sacado_PCE_OrthogPoly.hpp
    sacado/Sacado_PCE_OrthogPolyImp.hpp
    sacado/Sacado_PCE_OrthogPolyTraits.hpp
    sacado/Sacado_PCE_ScalarTraitsImp.hpp
    sacado/Sacado_ETPCE_OrthogPoly.hpp
    sacado/Sacado_ETPCE_OrthogPolyImp.hpp
    sacado/Sacado_ETPCE_OrthogPolyTraits.hpp
    sacado/Sacado_ETPCE_OrthogPolyOps.hpp
    sacado/Sacado_ETV_Vector.hpp
    sacado/Sacado_ETV_VectorImp.hpp
    sacado/Sacado_ETV_VectorOps.hpp
    sacado/Sacado_ETV_ScalarTraitsImp.hpp
    sacado/Sacado_ETV_VectorTraits.hpp
    sacado/Sacado_ETV_Vector2.hpp
    sacado/Sacado_ETV_Vector2Imp.hpp
    sacado/Sacado_ETV_Vector2Traits.hpp
    )
  SET(SOURCES ${SOURCES}
    sacado/Sacado_PCE_OrthogPoly.cpp
    sacado/Sacado_ETPCE_OrthogPoly.cpp
    )
  IF (Stokhos_ENABLE_KokkosArray)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos)
    SET(HEADERS ${HEADERS}
      sacado/kokkos/Stokhos_Sacado_Kokkos.hpp
      sacado/kokkos/Stokhos_Sacado_Kokkos_MathFunctions.hpp
      sacado/kokkos/Stokhos_StaticArrayTraits.hpp
      sacado/kokkos/Stokhos_DynArrayTraits.hpp
      sacado/kokkos/Stokhos_LocalStorage.hpp
      sacado/kokkos/Stokhos_StaticStorage.hpp
      sacado/kokkos/Stokhos_StaticFixedStorage.hpp
      sacado/kokkos/Stokhos_DynamicStorage.hpp
      sacado/kokkos/Stokhos_DynamicStridedStorage.hpp
      sacado/kokkos/Stokhos_DynamicThreadedStorage.hpp
      sacado/kokkos/Stokhos_DynamicThreadedStorage_cuda.hpp
      sacado/kokkos/Sacado_MP_Vector.hpp
      sacado/kokkos/Sacado_MP_VectorTraits.hpp
      sacado/kokkos/Sacado_MP_Vector_ops.hpp
    )
  ENDIF()
  IF (Stokhos_ENABLE_MueLu)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/muelu)
    TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/muelu)
    SET(HEADERS ${HEADERS}
      muelu/Stokhos_MueLu.hpp
      muelu/Stokhos_MueLu_QR_Interface_decl.hpp
      muelu/Stokhos_MueLu_QR_Interface_def.hpp
      )
    IF(MueLu_ENABLE_EXPLICIT_INSTANTIATION)
      APPEND_GLOB(SOURCES muelu/explicit_instantiation/*.cpp)
    ENDIF()
  ENDIF()
  IF (Stokhos_ENABLE_Tpetra)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpetra)
    IF(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
      SET(SOURCES ${SOURCES}
	tpetra/tpetra_pce_instantiations.cpp
	)
      SET(HEADERS ${HEADERS}
	tpetra/kokkos_pce_specializations.hpp
	)
    ENDIF()
  ENDIF()
  IF (Stokhos_ENABLE_Ifpack2)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ifpack2)
    IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
      SET(SOURCES ${SOURCES}
	ifpack2/ifpack2_pce_instantiations.cpp
	)
    ENDIF()
  ENDIF()
ENDIF()

# NOX
IF (Stokhos_ENABLE_NOX AND NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND
    NOX_ENABLE_EpetraExt)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/nox)
  SET(HEADERS ${HEADERS}
    Stokhos_GaussSeidelPreconditioner.hpp
    nox/NOX_Epetra_LinearSystem_SGJacobi.hpp
    nox/NOX_Epetra_LinearSystem_SGGS.hpp
    nox/NOX_Epetra_LinearSystem_MPBD.hpp
    )
  SET(SOURCES ${SOURCES}
    Stokhos_GaussSeidelPreconditioner.cpp
    nox/NOX_Epetra_LinearSystem_SGJacobi.cpp
    nox/NOX_Epetra_LinearSystem_SGGS.cpp
    nox/NOX_Epetra_LinearSystem_MPBD.cpp
    )
ENDIF()

# CUDA/Thrust
IF (Stokhos_ENABLE_Thrust)
  # this is a manual specification of include dirs, necessary because 
  # CUDA_COMPILE doesn't recognize the package architecture system.
  GET_PROPERTY(Stokhos_CUDA_INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
  LIST(APPEND Stokhos_CUDA_INC_DIRS ${Teuchos_INCLUDE_DIRS})
  IF (${PACKAGE_NAME}_ENABLE_TriKota)
    LIST(APPEND Stokhos_CUDA_INC_DIRS ${TriKota_INCLUDE_DIRS})
  ENDIF()
  #message("Stokhos_CUDA_INC_DIRS = ${Stokhos_CUDA_INC_DIRS}")
  CUDA_INCLUDE_DIRECTORIES(${Stokhos_CUDA_INC_DIRS})
  CUDA_COMPILE(CUDA_SOURCES
    Stokhos_CUDAStorage.cu 
    Stokhos_CUDAQuadOrthogPolyExpansion.cu
  )
  SET(HEADERS ${HEADERS}
    Stokhos_CUDAStorage.hpp
    Stokhos_CUDAStorageImp.hpp
    Stokhos_CUDAQuadOrthogPolyExpansion.hpp
    Stokhos_CUDAQuadOrthogPolyExpansionImp.hpp
    )
  SET(SOURCES ${SOURCES} ${CUDA_SOURCES})
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stokhos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
 