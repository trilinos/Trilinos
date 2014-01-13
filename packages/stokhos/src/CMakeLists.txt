INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

SET(HAVE_STOKHOS_DAKOTA ${${PACKAGE_NAME}_ENABLE_TriKota})

# OpenMP support -- Add our own option to turn it off
# without turning it off globally
SET(Stokhos_ENABLE_OpenMP_Default OFF)
IF(Trilinos_ENABLE_OpenMP)
  SET(Stokhos_ENABLE_OpenMP_Default ON)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Stokhos_ENABLE_OpenMP
  HAVE_STOKHOS_OPENMP
  "Enable OpenMP support in Stokhos."
  ${Stokhos_ENABLE_OpenMP_Default}
  )

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
  Stokhos_GaussPattersonLegendreBasis.hpp
  Stokhos_GaussPattersonLegendreBasisImp.hpp
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
  Stokhos_GrowthRules.hpp
  Stokhos_ProductBasisUtils.hpp
  Stokhos_CompletePolynomialBasis.hpp
  Stokhos_CompletePolynomialBasisImp.hpp
  Stokhos_TensorProductBasis.hpp
  Stokhos_TensorProductBasisImp.hpp
  Stokhos_TotalOrderBasis.hpp
  Stokhos_TotalOrderBasisImp.hpp
  Stokhos_SmolyakBasis.hpp
  Stokhos_SmolyakBasisImp.hpp
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
  Stokhos_Sparse3TensorPartition.hpp
  Stokhos_LTBSparse3Tensor.hpp
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
  Stokhos_PseudoSpectralOperator.hpp
  Stokhos_TensorProductPseudoSpectralOperator.hpp
  Stokhos_SmolyakPseudoSpectralOperator.hpp
  Stokhos_SmolyakPseudoSpectralOperatorImp.hpp
  Stokhos_QuadraturePseudoSpectralOperator.hpp
  Stokhos_PseudoSpectralOperatorFactory.hpp
  Stokhos_PseudoSpectralOperatorFactoryImp.hpp
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
  Stokhos_PseudoSpectralOrthogPolyExpansion.hpp
  Stokhos_PseudoSpectralOrthogPolyExpansionImp.hpp
  Stokhos_DerivOrthogPolyExpansion.hpp
  Stokhos_DerivOrthogPolyExpansionImp.hpp
  Stokhos_ForUQTKOrthogPolyExpansion.hpp
  Stokhos_ForUQTKOrthogPolyExpansionImp.hpp
  Stokhos_ExpansionFactory.hpp
  Stokhos_ExpansionFactoryImp.hpp
  Stokhos_StieltjesGramSchmidtBuilder.hpp
  Stokhos_StieltjesGramSchmidtBuilderImp.hpp
  Stokhos_KL_OneDEigenPair.hpp
  Stokhos_KL_ProductEigenPair.hpp
  Stokhos_KL_OneDExponentialCovarianceFunction.hpp
  Stokhos_KL_OneDExponentialCovarianceFunctionImp.hpp
  Stokhos_KL_ExponentialRandomField.hpp
  Stokhos_KL_ExponentialRandomFieldImp.hpp
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
)

# Epetra
IF (Stokhos_ENABLE_EpetraExt)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)
  SET(HEADERS ${HEADERS}
    epetra/Stokhos_Sparse3TensorUtilities.hpp
    epetra/Stokhos_EpetraSparse3Tensor.hpp
    epetra/Stokhos_SGOperator.hpp
    epetra/Stokhos_MatrixFreeOperator.hpp
    epetra/Stokhos_KLMatrixFreeOperator.hpp
    epetra/Stokhos_KLReducedMatrixFreeOperator.hpp
    epetra/Stokhos_FullyAssembledOperator.hpp
    epetra/Stokhos_SGOperatorFactory.hpp
    epetra/Stokhos_SGPreconditioner.hpp
    epetra/Stokhos_MeanBasedPreconditioner.hpp
    epetra/Stokhos_ApproxGaussSeidelPreconditioner.hpp
    epetra/Stokhos_ApproxJacobiPreconditioner.hpp
    epetra/Stokhos_ApproxSchurComplementPreconditioner.hpp
    epetra/Stokhos_KroneckerProductPreconditioner.hpp
    epetra/Stokhos_FullyAssembledPreconditioner.hpp
    epetra/Stokhos_SGPreconditionerFactory.hpp
    epetra/Stokhos_AbstractPreconditionerFactory.hpp
    epetra/Stokhos_IfpackPreconditionerFactory.hpp
    epetra/Stokhos_MLPreconditionerFactory.hpp
    epetra/Stokhos_PreconditionerFactory.hpp
    epetra/Stokhos_MLPrecOp.hpp
    epetra/Stokhos_SGModelEvaluator.hpp
    epetra/Stokhos_SGQuadModelEvaluator.hpp
    epetra/Stokhos_SGInverseModelEvaluator.hpp
    epetra/Stokhos_ResponseStatisticModelEvaluator.hpp
    epetra/Stokhos_Epetra.hpp
    epetra/Stokhos_ParallelData.hpp
    epetra/Stokhos_PCECovarianceOp.hpp
    epetra/Stokhos_PCEAnasaziKL.hpp
    epetra/Stokhos_ProductContainer.hpp
    epetra/Stokhos_ProductContainerImp.hpp
    epetra/Stokhos_VectorOrthogPoly.hpp
    epetra/Stokhos_VectorOrthogPolyImp.hpp
    epetra/Stokhos_VectorOrthogPolyTraitsEpetra.hpp
    epetra/Stokhos_ProductEpetraVector.hpp
    epetra/Stokhos_ProductEpetraMultiVector.hpp
    epetra/Stokhos_ProductEpetraOperator.hpp
    epetra/Stokhos_ProductEpetraMultiVectorOperator.hpp
    epetra/Stokhos_EpetraVectorOrthogPoly.hpp
    epetra/Stokhos_EpetraMultiVectorOrthogPoly.hpp
    epetra/Stokhos_EpetraOperatorOrthogPoly.hpp
    epetra/Stokhos_EpetraMultiVectorOperatorOrthogPoly.hpp
    epetra/Stokhos_EpetraMultiVectorOperator.hpp
    epetra/Stokhos_BlockDiagonalOperator.hpp
    epetra/Stokhos_MPPreconditioner.hpp
    epetra/Stokhos_MPBlockDiagonalPreconditioner.hpp
    epetra/Stokhos_MPMeanBasedPreconditioner.hpp
    epetra/Stokhos_MPPreconditionerFactory.hpp
    epetra/Stokhos_MPModelEvaluator.hpp
    epetra/Stokhos_SGQuadMPModelEvaluator.hpp
    epetra/Stokhos_MPModelEvaluatorAdapter.hpp
    epetra/Stokhos_MPInverseModelEvaluator.hpp
    epetra/Stokhos_InterlacedOperator.hpp
    epetra/Stokhos_SGModelEvaluator_Interlaced.hpp
    epetra/Stokhos_SGModelEvaluator_Adaptive.hpp
    epetra/Stokhos_BasisInteractionGraph.hpp
    epetra/Stokhos_AdaptivityUtils.hpp
    epetra/Stokhos_AdaptivityManager.hpp
    )

  SET(SOURCES ${SOURCES}
    epetra/Stokhos_EpetraSparse3Tensor.cpp
    epetra/Stokhos_MatrixFreeOperator.cpp
    epetra/Stokhos_KLMatrixFreeOperator.cpp
    epetra/Stokhos_KLReducedMatrixFreeOperator.cpp
    epetra/Stokhos_FullyAssembledOperator.cpp
    epetra/Stokhos_SGOperatorFactory.cpp
    epetra/Stokhos_MeanBasedPreconditioner.cpp
    epetra/Stokhos_ApproxGaussSeidelPreconditioner.cpp
    epetra/Stokhos_ApproxJacobiPreconditioner.cpp
    epetra/Stokhos_ApproxSchurComplementPreconditioner.cpp
    epetra/Stokhos_KroneckerProductPreconditioner.cpp
    epetra/Stokhos_FullyAssembledPreconditioner.cpp
    epetra/Stokhos_SGPreconditionerFactory.cpp
    epetra/Stokhos_IfpackPreconditionerFactory.cpp
    epetra/Stokhos_MLPreconditionerFactory.cpp
    epetra/Stokhos_PreconditionerFactory.cpp
    epetra/Stokhos_MLPrecOp.cpp
    epetra/Stokhos_ParallelData.cpp
    epetra/Stokhos_PCECovarianceOp.cpp
    epetra/Stokhos_PCEAnasaziKL.cpp
    epetra/Stokhos_ProductEpetraVector.cpp
    epetra/Stokhos_ProductEpetraMultiVector.cpp
    epetra/Stokhos_ProductEpetraOperator.cpp
    epetra/Stokhos_ProductEpetraMultiVectorOperator.cpp
    epetra/Stokhos_EpetraVectorOrthogPoly.cpp
    epetra/Stokhos_EpetraMultiVectorOrthogPoly.cpp
    epetra/Stokhos_EpetraOperatorOrthogPoly.cpp
    epetra/Stokhos_EpetraMultiVectorOperatorOrthogPoly.cpp
    epetra/Stokhos_SGModelEvaluator.cpp
    epetra/Stokhos_SGQuadModelEvaluator.cpp
    epetra/Stokhos_SGInverseModelEvaluator.cpp
    epetra/Stokhos_ResponseStatisticModelEvaluator.cpp
    epetra/Stokhos_EpetraMultiVectorOperator.cpp
    epetra/Stokhos_BlockDiagonalOperator.cpp
    epetra/Stokhos_MPBlockDiagonalPreconditioner.cpp
    epetra/Stokhos_MPMeanBasedPreconditioner.cpp
    epetra/Stokhos_MPPreconditionerFactory.cpp
    epetra/Stokhos_MPModelEvaluator.cpp
    epetra/Stokhos_SGQuadMPModelEvaluator.cpp
    epetra/Stokhos_MPModelEvaluatorAdapter.cpp
    epetra/Stokhos_MPInverseModelEvaluator.cpp
    epetra/Stokhos_InterlacedOperator.cpp
    epetra/Stokhos_SGModelEvaluator_Interlaced.cpp
    epetra/Stokhos_SGModelEvaluator_Adaptive.cpp
    epetra/Stokhos_BasisInteractionGraph.cpp
    epetra/Stokhos_AdaptivityUtils.cpp
    epetra/Stokhos_AdaptivityManager.cpp
    )
ENDIF()

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
    sacado/Sacado_ETPCE_ExpressionTraits.hpp
    sacado/Sacado_ETV_Vector.hpp
    sacado/Sacado_ETV_VectorImp.hpp
    sacado/Sacado_ETV_VectorOps.hpp
    sacado/Sacado_ETV_ScalarTraitsImp.hpp
    sacado/Sacado_ETV_VectorTraits.hpp
    sacado/Sacado_ETV_ExpressionTraits.hpp
    sacado/Sacado_ETV_Vector2.hpp
    sacado/Sacado_ETV_Vector2Imp.hpp
    sacado/Sacado_ETV_Vector2Traits.hpp
    )
  SET(SOURCES ${SOURCES}
    sacado/Sacado_PCE_OrthogPoly.cpp
    sacado/Sacado_ETPCE_OrthogPoly.cpp
    )
  IF (Stokhos_ENABLE_KokkosCore)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector)
    SET(HEADERS ${HEADERS}
      sacado/kokkos/Stokhos_Sacado_Kokkos.hpp
      sacado/kokkos/Stokhos_Sacado_Kokkos_MathFunctions.hpp
      sacado/kokkos/Stokhos_KokkosTraits.hpp
      sacado/kokkos/Stokhos_StaticArrayTraits.hpp
      sacado/kokkos/Stokhos_DynArrayTraits.hpp
      sacado/kokkos/Stokhos_LocalStorage.hpp
      sacado/kokkos/Stokhos_StaticStorage.hpp
      sacado/kokkos/Stokhos_StaticFixedStorage.hpp
      sacado/kokkos/Stokhos_DynamicStorage.hpp
      sacado/kokkos/Stokhos_DynamicStridedStorage.hpp
      sacado/kokkos/Stokhos_DynamicThreadedStorage.hpp
      sacado/kokkos/Stokhos_DynamicThreadedStorage_cuda.hpp
      sacado/kokkos/Stokhos_ViewStorage.hpp
      sacado/kokkos/Stokhos_ViewStridedStorage.hpp
      sacado/kokkos/vector/Sacado_MP_Vector.hpp
      sacado/kokkos/vector/Sacado_MP_VectorTraits.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_ops.hpp
      sacado/kokkos/vector/Sacado_MP_ExpressionTraits.hpp
      sacado/kokkos/vector/Kokkos_View_MP_Vector.hpp
      sacado/kokkos/vector/Kokkos_Atomic_MP_Vector.hpp
      )
    IF (Stokhos_ENABLE_KokkosLinAlg)
      INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/linalg)
      SET(HEADERS ${HEADERS}
        sacado/kokkos/vector/linalg/Kokkos_CrsMatrix_MP_Vector.hpp
        )
      IF(Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_Cuda)
        SET(HEADERS ${HEADERS}
          sacado/kokkos/vector/linalg/Kokkos_CrsMatrix_MP_Vector_Cuda.hpp
          )
      ENDIF()
    ENDIF()
    IF (Stokhos_ENABLE_Tpetra AND
        NOT Tpetra_ENABLE_Kokkos_DistObject AND
        NOT Tpetra_ENABLE_Kokkos_Refactor)
      INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/tpetra)
      IF(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
        SET(HEADERS ${HEADERS}
          sacado/kokkos/vector/tpetra/Stokhos_Tpetra_ETI_Helpers_MP_Vector.hpp
          )
        SET(SOURCES ${SOURCES}
          sacado/kokkos/vector/tpetra/Tpetra_MultiVector.cpp
          sacado/kokkos/vector/tpetra/Tpetra_Vector.cpp
          sacado/kokkos/vector/tpetra/Tpetra_CrsMatrix.cpp
          sacado/kokkos/vector/tpetra/Tpetra_RowMatrix.cpp
          sacado/kokkos/vector/tpetra/Tpetra_DistObject.cpp
          )
      ENDIF()
    ENDIF()
    IF (Stokhos_ENABLE_Belos)
      INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/belos)
      SET(HEADERS ${HEADERS}
        sacado/kokkos/vector/belos/Belos_TpetraAdapter_MP_Vector.hpp
        )
    ENDIF()
  ENDIF()
  IF (Stokhos_ENABLE_Experimental AND Stokhos_ENABLE_MueLu)
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
  IF (Stokhos_ENABLE_Experimental AND Stokhos_ENABLE_Tpetra)
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
  IF (Stokhos_ENABLE_Experimental AND Stokhos_ENABLE_Belos)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/belos)
    SET(HEADERS ${HEADERS}
      belos/BelosPCETpetraAdapter.hpp
      )
  ENDIF()
  IF (Stokhos_ENABLE_Experimental AND Stokhos_ENABLE_Ifpack2 AND Stokhos_ENABLE_Belos)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ifpack2)
    SET(HEADERS ${HEADERS}
      ifpack2/Stokhos_Ifpack2.hpp
      )
    IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
      SET(SOURCES ${SOURCES}
        ifpack2/ifpack2_pce_instantiations.cpp
        )
    ENDIF()
  ENDIF()
ENDIF()

IF (Stokhos_ENABLE_KokkosCore)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos)
  SET(HEADERS ${HEADERS}
    kokkos/Stokhos_Multiply.hpp
    kokkos/Stokhos_CrsMatrix.hpp
    kokkos/Stokhos_BlockCrsMatrix.hpp
    kokkos/Stokhos_StochasticProductTensor.hpp
    kokkos/Stokhos_SymmetricDiagonalSpec.hpp
    kokkos/Stokhos_CrsProductTensor.hpp
    kokkos/Stokhos_TiledCrsProductTensor.hpp
    kokkos/Stokhos_SimpleTiledCrsProductTensor.hpp
    kokkos/Stokhos_CooProductTensor.hpp
    kokkos/Stokhos_FlatSparse3Tensor.hpp
    kokkos/Stokhos_FlatSparse3Tensor_kji.hpp
    kokkos/Stokhos_LinearSparse3Tensor.hpp
    kokkos/Stokhos_LexicographicBlockSparse3Tensor.hpp
    kokkos/Stokhos_TinyVec.hpp
    kokkos/Stokhos_Update.hpp
    kokkos/Stokhos_MatrixMarket.hpp
    )
  IF(Kokkos_ENABLE_Pthread)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/Threads)
    SET(HEADERS ${HEADERS}
      kokkos/Threads/Stokhos_Threads_CrsProductTensor.hpp
      )
  ENDIF()
  IF(Kokkos_ENABLE_OpenMP)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/OpenMP)
    SET(HEADERS ${HEADERS}
      kokkos/OpenMP/Stokhos_OpenMP_CrsProductTensor.hpp
      kokkos/OpenMP/Stokhos_OpenMP_MKL_CrsMatrix.hpp
      )
  ENDIF()
  IF(Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_Cuda)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/Cuda)
    SET(HEADERS ${HEADERS}
      kokkos/Cuda/Stokhos_Cuda_DeviceProp.hpp
      kokkos/Cuda/Stokhos_Cuda_WarpShuffle.hpp
      kokkos/Cuda/Stokhos_Cuda_CrsMatrix.hpp
      kokkos/Cuda/Stokhos_Cuda_BlockCrsMatrix.hpp
      kokkos/Cuda/Stokhos_Cuda_StochasticProductTensor.hpp
      kokkos/Cuda/Stokhos_Cuda_SymmetricDiagonalSpec.hpp
      kokkos/Cuda/Stokhos_Cuda_CrsProductTensor.hpp
      kokkos/Cuda/Stokhos_Cuda_TiledCrsProductTensor.hpp
      kokkos/Cuda/Stokhos_Cuda_SimpleTiledCrsProductTensor.hpp
      kokkos/Cuda/Stokhos_Cuda_CooProductTensor.hpp
      kokkos/Cuda/Stokhos_Cuda_FlatSparse3Tensor.hpp
      kokkos/Cuda/Stokhos_Cuda_FlatSparse3Tensor_kji.hpp
      kokkos/Cuda/Stokhos_Cuda_LinearSparse3Tensor.hpp
      kokkos/Cuda/Stokhos_Cuda_LexicographicBlockSparse3Tensor.hpp
      )
  ENDIF()
ENDIF()

# NOX
IF (Stokhos_ENABLE_NOX AND NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND
    NOX_ENABLE_EpetraExt)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/nox)
  SET(HEADERS ${HEADERS}
    epetra/Stokhos_GaussSeidelPreconditioner.hpp
    nox/NOX_Epetra_LinearSystem_SGJacobi.hpp
    nox/NOX_Epetra_LinearSystem_SGGS.hpp
    nox/NOX_Epetra_LinearSystem_MPBD.hpp
    )
  SET(SOURCES ${SOURCES}
    epetra/Stokhos_GaussSeidelPreconditioner.cpp
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

# Include cusp
IF (Stokhos_ENABLE_Cusp)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cusp)
  SET(HEADERS ${HEADERS}
    cusp/krylov/blockcg.h
    cusp/block_monitor.h
    cusp/block_multilevel.h
    cusp/MVmultiply.h
    cusp/detail/block_lu.h
    cusp/detail/device/MVmultiply.h
    cusp/detail/device/spmm/array2d.h
    cusp/detail/device/spmm/csr_vector.h
    cusp/detail/dispatch/MVmultiply.h
    cusp/precond/block_smoothed_aggregation.h
    cusp/relaxation/block_polynomial.h
    cusp/relaxation/block_jacobi.h
    )
  SET(SOURCES ${SOURCES} ${CUDA_SOURCES}
    cusp/krylov/blockcg.inl
    cusp/block_multilevel.inl
    cusp/detail/MVmultiply.inl
    cusp/precond/block_smoothed_aggregation.inl
    cusp/relaxation/block_polynomial.inl
    cusp/relaxation/block_jacobi.inl
    )
ENDIF()


#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stokhos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
