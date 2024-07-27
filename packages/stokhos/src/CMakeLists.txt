INCLUDE(TrilinosCreateClientTemplateHeaders)

FUNCTION(PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(PROCESS_ETI_TEMPLATE)

FUNCTION(PROCESS_ENSEMBLE_ETI ENSEMBLE_SIZES ETI_CLASSES SRC_DIR BASE_NAME APPEND_BASE DEPLIBS DEPLIBS_OUT)
  # Template file name for each cpp file
  string(CONCAT TMPL_FILE ${SRC_DIR} "/" ${BASE_NAME} "_ETI_MP_Vector.tmpl")
  string(TOLOWER "${BASE_NAME}" LOWER_BASE_NAME)

  # List of devices to instantiate
  set(DEVICES "Serial" "Threads" "OpenMP" "Cuda")
  set(ADDED_LIBS "")
  foreach(DEVICE ${DEVICES})
    string(TOUPPER "${DEVICE}" UPPER_DEVICE)
    string(TOLOWER "${DEVICE}" LOWER_DEVICE)
    string(CONCAT INST_VAR "Tpetra_INST_" "${UPPER_DEVICE}")

    # Check if we are a device node and set a variable suitable for the
    # preprocessor accordingly
    string(COMPARE EQUAL "${DEVICE}" "Cuda" IS_DEVICE_NODE)

    # Create libraries for each enabled device instantation
    if(${${INST_VAR}})

      # Create library for each ensemble size
      foreach(ENSEMBLE_SIZE ${ENSEMBLE_SIZES})

        # Create header file for each ensemble size
        string(CONCAT BASE_HEADER "Stokhos_Tpetra_ETI_Helpers_MP_Vector_ENSEMBLE_SIZE_DEVICE.hpp")
        string(CONCAT BASE_HEADER_FULL "sacado/kokkos/vector/tpetra/" "${BASE_HEADER}")
        string(REPLACE "ENSEMBLE_SIZE" "${ENSEMBLE_SIZE}" HEADER_NAME "${BASE_HEADER}")
        string(REPLACE "DEVICE" "${DEVICE}" HEADER_NAME "${HEADER_NAME}")
        configure_file(${BASE_HEADER_FULL} ${HEADER_NAME})

        # Generate cpp file for each CLASS
        set(SRCS "")
        foreach(CLASS ${ETI_CLASSES})

          # find spearator of class name and (optional) conditions
          string(FIND ${CLASS} "-" DASHPOS)
          if ('${DASHPOS}' STREQUAL '-1')
            string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
            #MESSAGE("New class name = ${CLASS_NAME}")
            set (CONDITION_NAME "")
            set (CONDITION_NAME_END "")
          else()
            string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
            #MESSAGE("New class name = ${CLASS_NAME}")
            string(SUBSTRING ${CLASS} ${DASHPOS} -1 CONDITION_NAME)
            string(REPLACE "[" "(" CONDITION_NAME ${CONDITION_NAME})
            string(REPLACE "]" ")" CONDITION_NAME ${CONDITION_NAME})
            string(REPLACE "." " " CONDITION_NAME ${CONDITION_NAME})
            string(REPLACE "-" "" CONDITION_NAME ${CONDITION_NAME})
            string(REPLACE "?" "#" CONDITION_NAME ${CONDITION_NAME})
            string(STRIP CONDITION_NAME ${CONDITION_NAME})
            #MESSAGE("New condition name = ${CONDITION_NAME}")
            set (CONDITION_NAME_END "#endif")
          endif()

          string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS_NAME}")
          string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)

          if(APPEND_BASE)
            string(CONCAT CPP_FILE ${BASE_NAME} "_" "${CLASS_FILE_NAME}" "_MP_Vector_" "${ENSEMBLE_SIZE}" "_" "${DEVICE}" ".cpp")
          else()
            string(CONCAT CPP_FILE "${CLASS_FILE_NAME}" "_MP_Vector_" "${ENSEMBLE_SIZE}" "_" "${DEVICE}" ".cpp")
          endif()

          configure_file(${TMPL_FILE} ${CPP_FILE})
          list(APPEND SRCS ${CPP_FILE})
        endforeach()

        # Create library for each ensemble size
        string(CONCAT LIB_NAME "stokhos_" ${LOWER_BASE_NAME} "_mp_" "${ENSEMBLE_SIZE}" "_" "${LOWER_DEVICE}")
        TRIBITS_ADD_LIBRARY(
          "${LIB_NAME}"
          SOURCES ${SRCS}
          DEPLIBS ${DEPLIBS}
          )
        list(APPEND ADDED_LIBS "${LIB_NAME}")

      endforeach()
    endif()
  endforeach()
  set(${DEPLIBS_OUT} ${ADDED_LIBS} PARENT_SCOPE)
ENDFUNCTION(PROCESS_ENSEMBLE_ETI)

FUNCTION(PROCESS_PCE_ETI ETI_CLASSES EXTRA_ETI_CLASSES SRC_DIR BASE_NAME APPEND_BASE DEPLIBS DEPLIBS_OUT)
  # Template file name for each cpp file
  string(CONCAT TMPL_FILE ${SRC_DIR} "/" ${BASE_NAME} "_ETI_UQ_PCE.tmpl")
  string(TOLOWER "${BASE_NAME}" LOWER_BASE_NAME)

  # List of devices to instantiate
  set(DEVICES "Serial" "Threads" "OpenMP" "Cuda")
  set(ADDED_LIBS "")
  foreach(DEVICE ${DEVICES})
    string(TOUPPER "${DEVICE}" UPPER_DEVICE)
    string(TOLOWER "${DEVICE}" LOWER_DEVICE)
    string(CONCAT INST_VAR "Tpetra_INST_" "${UPPER_DEVICE}")

    # Check if we are a device node and set a variable suitable for the
    # preprocessor accordingly
    string(COMPARE EQUAL "${DEVICE}" "Cuda" IS_DEVICE_NODE)

    # Create libraries for each enabled device instantation
    if(${${INST_VAR}})

      # Create header file for device
      string(CONCAT BASE_HEADER "Stokhos_Tpetra_ETI_Helpers_UQ_PCE_DEVICE.hpp")
      string(CONCAT BASE_HEADER_FULL "sacado/kokkos/pce/tpetra/" "${BASE_HEADER}")
      string(REPLACE "DEVICE" "${DEVICE}" HEADER_NAME "${BASE_HEADER}")
      configure_file(${BASE_HEADER_FULL} ${HEADER_NAME})

      # Generate cpp file for each CLASS
      set(SRCS "")
      foreach(CLASS ${ETI_CLASSES})

        # find spearator of class name and (optional) conditions
        string(FIND ${CLASS} "-" DASHPOS)
        if ('${DASHPOS}' STREQUAL '-1')
          string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
          #MESSAGE("New class name = ${CLASS_NAME}")
          set (CONDITION_NAME "")
          set (CONDITION_NAME_END "")
        else()
          string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
          #MESSAGE("New class name = ${CLASS_NAME}")
          string(SUBSTRING ${CLASS} ${DASHPOS} -1 CONDITION_NAME)
          string(REPLACE "[" "(" CONDITION_NAME ${CONDITION_NAME})
          string(REPLACE "]" ")" CONDITION_NAME ${CONDITION_NAME})
          string(REPLACE "." " " CONDITION_NAME ${CONDITION_NAME})
          string(REPLACE "-" "" CONDITION_NAME ${CONDITION_NAME})
          string(REPLACE "?" "#" CONDITION_NAME ${CONDITION_NAME})
          string(STRIP CONDITION_NAME ${CONDITION_NAME})
          #MESSAGE("New condition name = ${CONDITION_NAME}")
          set (CONDITION_NAME_END "#endif")
        endif()

        string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS_NAME}")
        string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)

        if(APPEND_BASE)
          string(CONCAT CPP_FILE ${BASE_NAME} "_" "${CLASS_FILE_NAME}" "_UQ_PCE_" "${DEVICE}" ".cpp")
        else()
          string(CONCAT CPP_FILE "${CLASS_FILE_NAME}" "_UQ_PCE_" "${DEVICE}" ".cpp")
        endif()

        configure_file(${TMPL_FILE} ${CPP_FILE})
        list(APPEND SRCS ${CPP_FILE})
      endforeach()

      # Add in extra cpp files that we don't generate
      foreach(CLASS ${EXTRA_ETI_CLASSES})
        string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS}")
        string(CONCAT CPP_FILE "${SRC_DIR}" "/" ${BASE_NAME} "_" "${CLASS_FILE_NAME}" "_UQ_PCE_" "${DEVICE}" ".cpp")
        list(APPEND SRCS ${CPP_FILE})
      endforeach()

      # Create library for each device
      string(CONCAT LIB_NAME "stokhos_" ${LOWER_BASE_NAME} "_pce_" "${LOWER_DEVICE}")
      TRIBITS_ADD_LIBRARY(
        "${LIB_NAME}"
        SOURCES ${SRCS}
        DEPLIBS ${DEPLIBS}
        )
      list(APPEND ADDED_LIBS "${LIB_NAME}")
    endif()
  endforeach()
  set(${DEPLIBS_OUT} ${ADDED_LIBS} PARENT_SCOPE)
ENDFUNCTION(PROCESS_PCE_ETI)

#
# A) Package-specific configuration options
#

SET(HAVE_STOKHOS_DAKOTA ${${PACKAGE_NAME}_ENABLE_TriKota})

# OpenMP support -- Add our own option to turn it off
# without turning it off globally
SET(Stokhos_ENABLE_OpenMP_Default OFF)
IF(${PROJECT_NAME}_ENABLE_OpenMP)
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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add path to Dakota headers in installed directory
IF (${PACKAGE_NAME}_ENABLE_TriKota)
TRIBITS_INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_BINARY_DIR}/../../TriKota/Dakota/install/include
  )
ENDIF()


SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
  Stokhos_SmolyakSparseGridQuadrature.hpp
  Stokhos_SmolyakSparseGridQuadratureImp.hpp
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
  Stokhos_KL_ProductEigenPair.hpp
  Stokhos_KL_OneDExponentialEigenPair.hpp
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
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)
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
    epetra/Stokhos_SGModelEvaluatorBase.hpp
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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/kl)
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
IF(Kokkos_ENABLE_SERIAL)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/Serial)
  SET(HEADERS ${HEADERS}
    kokkos/Serial/Stokhos_Serial_CrsProductTensor.hpp
    )
ENDIF()
IF(Kokkos_ENABLE_THREADS)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/Threads)
  SET(HEADERS ${HEADERS}
    kokkos/Threads/Stokhos_Threads_CrsProductTensor.hpp
    )
ENDIF()
IF(Kokkos_ENABLE_OPENMP)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/OpenMP)
  SET(HEADERS ${HEADERS}
    kokkos/OpenMP/Stokhos_OpenMP_CrsProductTensor.hpp
    kokkos/OpenMP/Stokhos_OpenMP_MKL_CrsMatrix.hpp
    )
ENDIF()
IF(Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_CUDA)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kokkos/Cuda)
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

# NOX
IF (Stokhos_ENABLE_NOX AND NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND
    NOX_ENABLE_EpetraExt)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/nox)
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

# Include cusp
IF (Stokhos_ENABLE_Cusp)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cusp)
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

TRIBITS_ADD_LIBRARY(
  stokhos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# Sacado
IF (Stokhos_ENABLE_Sacado)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado)
  SET(SACADO_HEADERS "")
  SET(SACADO_SOURCES "")
  SET(STOKHOS_DEPLIBS stokhos)
  IF(Stokhos_ENABLE_PCE_Scalar_Type)
    LIST(APPEND SACADO_HEADERS
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
      )
    LIST(APPEND SACADO_SOURCES
      sacado/Sacado_PCE_OrthogPoly.cpp
      sacado/Sacado_ETPCE_OrthogPoly.cpp
      )
  ENDIF()
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos)
  LIST(APPEND SACADO_HEADERS
    sacado/kokkos/Stokhos_Sacado_Kokkos.hpp
    sacado/kokkos/Stokhos_Sacado_Kokkos_MathFunctions.hpp
    sacado/kokkos/Stokhos_MemoryTraits.hpp
    sacado/kokkos/Stokhos_KokkosTraits.hpp
    sacado/kokkos/Stokhos_StaticArrayTraits.hpp
    sacado/kokkos/Stokhos_DynArrayTraits.hpp
    sacado/kokkos/Stokhos_StaticStorage.hpp
    sacado/kokkos/Stokhos_StaticFixedStorage.hpp
    sacado/kokkos/Stokhos_DynamicStorage.hpp
    sacado/kokkos/Stokhos_DynamicStridedStorage.hpp
    sacado/kokkos/Stokhos_DynamicThreadedStorage.hpp
    sacado/kokkos/Stokhos_DynamicThreadedStorage_cuda.hpp
    sacado/kokkos/Stokhos_StorageHelpers.hpp
    sacado/kokkos/Stokhos_mpl_for_each.hpp
    sacado/kokkos/Stokhos_ViewStorage.hpp
    sacado/kokkos/Stokhos_Is_Constant.hpp
    sacado/kokkos/Kokkos_AnalyzeStokhosShape.hpp
    sacado/kokkos/Kokkos_View_Utils.hpp
    sacado/kokkos/Kokkos_View_Utils_Def.hpp
    )
  LIST(APPEND SACADO_SOURCES
    sacado/kokkos/Stokhos_MemoryTraits.cpp
    )
  IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector)
    LIST(APPEND SACADO_HEADERS
      sacado/kokkos/vector/Stokhos_Sacado_Kokkos_MP_Vector.hpp
      sacado/kokkos/vector/Stokhos_MP_Vector_MaskTraits.hpp
      sacado/kokkos/vector/Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp
      sacado/kokkos/vector/Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp
      sacado/kokkos/vector/Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp
      sacado/kokkos/vector/Sacado_MP_Vector.hpp
      sacado/kokkos/vector/Sacado_MP_VectorTraits.hpp
      sacado/kokkos/vector/Sacado_MP_ScalarTraitsImp.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_ops.hpp
      sacado/kokkos/vector/Sacado_MP_ExpressionTraits.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_SFS.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_SFS_unary_op_tmpl.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_SFS_unary_func_tmpl.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_SFS_binary_op_tmpl.hpp
      sacado/kokkos/vector/Sacado_MP_Vector_SFS_binary_func_tmpl.hpp
      sacado/kokkos/vector/Kokkos_View_MP_Vector_Fwd.hpp
      sacado/kokkos/vector/Kokkos_View_MP_Vector.hpp
      sacado/kokkos/vector/KokkosExp_View_MP_Vector_Contiguous.hpp
      sacado/kokkos/vector/Kokkos_View_MP_Vector_Interlaced.hpp
      sacado/kokkos/vector/Kokkos_View_MP_Vector_Utils.hpp
      sacado/kokkos/vector/Kokkos_Atomic_MP_Vector.hpp
      sacado/kokkos/vector/Kokkos_Parallel_MP_Vector.hpp
      sacado/kokkos/vector/KokkosBlas_MP_Vector.hpp
      sacado/kokkos/vector/KokkosBlas3_gemm_MP_Vector.hpp
      sacado/kokkos/vector/Teuchos_BLAS_MP_Vector.hpp
      sacado/kokkos/vector/Teuchos_LAPACK_MP_Vector.hpp
      sacado/kokkos/vector/Teuchos_SerialQRDenseSolver_MP_Vector.hpp
      ${CMAKE_CURRENT_BINARY_DIR}/Stokhos_Ensemble_Sizes.hpp
      )
    IF(HAVE_STOKHOS_ENSEMBLE_GEMV)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/vector/KokkosBlas2_gemv_MP_Vector.hpp
        )
    ENDIF()
    LIST(APPEND SACADO_SOURCES
      sacado/kokkos/vector/Kokkos_View_MP_Vector.cpp
      )

    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/Fad)
    LIST(APPEND SACADO_HEADERS
      sacado/kokkos/vector/Fad/Sacado_Fad_ExprSpec_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_Expr_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_GeneralFad_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_Ops_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_DFad_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_SFad_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_SLFad_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_ViewFad_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_Exp_MP_Vector.hpp
      sacado/kokkos/vector/Fad/Sacado_Fad_MP_Vector.hpp
      )
  ENDIF()
  IF(Stokhos_ENABLE_PCE_Scalar_Type)
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce)
    LIST(APPEND SACADO_HEADERS
      sacado/kokkos/pce/Stokhos_Sacado_Kokkos_UQ_PCE.hpp
      sacado/kokkos/pce/Sacado_UQ_PCE.hpp
      sacado/kokkos/pce/Sacado_UQ_PCE_Imp.hpp
      sacado/kokkos/pce/Sacado_UQ_PCE_Traits.hpp
      sacado/kokkos/pce/Sacado_UQ_PCE_ScalarTraitsImp.hpp
      sacado/kokkos/pce/Kokkos_View_UQ_PCE_Fwd.hpp
      sacado/kokkos/pce/Kokkos_View_UQ_PCE.hpp
      sacado/kokkos/pce/KokkosExp_View_UQ_PCE_Contiguous.hpp
      sacado/kokkos/pce/Kokkos_View_UQ_PCE_Utils.hpp
      sacado/kokkos/pce/Kokkos_Atomic_UQ_PCE.hpp
      sacado/kokkos/pce/Teuchos_BLAS_UQ_PCE.hpp
      sacado/kokkos/pce/Teuchos_LAPACK_UQ_PCE.hpp
      sacado/kokkos/pce/Teuchos_SerialQRDenseSolver_UQ_PCE.hpp
      )
  ENDIF()

  IF (Stokhos_ENABLE_KokkosKernels)
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/linalg)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/vector/linalg/Kokkos_MV_MP_Vector.hpp
        sacado/kokkos/vector/linalg/Kokkos_Blas1_MP_Vector.hpp
        sacado/kokkos/vector/linalg/Kokkos_CrsMatrix_MP_Vector.hpp
        sacado/kokkos/vector/linalg/Kokkos_CrsMatrix_MP_Vector_Cuda.hpp
        sacado/kokkos/vector/linalg/Kokkos_ArithTraits_MP_Vector.hpp
        sacado/kokkos/vector/linalg/Kokkos_InnerProductSpaceTraits_MP_Vector.hpp
        )
    ENDIF()
    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/linalg)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/pce/linalg/Kokkos_MV_UQ_PCE.hpp
        sacado/kokkos/pce/linalg/Kokkos_Blas1_UQ_PCE.hpp
        sacado/kokkos/pce/linalg/Kokkos_CrsMatrix_UQ_PCE.hpp
        sacado/kokkos/pce/linalg/Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp
        sacado/kokkos/pce/linalg/Kokkos_ArithTraits_UQ_PCE.hpp
        sacado/kokkos/pce/linalg/Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp
        )
    ENDIF()
  ENDIF()

  IF (Stokhos_ENABLE_TeuchosKokkosComm)
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/mpicomm)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/vector/mpicomm/Kokkos_TeuchosCommAdapters_MP_Vector.hpp
        )
    ENDIF()
    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/mpicomm)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/pce/mpicomm/Kokkos_TeuchosCommAdapters_UQ_PCE.hpp
        )
    ENDIF()
  ENDIF()

  IF (Stokhos_ENABLE_Kokkos)
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/algorithms)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/vector/algorithms/Kokkos_Random_MP_Vector.hpp
        )
    ENDIF()
    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/algorithms)
      LIST(APPEND SACADO_HEADERS
        sacado/kokkos/pce/algorithms/Kokkos_Random_UQ_PCE.hpp
        )
    ENDIF()
  ENDIF()

  TRIBITS_ADD_LIBRARY(
    stokhos_sacado
    HEADERS ${SACADO_HEADERS}
    SOURCES ${SACADO_SOURCES}
    DEPLIBS ${STOKHOS_DEPLIBS}
    )
  LIST(APPEND STOKHOS_DEPLIBS stokhos_sacado)

  IF (Stokhos_ENABLE_Tpetra)
    SET(TPETRA_HEADERS
      sacado/kokkos/Stokhos_Tpetra_Utilities.hpp)
    SET(TPETRA_CPP_SOURCES
      Stokhos_Dummy.cpp)
    SET(STOKHOS_TPETRA_DEPLIBS "")
    SET(TPETRA_ETI_CLASSES
        MultiVector
        Vector
        RowMatrix
        RowMatrixTransposer
        DistObject
        CrsMatrix
        FECrsMatrix
        FEMultiVector
        BlockVector
        BlockMultiVector
        BlockCrsMatrix
        Details::getDiagCopyWithoutOffsets
        Details::packCrsMatrix
        Details::unpackCrsMatrixAndCombine
        replaceDiagonalCrsMatrix
        BlockCrsMatrix_Helpers
        )
      SET(TPETRAEXT_ETI_CLASSES
        MatrixMatrix
        TripleMatrixMultiply
        )
      SET(TPETRA_SD_ETI_CLASSES
        LocalCrsMatrixOperator
        )
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/tpetra)
      APPEND_GLOB(TPETRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/tpetra/*.hpp)
      IF(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
        IF (NOT ("${TPETRA_SD_ETI_CLASSES}" STREQUAL ""))
          PROCESS_ENSEMBLE_ETI(
            "${Stokhos_ETI_Ensemble_Sizes}"
            "${TPETRA_SD_ETI_CLASSES}"
            "sacado/kokkos/vector/tpetra"
            "Tpetra_SD"
            TRUE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()

        # Tpetra depends on Tpetra_SD
        SET(STOKHOS_TPETRA_DEPS ${STOKHOS_DEPLIBS} ${STOKHOS_TPETRA_DEPLIBS})
        PROCESS_ENSEMBLE_ETI(
          "${Stokhos_ETI_Ensemble_Sizes}"
          "${TPETRA_ETI_CLASSES}"
          "sacado/kokkos/vector/tpetra"
          "Tpetra"
          TRUE
          "${STOKHOS_TPETRA_DEPS}"
          DEPLIBS_OUT)
        LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
        # TpetraExt depends on Tpetra
        SET(STOKHOS_TPETRAEXT_DEPS ${STOKHOS_DEPLIBS} ${STOKHOS_TPETRA_DEPLIBS})
        PROCESS_ENSEMBLE_ETI(
          "${Stokhos_ETI_Ensemble_Sizes}"
          "${TPETRAEXT_ETI_CLASSES}"
          "sacado/kokkos/vector/tpetra"
          "TpetraExt"
          TRUE
          "${STOKHOS_TPETRAEXT_DEPS}"
          DEPLIBS_OUT)
        LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
      ENDIF()
    ENDIF()
    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/tpetra)
      APPEND_GLOB(TPETRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/tpetra/*.hpp)
      IF(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
        # We have manually written cpp files for CrsMatrix
        LIST(REMOVE_ITEM TPETRA_ETI_CLASSES CrsMatrix)
        SET(EXTRA_TPETRA_ETI_CLASSES
          CrsMatrix
          )
        IF (NOT ("${TPETRA_SD_ETI_CLASSES}" STREQUAL ""))
          PROCESS_PCE_ETI(
            "${TPETRA_SD_ETI_CLASSES}"
            ""
            "sacado/kokkos/pce/tpetra"
            "Tpetra_SD"
            TRUE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()

        # Tpetra depends on Tpetra_SD
        SET(STOKHOS_TPETRA_DEPS ${STOKHOS_DEPLIBS} ${STOKHOS_TPETRA_DEPLIBS})
        PROCESS_PCE_ETI(
          "${TPETRA_ETI_CLASSES}"
          "${EXTRA_TPETRA_ETI_CLASSES}"
          "sacado/kokkos/pce/tpetra"
          "Tpetra"
          TRUE
          "${STOKHOS_TPETRA_DEPS}"
          DEPLIBS_OUT)
        LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
        # TpetraExt depends on Tpetra
        SET(STOKHOS_TPETRAEXT_DEPS ${STOKHOS_DEPLIBS} ${STOKHOS_TPETRA_DEPLIBS})
        PROCESS_PCE_ETI(
          "${TPETRAEXT_ETI_CLASSES}"
          ""
          "sacado/kokkos/pce/tpetra"
          "TpetraExt"
          TRUE
          "${STOKHOS_TPETRAEXT_DEPS}"
          DEPLIBS_OUT)
        LIST(APPEND STOKHOS_TPETRA_DEPLIBS ${DEPLIBS_OUT})
      ENDIF()
    ENDIF()

    # Belos is a required dependency for Ifpack2, so we need the Belos
    # adapters if either is enabled.
    IF (Stokhos_ENABLE_Belos OR Stokhos_ENABLE_Ifpack2)
      IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/belos)
        LIST(APPEND TPETRA_HEADERS
          sacado/kokkos/vector/belos/Belos_Tpetra_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_TpetraAdapter_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_SolverManager_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_PseudoBlockCGIter_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_StatusTest_GenResNorm_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_StatusTest_ImpResNorm_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_PseudoBlockGmresIter_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_DGKS_OrthoManager_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_ICGS_OrthoManager_MP_Vector.hpp
          sacado/kokkos/vector/belos/Belos_IMGS_OrthoManager_MP_Vector.hpp
          )
      ENDIF()
      IF(Stokhos_ENABLE_PCE_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/belos)
        LIST(APPEND TPETRA_HEADERS
          sacado/kokkos/pce/belos/Belos_TpetraAdapter_UQ_PCE.hpp
          )
      ENDIF()
    ENDIF()

    TRIBITS_ADD_LIBRARY(
      stokhos_tpetra
      HEADERS ${TPETRA_HEADERS}
      SOURCES ${TPETRA_CPP_SOURCES}
      DEPLIBS ${STOKHOS_DEPLIBS}
      ADDED_LIB_TARGET_NAME_OUT STOKHOS_TPETRA_LIBNAME
      )
    SET_PROPERTY(TARGET ${STOKHOS_TPETRA_LIBNAME} APPEND PROPERTY LINKER_LANGUAGE CXX)
    LIST(APPEND STOKHOS_TPETRA_DEPLIBS stokhos_tpetra)
    LIST(APPEND STOKHOS_DEPLIBS ${STOKHOS_TPETRA_DEPLIBS})

    IF (Stokhos_ENABLE_Xpetra)
      SET(XPETRA_HEADERS "")
      SET(XPETRA_CPP_SOURCES Stokhos_Dummy.cpp)
      SET(STOKHOS_XPETRA_DEPLIBS "")

      # Load list of Xpetra ETI classes
      #
      # CAG: We make sure that no Thyra classes are used since Stokhos currently does not ETI them.
      #      We disabled the Xpetra-Thyra interface in
      #       src/sacado/kokkos/pce/xpetra
      #       src/sacado/kokkos/vector/xpetra
      GLOBAL_SET(XPETRA_SC_LO_GO_NO_ETI_CLASSES  )
      IF (Xpetra_ENABLE_EXPLICIT_INSTANTIATION)
        INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../../xpetra/src/Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)
        INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../../xpetra/src/Utils/ExplicitInstantiation/ETI_MULTIVECTOR_classes.cmake)
        LIST(REMOVE_ITEM XPETRA_SC_LO_GO_NO_ETI_CLASSES Xpetra::ThyraUtils-.?if.defined[HAVE_XPETRA_THYRA])
        SET(XPETRA_ALL_ETI_CLASSES ${XPETRA_SC_LO_GO_NO_ETI_CLASSES} ${XPETRA_MULTIVECTOR_ETI_CLASSES})
      ENDIF()

      IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/xpetra)
        APPEND_GLOB(XPETRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/xpetra/*.hpp)
        IF (Xpetra_ENABLE_EXPLICIT_INSTANTIATION)
          PROCESS_ENSEMBLE_ETI(
            "${Stokhos_ETI_Ensemble_Sizes}"
            "${XPETRA_ALL_ETI_CLASSES}"
            "sacado/kokkos/vector/xpetra"
            "Xpetra"
            FALSE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_XPETRA_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      IF(Stokhos_ENABLE_PCE_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/xpetra)
        APPEND_GLOB(XPETRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/xpetra/*.hpp)
        IF(Xpetra_ENABLE_EXPLICIT_INSTANTIATION)
          PROCESS_PCE_ETI(
            "${XPETRA_ALL_ETI_CLASSES}"
            ""
            "sacado/kokkos/pce/xpetra"
            "Xpetra"
            FALSE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_XPETRA_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      TRIBITS_ADD_LIBRARY(
        stokhos_xpetra
        HEADERS ${XPETRA_HEADERS}
        SOURCES ${XPETRA_CPP_SOURCES}
        DEPLIBS ${STOKHOS_DEPLIBS}
        ADDED_LIB_TARGET_NAME_OUT STOKHOS_XPETRA_LIBNAME
        )
      SET_PROPERTY(TARGET ${STOKHOS_XPETRA_LIBNAME} APPEND PROPERTY LINKER_LANGUAGE CXX)
      LIST(APPEND STOKHOS_DEPLIBS stokhos_xpetra)
      LIST(APPEND STOKHOS_DEPLIBS ${STOKHOS_XPETRA_DEPLIBS})
    ENDIF()

    IF (Stokhos_ENABLE_Amesos2)
      SET(AMESOS2_HEADERS "")
      SET(AMESOS2_SOURCES Stokhos_Dummy.cpp)
      SET(AMESOS2_ETI_CLASSES "")
      IF(Amesos2_ENABLE_Basker)
        LIST(APPEND AMESOS2_ETI_CLASSES
          Basker
          )
      ENDIF()
      SET(STOKHOS_AMESOS2_DEPLIBS "")
      IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/amesos2)
        APPEND_GLOB(AMESOS2_HEADERS
          sacado/kokkos/vector/amesos2/*.hpp
          )
        APPEND_GLOB(AMESOS2_SOURCES
          sacado/kokkos/vector/amesos2/*.cpp
          )
        IF(AMESOS2_ETI_CLASSES) # Only process if we have classes to instantiate
          PROCESS_ENSEMBLE_ETI(
            "${Stokhos_ETI_Ensemble_Sizes}"
            "${AMESOS2_ETI_CLASSES}"
            "sacado/kokkos/vector/amesos2"
            "Amesos2"
            TRUE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_AMESOS2_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      IF(Stokhos_ENABLE_PCE_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/amesos2)
        APPEND_GLOB(AMESOS2_HEADERS
          sacado/kokkos/pce/amesos2/*.hpp
          )
        APPEND_GLOB(AMESOS2_SOURCES
          sacado/kokkos/pce/amesos2/*.cpp
          )
      ENDIF()
      TRIBITS_ADD_LIBRARY(
        stokhos_amesos2
        HEADERS ${AMESOS2_HEADERS}
        SOURCES ${AMESOS2_SOURCES}
        DEPLIBS ${STOKHOS_DEPLIBS}
        ADDED_LIB_TARGET_NAME_OUT STOKHOS_AMESOS2_LIBNAME
        )
      SET_PROPERTY(TARGET ${STOKHOS_AMESOS2_LIBNAME} APPEND PROPERTY LINKER_LANGUAGE CXX)
      LIST(APPEND STOKHOS_DEPLIBS stokhos_amesos2)
      LIST(APPEND STOKHOS_DEPLIBS ${STOKHOS_AMESOS2_DEPLIBS})
    ENDIF()

    IF (Stokhos_ENABLE_Ifpack2)
      SET(IFPACK2_HEADERS "")
      SET(IFPACK2_CPP_SOURCES Stokhos_Dummy.cpp)
      SET(STOKHOS_IFPACK2_DEPLIBS "")
      IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/ifpack2)
        APPEND_GLOB(IFPACK2_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/ifpack2/*.hpp)
        IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
          APPEND_GLOB(IFPACK2_CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/ifpack2/*.cpp)
          PROCESS_ENSEMBLE_ETI(
            "${Stokhos_ETI_Ensemble_Sizes}"
            "${IFPACK2_ETI_CLASSES}"
            "sacado/kokkos/vector/ifpack2"
            "Ifpack2"
            TRUE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_IFPACK2_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      IF(Stokhos_ENABLE_PCE_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/ifpack2)
        APPEND_GLOB(IFPACK2_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/ifpack2/*.hpp)
        IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
          APPEND_GLOB(IFPACK2_CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/ifpack2/*.cpp)
          PROCESS_PCE_ETI(
            "${IFPACK2_ETI_CLASSES}"
            ""
            "sacado/kokkos/pce/ifpack2"
            "Ifpack2"
            TRUE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_IFPACK2_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      TRIBITS_ADD_LIBRARY(
        stokhos_ifpack2
        HEADERS ${IFPACK2_HEADERS}
        SOURCES ${IFPACK2_CPP_SOURCES}
        DEPLIBS ${STOKHOS_DEPLIBS}
        ADDED_LIB_TARGET_NAME_OUT STOKHOS_IFPACK2_LIBNAME
        )
      SET_PROPERTY(TARGET ${STOKHOS_IFPACK2_LIBNAME} APPEND PROPERTY LINKER_LANGUAGE CXX)
      LIST(APPEND STOKHOS_DEPLIBS stokhos_ifpack2)
      LIST(APPEND STOKHOS_DEPLIBS ${STOKHOS_IFPACK2_DEPLIBS})
    ENDIF()

    IF (Stokhos_ENABLE_MueLu)
      IF (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
        message(WARNING "Enabling MueLu in Stokhos when compiling with CUDA relocatable device code may cause infinite compile times")
      ENDIF()
      SET(MUELU_HEADERS "")
      SET(MUELU_CPP_SOURCES Stokhos_Dummy.cpp)
      SET(STOKHOS_MUELU_DEPLIBS "")

      # Load list of MueLu ETI classes
      GLOBAL_SET(MUELU_SC_LO_GO_NO_ETI_CLASSES  )
      IF (MueLu_ENABLE_EXPLICIT_INSTANTIATION)
        INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../../muelu/src/Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

        # We've deprecated MueLu::MLParameterListInterpreter, so we'll drop Stokhos Support
        LIST(REMOVE_ITEM MUELU_SC_LO_GO_NO_ETI_CLASSES "MueLu::MLParameterListInterpreter-.?if.defined[HAVE_MUELU_DEPRECATED_CODE]" )

      ENDIF()

      IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/muelu)
        APPEND_GLOB(MUELU_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/vector/muelu/*.hpp)
        IF (MueLu_ENABLE_EXPLICIT_INSTANTIATION)
          PROCESS_ENSEMBLE_ETI(
            "${Stokhos_ETI_Ensemble_Sizes}"
            "${MUELU_SC_LO_GO_NO_ETI_CLASSES}"
            "sacado/kokkos/vector/muelu"
            "MueLu"
            FALSE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_MUELU_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      IF(Stokhos_ENABLE_PCE_Scalar_Type)
        TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/muelu)
        APPEND_GLOB(MUELU_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/sacado/kokkos/pce/muelu/*.hpp)
        IF(MueLu_ENABLE_EXPLICIT_INSTANTIATION)
          PROCESS_PCE_ETI(
            "${MUELU_SC_LO_GO_NO_ETI_CLASSES}"
            ""
            "sacado/kokkos/pce/muelu"
            "MueLu"
            FALSE
            "${STOKHOS_DEPLIBS}"
            DEPLIBS_OUT)
          LIST(APPEND STOKHOS_MUELU_DEPLIBS ${DEPLIBS_OUT})
        ENDIF()
      ENDIF()
      TRIBITS_ADD_LIBRARY(
        stokhos_muelu
        HEADERS ${MUELU_HEADERS}
        SOURCES ${MUELU_CPP_SOURCES}
        DEPLIBS ${STOKHOS_DEPLIBS}
        ADDED_LIB_TARGET_NAME_OUT STOKHOS_MUELU_LIBNAME
        )
      SET_PROPERTY(TARGET ${STOKHOS_MUELU_LIBNAME} APPEND PROPERTY LINKER_LANGUAGE CXX)
      LIST(APPEND STOKHOS_DEPLIBS stokhos_muelu)
      LIST(APPEND STOKHOS_DEPLIBS ${STOKHOS_MUELU_DEPLIBS})
    ENDIF()

  ENDIF()

ENDIF()
