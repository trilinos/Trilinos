
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LegendreBasisUnitTest
  SOURCES Stokhos_LegendreBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NormalizedLegendreBasisUnitTest
  SOURCES Stokhos_NormalizedLegendreBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HermiteBasisUnitTest
  SOURCES Stokhos_HermiteBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NormalizedHermiteBasisUnitTest
  SOURCES Stokhos_NormalizedHermiteBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  JacobiBasisUnitTest
  SOURCES Stokhos_JacobiBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  DerivExpansionUnitTest
#  SOURCES Stokhos_DerivExpansionUnitTest.cpp
#  COMM serial mpi
#  STANDARD_PASS_OUTPUT
#  NUM_MPI_PROCS 1
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  QuadExpansionUnitTest
  SOURCES Stokhos_QuadExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AlgebraicExpansionUnitTest
  SOURCES Stokhos_AlgebraicExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConstantExpansionUnitTest
  SOURCES Stokhos_ConstantExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DivisionOperatorUnitTest
  SOURCES Stokhos_DivisionOperatorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

IF(TPL_ENABLE_ForUQTK)
 
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ForUQTKExpansionUnitTest
    SOURCES Stokhos_ForUQTKExpansionUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StieltjesUnitTest
  SOURCES Stokhos_StieltjesUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LanczosUnitTest
  SOURCES Stokhos_LanczosUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GramSchmidtUnitTest
  SOURCES Stokhos_GramSchmidtUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Sparse3TensorUnitTest
  SOURCES Stokhos_Sparse3TensorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExponentialRandomFieldUnitTest
  SOURCES Stokhos_ExponentialRandomFieldUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixFreeOperatorUnitTest
  SOURCES Stokhos_MatrixFreeOperatorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InterlacedOpUnitTest 
  SOURCES Stokhos_InterlacedOpUnitTest.cpp
          Stokhos_InterlacedTestSupport.cpp
          Stokhos_InterlacedTestSupport.hpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasisInteractionGraphUnitTest 
  SOURCES Stokhos_BasisInteractionGraphUnitTest.cpp
          Stokhos_InterlacedTestSupport.cpp
          Stokhos_InterlacedTestSupport.hpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AdaptivityToolsUnitTest 
  SOURCES Stokhos_AdaptivityToolsUnitTest.cpp
          Stokhos_InterlacedTestSupport.cpp
          Stokhos_InterlacedTestSupport.hpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BuildColBasisUnitTest 
  SOURCES Stokhos_BuildColBasisUnitTest.cpp
          Stokhos_InterlacedTestSupport.cpp
          Stokhos_InterlacedTestSupport.hpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InterlacedMapUnitTest 
  SOURCES Stokhos_InterlacedMapUnitTest.cpp
          Stokhos_InterlacedTestSupport.cpp
          Stokhos_InterlacedTestSupport.hpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )

IF(Stokhos_ENABLE_TriKota)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SparseGridQuadratureUnitTest
    SOURCES Stokhos_SparseGridQuadratureUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

ENDIF()

IF(Stokhos_ENABLE_Sacado)

  INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../sacado/test/utils)

  SET(UTILS_SOURCES 
    ${PACKAGE_SOURCE_DIR}/../sacado/test/utils/Sacado_Random.hpp
    ${PACKAGE_SOURCE_DIR}/../sacado/test/utils/Sacado_RandomImp.hpp
  )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCEUnitTest
    SOURCES Stokhos_SacadoPCEUnitTest.cpp Stokhos_SacadoPCEUnitTest.hpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoETPCEUnitTest
    SOURCES Stokhos_SacadoETPCEUnitTest.cpp Stokhos_SacadoPCEUnitTest.hpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorUnitTest
    SOURCES Stokhos_SacadoVectorUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCESerializationTests
    SOURCES Stokhos_SacadoPCESerializationTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCECommTests
    SOURCES Stokhos_SacadoPCECommTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorSerializationTests
    SOURCES Stokhos_SacadoVectorSerializationTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorCommTests
    SOURCES Stokhos_SacadoVectorCommTests.cpp ${UTILS_SOURCES}
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  IF (Stokhos_ENABLE_KokkosArray)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SacadoMPVectorUnitTest
      SOURCES Stokhos_SacadoMPVectorUnitTest.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SacadoMPVectorSerializationTests
      SOURCES Stokhos_SacadoMPVectorSerializationTests.cpp ${UTILS_SOURCES}
      ARGS 
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SacadoMPVectorCommTests
      SOURCES Stokhos_SacadoMPVectorCommTests.cpp ${UTILS_SOURCES}
      ARGS 
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
  ENDIF()

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LogNormalUnitTest
  SOURCES Stokhos_LogNormalUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SDMUtilsUnitTest
  SOURCES Stokhos_SDMUtilsUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ProductBasisUtilsUnitTest
  SOURCES Stokhos_ProductBasisUtilsUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorProductBasisUnitTest
  SOURCES Stokhos_TensorProductBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TotalOrderBasisUnitTest
  SOURCES Stokhos_TotalOrderBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   SmolyakBasisUnitTest
#   SOURCES Stokhos_SmolyakBasisUnitTest.cpp
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   NUM_MPI_PROCS 1
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorProductPseudoSpectralOperatorUnitTest
  SOURCES Stokhos_TensorProductPseudoSpectralOperatorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )
