
INCLUDE(TribitsAddExecutableAndTest)

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(TU_UT_OBJECTS_CUDA TpetraUtils_UnitTests_CU.cu)
  CUDA_COMPILE(TU_M_OBJECTS_CUDA TpetraUtils_merge_CU.cu)
  CUDA_COMPILE(TU_C_OBJECTS_CUDA TpetraUtils_congruent_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_UnitTests
  SOURCES
    TpetraUtils_UnitTests
    ${TU_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_merge
  SOURCES
    TpetraUtils_merge
    ${TU_M_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TpetraUtils_congruent
  SOURCES
    TpetraUtils_congruent
    ${TU_C_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
