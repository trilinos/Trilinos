

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(DP_T_OBJECTS_CUDA DefaultPlatform_test_CU.cu)
  CUDA_COMPILE(P_UT_OBJECTS_CUDA Platform_UnitTests_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultPlatform_test
  SOURCES
    DefaultPlatform_test.cpp
    ${DP_T_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Platform_UnitTests
  SOURCES
    Platform_UnitTests.cpp
    ${P_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
