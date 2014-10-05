
INCLUDE(TribitsAddExecutableAndTest)

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(DIRECTORY_UT_OBJECTS_CUDA Directory_UnitTests_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Directory_UnitTests
  SOURCES
    Directory_UnitTests
    ${DIRECTORY_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
