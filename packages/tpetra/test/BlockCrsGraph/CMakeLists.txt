
INCLUDE(TribitsAddExecutableAndTest)

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(BLOCKCRSGRAPH_UNIT_TEST_OBJECTS_CUDA BlockCrsGraph_UnitTests_CU.cu)
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockCrsGraph_UnitTests
  SOURCES
    BlockCrsGraph_UnitTests
    ${BLOCKCRSGRAPH_UNIT_TEST_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
