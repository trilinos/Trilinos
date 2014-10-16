IF (${PACKAGE_NAME}_ENABLE_CLASSIC_VBR)
  IF(Kokkos_ENABLE_Cuda)
    # Add subpackage include path dependencies for CUDA_COMPILE().
    CUDA_INCLUDE_DIRECTORIES( "." )
    CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
    CUDA_COMPILE(BLOCKCRSGRAPH_UNIT_TEST_OBJECTS_CUDA BlockCrsGraph_UnitTests_CU.cu)
  ENDIF()

  # FIXME (mfh 10 Jul 2014) The two TEST_COMPARE_ARRAYS in these tests
  # don't pass on all the MPI processes, but because they pass on
  # Process 0, the tests report as passed.  Thus, I am disabling the
  # tests for now.  VBR in general could use a lot of work and shouldn't
  # be trusted until we revamp it.
  #
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   BlockCrsGraph_UnitTests
  #   SOURCES
  #     BlockCrsGraph_UnitTests
  #     ${BLOCKCRSGRAPH_UNIT_TEST_OBJECTS_CUDA}
  #     ${TEUCHOS_STD_UNIT_TEST_MAIN}
  #   COMM serial mpi
  #   STANDARD_PASS_OUTPUT
  #   )
ENDIF ()
