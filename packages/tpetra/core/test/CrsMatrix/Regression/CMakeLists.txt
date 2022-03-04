IF (Tpetra_INST_DOUBLE)
  # This test explicitly checks for Scalar == double
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Albany182
    SOURCES
      Albany182.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    ARGS "--globally-reduce-test-result"
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
    )
ENDIF()
# --output-show-proc-rank --output-to-root-rank-only=-1
