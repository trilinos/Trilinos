
IF(Tpetra_ENABLE_DEPRECATED_CODE)
  # MultiVectorFiller may use Kokkos data structures.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVectorFiller
    SOURCES MultiVectorFiller.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  # This test is only informative with multiple MPI processes,
  # but it should at least compile and run correctly without MPI.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVectorFiller2
    SOURCES
    MultiVectorFiller2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 2-4
    ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
    STANDARD_PASS_OUTPUT
    )
ENDIF()

