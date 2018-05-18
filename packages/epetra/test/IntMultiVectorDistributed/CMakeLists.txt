
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IntMultiVectorDistributed_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM mpi
  NUM_MPI_PROCS 4
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
