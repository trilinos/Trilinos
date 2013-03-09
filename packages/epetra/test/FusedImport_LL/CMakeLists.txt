INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FusedImport_test_LL
  SOURCES cxx_main.cpp 
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  PASS_REGULAR_EXPRESSION "tests PASSED"
  )
