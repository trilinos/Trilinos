
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug_5988_FECrs_buildNonlocalGraph_LL
  SOURCES cxx_main.cpp
  ARGS -v
  COMM mpi
  NUM_MPI_PROCS 2-4
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
