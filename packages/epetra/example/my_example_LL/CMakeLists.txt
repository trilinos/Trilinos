INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  my_example_LL
  SOURCES cxx_main.cpp
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1-4
#  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
