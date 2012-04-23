INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SimpleLongLongTest
  SOURCES cxx_main.cpp 
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  FAIL_REGULAR_EXPRESSION "Non zero error code"
  )
