
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vector_random
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 1-3
  ARGS "1000 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random
  NUM_MPI_PROCS 4
  ARGS "1000 2"
  COMM mpi
  )
