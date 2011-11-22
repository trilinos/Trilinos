
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vector_random_fs
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 2-3
  ARGS "20 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NUM_MPI_PROCS 4
  ARGS "20 2"
  COMM mpi
  )
