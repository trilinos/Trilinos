
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Communicator
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
