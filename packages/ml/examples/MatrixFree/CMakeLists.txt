
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ml_1D
  SOURCES ml_1D.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
