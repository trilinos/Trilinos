
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ml_multiple_rhs
  SOURCES ml_multiple_rhs.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
