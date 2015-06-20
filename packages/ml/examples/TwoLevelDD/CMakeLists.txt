

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ml_2level_DD
  SOURCES ml_2level_DD.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
