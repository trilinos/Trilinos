# -*- cmake -*-

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DrivenCavity
  SOURCES DrivenCavity.cpp
  COMM mpi
  NUM_MPI_PROCS 8
)
