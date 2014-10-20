

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MLP_Maxwell
  SOURCES MultiLevelPreconditioner_Maxwell.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
