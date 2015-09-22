

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MLP_Restart
  SOURCES MultiLevelPreconditioner_Restart.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
