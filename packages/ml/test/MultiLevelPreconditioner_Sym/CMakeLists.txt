

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MLP_Sym
  SOURCES MultiLevelPreconditioner_Sym.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
