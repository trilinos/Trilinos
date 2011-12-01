
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Workspace_test  
  SOURCES cxx_main.cpp 
  ARGS --verbose
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
