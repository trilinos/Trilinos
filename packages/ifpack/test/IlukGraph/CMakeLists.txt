
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  IlukGraph_test
  SOURCES cxx_main.cpp
  ARGS "-v"
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
