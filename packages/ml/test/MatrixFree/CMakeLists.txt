

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixFree
  SOURCES MatrixFree.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
