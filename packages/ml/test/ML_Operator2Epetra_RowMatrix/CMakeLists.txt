

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ML_Operator2Epetra_RowMatrix
  SOURCES ML_Operator2Epetra_RowMatrix.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
