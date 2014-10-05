
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEVbrMatrix_test
  SOURCES cxx_main.cpp ExecuteTestProblems.cpp
  DEPLIBS  epetratest
  ARGS -v
  COMM mpi serial
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
