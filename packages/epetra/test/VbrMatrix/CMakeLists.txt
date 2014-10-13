
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  VbrMatrix_test
  SOURCES cxx_main.cpp
  DEPLIBS  epetratest
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
