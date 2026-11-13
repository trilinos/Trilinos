

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TriDiMatrix_test 
  SOURCES TriDiMatrix_UnitTest.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
