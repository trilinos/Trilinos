TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExprEvalTest
  SOURCES ExprEvalTest.cpp UnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
