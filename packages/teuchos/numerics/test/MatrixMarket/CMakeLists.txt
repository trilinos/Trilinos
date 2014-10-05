INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Raw_InOutTest
  SOURCES
    MatrixMarket_Raw_InOutTest.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
