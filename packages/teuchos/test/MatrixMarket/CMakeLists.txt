INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Raw_InOutTest
  SOURCES
    MatrixMarket_Raw_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS teuchos
  )
