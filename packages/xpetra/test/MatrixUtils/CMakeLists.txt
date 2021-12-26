
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixUtils_UnitTests
  SOURCES
    MatrixUtils_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
