
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixUtils_UnitTests
  SOURCES
    MatrixUtils_UnitTests
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
