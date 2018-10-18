


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Matrix_UnitTests
  SOURCES
    Matrix_UnitTests
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
