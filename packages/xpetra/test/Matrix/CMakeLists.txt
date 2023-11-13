


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Matrix_UnitTests
  SOURCES
    Matrix_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
