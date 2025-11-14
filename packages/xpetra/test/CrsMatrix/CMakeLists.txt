



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
