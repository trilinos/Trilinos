



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests.cpp
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixUtils_UnitTests
  SOURCES
    CrsMatrixUtils_UnitTests.cpp
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
