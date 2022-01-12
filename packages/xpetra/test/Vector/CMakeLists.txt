
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Vector_UnitTests2
  SOURCES
    Vector_UnitTests2.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  )
