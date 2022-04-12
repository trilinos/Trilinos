

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MapExtractorFactory_UnitTests
  SOURCES
    MapExtractorFactory_UnitTests.cpp
    ../Xpetra_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
