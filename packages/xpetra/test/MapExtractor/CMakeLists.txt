

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MapExtractorFactory_UnitTests
  SOURCES
    MapExtractorFactory_UnitTests
    ../Xpetra_UnitTests
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )