TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PointToBlockDiagPermute_UnitTests
  SOURCES
    PointToBlockDiagPermute_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )
