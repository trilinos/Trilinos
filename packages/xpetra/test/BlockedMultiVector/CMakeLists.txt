
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockedMultiVector_UnitTests
  SOURCES
    BlockedMultiVector_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
