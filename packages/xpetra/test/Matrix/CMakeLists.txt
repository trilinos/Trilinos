
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Matrix_UnitTests
  SOURCES
    Matrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
