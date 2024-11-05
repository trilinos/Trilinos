TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LinearProblem_UnitsTests
  SOURCES
    LinearProblem_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
