
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LinearSolverSetupFailure
  SOURCES
    LinearSolverSetupFailure.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
