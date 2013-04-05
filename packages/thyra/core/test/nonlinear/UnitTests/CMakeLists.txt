INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultNonlinearSolverBuilder
  SOURCES
    DefaultNonlinearSolverBuilder_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultFiniteDifferenceModelEvaluator
  SOURCES
    DefaultFiniteDifferenceModelEvaluator_UnitTests.cpp
    DirectionalFiniteDiffCalculator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScaledModelEvaluatorUnitTests
  SOURCES
    ScaledModelEvaluator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS thyra_nonlinear_test_models
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
