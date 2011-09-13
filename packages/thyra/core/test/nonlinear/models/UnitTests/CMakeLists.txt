INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Simple2DModelEvaluatorUnitTests
  SOURCES
    Simple2DModelEvaluator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS thyra_nonlinear_test_models
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE(
  DiagonalQuadraticResponseOnlyModelEvaluatorUnitTests
  SOURCES
    DiagonalQuadraticResponseOnlyModelEvaluator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS thyra_nonlinear_test_models
  COMM serial mpi
  )

PACKAGE_ADD_TEST(
  DiagonalQuadraticResponseOnlyModelEvaluatorUnitTests
  POSTFIX_AND_ARGS_0 ld1 --local-dim=1 --show-test-details=ALL
  POSTFIX_AND_ARGS_1 ld2 --local-dim=2 --show-test-details=ALL
  POSTFIX_AND_ARGS_2 ld3 --local-dim=3 --show-test-details=ALL
  POSTFIX_AND_ARGS_3 ld4_ts100 --local-dim=4 --show-test-details=ALL --tol-scale=100.0
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )
