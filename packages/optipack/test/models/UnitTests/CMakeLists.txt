INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DiagonalQuadraticResponseOnlyModelEvaluatorUnitTests
  SOURCES
    DiagonalQuadraticResponseOnlyModelEvaluator_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS optipack_test_models
  ARGS
    "--local-dim=1 --show-test-details=ALL"
    "--local-dim=2 --show-test-details=ALL"
    "--local-dim=3 --show-test-details=ALL"
    "--local-dim=4 --show-test-details=ALL --tol-scale=100.0 "
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )
