INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Simple2DModelEvaluatorUnitTests
  SOURCES
    Simple2DModelEvaluator_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS thyra_nonlinear_test_models
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
