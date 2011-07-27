INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TpetraThyraWrappersUnitTests
  SOURCES
    TpetraThyraWrappers_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS "--no-run-linear-op-tester"
  )
  # Above is a defect that needs to be fixed in Tpetra ...

PACKAGE_ADD_TEST(
  TpetraThyraWrappersUnitTests
  POSTFIX_AND_ARGS_0 serial
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Simple2DTpetraModelEvaluatorUnitTests
  NUM_MPI_PROCS 1
  SOURCES
    Simple2DTpetraModelEvaluator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )
