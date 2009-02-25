INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultNonlinearSolverBuilder
  SOURCES
    DefaultNonlinearSolverBuilder_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultFiniteDifferenceModelEvaluator
  SOURCES
    DefaultFiniteDifferenceModelEvaluator_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
