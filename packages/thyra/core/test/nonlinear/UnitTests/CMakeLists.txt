INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultNonlinearSolverBuilder
  SOURCES
    DefaultNonlinearSolverBuilder_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultFiniteDifferenceModelEvaluator
  SOURCES
    DefaultFiniteDifferenceModelEvaluator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
