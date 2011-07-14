INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraLinearOpUnitTests
  SOURCES
    EpetraLinearOp_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraLinearOpParallelUnitTests
  SOURCES
    EpetraLinearOpParallel_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 3
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraThyraWrappersUnitTests
  SOURCES
    EpetraThyraWrappers_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
