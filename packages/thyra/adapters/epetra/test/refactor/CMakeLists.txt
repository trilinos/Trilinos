TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraThyraWrappersUnitTests
  SOURCES
    EpetraThyraWrappers_UnitTests.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  )
