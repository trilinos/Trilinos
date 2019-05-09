
ADD_SUBDIRECTORY(EpetraOperatorWrapper)
ADD_SUBDIRECTORY(UnitTests)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_epetra_adapters
  SOURCES test_epetra_adapters.cpp
  ARGS --verbose
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestThyraDebugHang
  SOURCES
    TestThyraDebugHang.cpp
    ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )
