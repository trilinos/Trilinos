
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Distributor_UnitTests
  SOURCES
    Distributor_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense with > 1 MPI processes.
# See Github Issue #488.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Distributor_CreateFromSendsAndRecvs
  SOURCES
    createfromsendsandrecvs
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-8
  STANDARD_PASS_OUTPUT
  )
