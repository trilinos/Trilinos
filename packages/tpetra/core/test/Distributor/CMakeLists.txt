
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Distributor_UnitTests
  SOURCES
    Distributor_UnitTests.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense with > 1 MPI processes.
# See Github Issue #488.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Distributor_CreateFromSendsAndRecvs
  SOURCES
    createfromsendsandrecvs.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-8
  STANDARD_PASS_OUTPUT
  )

# See Github Issue #1454.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue1454
  SOURCES
    Issue1454.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

# See Github Issue #1752.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue1752
  SOURCES
    Issue1752.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )
