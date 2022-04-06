
# This test does not depend on MPI and only needs one process, but it
# still works if MPI is enabled, and we want to test it regardless.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue510
  SOURCES
    Issue510.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--numTimesToCreateNode=1"
  STANDARD_PASS_OUTPUT
  )

# This test does not depend on MPI and only needs one process, but it
# still works if MPI is enabled, and we want to test it regardless.
TRIBITS_ADD_TEST(
  Issue510
  NAME Issue510_Create2x
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--numTimesToCreateNode=2"
  STANDARD_PASS_OUTPUT
  )
