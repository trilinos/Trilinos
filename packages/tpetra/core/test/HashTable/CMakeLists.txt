TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FixedHashTableTest
  SOURCES
    FixedHashTableTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# We don't use the Teuchos standard unit test executable here, because
# we want control over Kokkos initialization and finalization.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  computeOffsetsFromCounts
  SOURCES computeOffsetsFromCounts
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# We don't use the Teuchos standard unit test executable here, because
# we want control over Kokkos initialization and finalization.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  computeOffsetsFromConstantCount
  SOURCES computeOffsetsFromConstantCount
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

