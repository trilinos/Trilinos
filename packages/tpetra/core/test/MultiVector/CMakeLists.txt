TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug5474
  SOURCES
    Bug5474
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test requires exactly 2 MPI processes.  It could be rewritten
# so that it will run with more, but then the test would just waste
# the remaining processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MV_subViewSomeZeroRows
  SOURCES
    MV_subViewSomeZeroRows
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not require or use MPI,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MV_Kokkos_DualView_subview
  SOURCES
  MV_Kokkos_DualView_subview
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test comes from zoltan2/test/temp/multivectorTest.cpp.  It
# mimics the recursive bisection algorithm in Zoltan2, and times data
# migration.
#
# The test itself does not refer to MPI.  Hoewver, it only makes sense
# to run this test in an MPI build.  Otherwise, no data migration is
# necessary.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  rcb
  SOURCES rcb.cpp ${MV_RCB_OBJECTS_CUDA}
  COMM mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

# This test relates to Github Issue #358.
# It does not need multiple MPI processes,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_get2dView
  SOURCES
    MultiVector_get2dView
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test relates to Github Issue #567.
# It does not need multiple MPI processes,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_TwoArgRandomize
  SOURCES
    MultiVector_TwoArgRandomize
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not need MPI,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue364
  SOURCES
    Issue364
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not need MPI,
# but it should still work if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue46
  SOURCES
    Issue46
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_MicroBenchmark
  SOURCES
    MultiVector_MicroBenchmark
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )
