
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultPlatform_test
  SOURCES
    DefaultPlatform_test.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Platform_UnitTests
  SOURCES
    Platform_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# The "issue 434" tests exercise both cases of Github Issue #434:
#
# https://github.com/trilinos/Trilinos/issues/434
#
# The first case is that Kokkos and MPI are already initialized before
# Tpetra::initialize is called.  The second case is that neither
# Kokkos nor MPI have been initialized before Tpetra::initialize is
# called.
#
# Neither of these tests may use the standard Teuchos unit test
# executable, because they both must control and test Kokkos and MPI
# initialization and finalization directly.
#
# It suffices to use only 1 MPI process for these tests.  That also
# simplifies output, since according to the MPI standard, only Process
# 0 in MPI_COMM_WORLD may write to stdout and stderr.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  issue_434_already_initialized
  SOURCES issue_434_already_initialized.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  issue_434_not_yet_initialized
  SOURCES issue_434_not_yet_initialized.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
