TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_tpetra_initializes_mpi
  SOURCES
    initialize_tpetra_inits_mpi
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_tpetra_initializes_mpi
  SOURCES
    ScopeGuard_tpetra_inits_mpi
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_user_initializes_mpi
  SOURCES
    initialize_user_inits_mpi
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_user_initializes_mpi
  SOURCES
    ScopeGuard_user_inits_mpi
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

# This test does an MPI_Comm_split,
# so it needs at least 2 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_user_initializes_mpi_and_provides_comm
  SOURCES
    initialize_user_inits_mpi_and_provides_comm
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

# This test does an MPI_Comm_split,
# so it needs at least 2 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_user_initializes_mpi_and_provides_comm
  SOURCES
    ScopeGuard_user_inits_mpi_and_provides_comm
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_tpetra_initializes_kokkos
  SOURCES
    initialize_tpetra_inits_kokkos
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_tpetra_initializes_kokkos
  SOURCES
    ScopeGuard_tpetra_inits_kokkos
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_user_initializes_kokkos
  SOURCES
    initialize_user_inits_kokkos
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_user_initializes_kokkos
  SOURCES
    ScopeGuard_user_inits_kokkos
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_tpetra_initializes_mpi_and_user_initializes_kokkos
  SOURCES
    initialize_tpetra_inits_mpi_user_inits_kokkos
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_tpetra_initializes_mpi_and_user_initializes_kokkos
  SOURCES
    ScopeGuard_tpetra_inits_mpi_user_inits_kokkos
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_user_initializes_mpi_and_tpetra_initializes_kokkos
  SOURCES
    initialize_user_inits_mpi_tpetra_inits_kokkos
  COMM mpi
  NUM_MPI_PROCS 2  
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_ScopeGuard_where_user_initializes_mpi_and_tpetra_initializes_kokkos
  SOURCES
    ScopeGuard_user_inits_mpi_tpetra_inits_kokkos
  COMM mpi
  NUM_MPI_PROCS 2  
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
