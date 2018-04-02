TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_tpetra_initializes_mpi
  SOURCES
    initialize_tpetra_inits_mpi
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Core_initialize_where_tpetra_initializes_kokkos
  SOURCES
    initialize_tpetra_inits_kokkos
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
  Core_initialize_where_tpetra_initializes_mpi_and_user_initializes_kokkos
  SOURCES
    initialize_tpetra_inits_mpi_user_inits_kokkos
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
