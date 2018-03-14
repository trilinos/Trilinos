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
