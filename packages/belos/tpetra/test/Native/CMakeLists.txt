# A diagonal matrix with varying positive diagonal entries is a good
# sanity test for all Krylov solvers.  It also is a sufficient test
# for CG and related algorithms.

TRIBITS_ADD_EXECUTABLE(
  Tpetra_Native_solvers_diagonal_matrix
  SOURCES diagonal_matrix.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_CG_diagonal_matrix
  ARGS "--solver=\"TPETRA CG\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_CG_PIPELINE_diagonal_matrix
  ARGS "--solver=\"TPETRA CG PIPELINE\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_CG_SINGLE_REDUCE_diagonal_matrix
  ARGS "--solver=\"TPETRA CG SINGLE REDUCE\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_GMRES_diagonal_matrix
  ARGS "--solver=\"TPETRA GMRES\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_GMRES_PIPELINE_diagonal_matrix
  ARGS "--solver=\"TPETRA GMRES PIPELINE\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_GMRES_SINGLE_REDUCE_diagonal_matrix
  ARGS "--solver=\"TPETRA GMRES SINGLE REDUCE\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_diagonal_matrix
  NAME Tpetra_Native_GMRES_S_STEP_diagonal_matrix
  ARGS "--solver=\"TPETRA GMRES S-STEP\""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# For GMRES and other algorithms, we also need to test their ability
# to solve nonsymmetric linear systems.

TRIBITS_ADD_EXECUTABLE(
  Tpetra_Native_solvers_nonsymm_matrix
  SOURCES nonsymm_matrix.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_nonsymm_matrix
  NAME Tpetra_Native_GMRES_nonsymm_matrix
  ARGS "--solver=\"TPETRA GMRES\" --restartLength=10"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_nonsymm_matrix
  NAME Tpetra_Native_GMRES_PIPELINE_nonsymm_matrix
  ARGS "--solver=\"TPETRA GMRES PIPELINE\" --restartLength=10"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_nonsymm_matrix
  NAME Tpetra_Native_GMRES_SINGLE_REDUCE_nonsymm_matrix
  ARGS "--solver=\"TPETRA GMRES SINGLE REDUCE\" --restartLength=10"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  Tpetra_Native_solvers_nonsymm_matrix
  NAME Tpetra_Native_GMRES_S_STEP_nonsymm_matrix
  ARGS "--solver=\"TPETRA GMRES S-STEP\" --restartLength=10"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
