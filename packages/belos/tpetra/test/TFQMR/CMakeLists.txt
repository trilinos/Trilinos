
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_tfqmr_diag
  SOURCES test_tfqmr_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

