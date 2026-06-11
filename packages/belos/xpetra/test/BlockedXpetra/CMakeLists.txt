
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Xpetra_BlockedCrsMatrixGMRES_test
  SOURCES test_xpetra_blocked_crs_matrix.cpp
  ARGS "--verbose"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
