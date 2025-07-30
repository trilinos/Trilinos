
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BiCGStab_diag
  SOURCES test_bicgstab_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

