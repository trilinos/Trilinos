
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tpetra_pseudo_gmres_multidiag
  SOURCES test_pseudo_gmres_multidiag.cpp BelosLinearMultiShiftProblem.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
