
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bl_gmres_multidiag
  SOURCES test_bl_gmres_multidiag.cpp BelosLinearMultiShiftProblem.hpp 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

