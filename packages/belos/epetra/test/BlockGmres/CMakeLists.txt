
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bl_gmres_diag
  SOURCES test_bl_gmres_diag.cpp 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

