

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BiCGStab_diag
  SOURCES test_bicgstab_diag.cpp 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

# Dependency on Triutils fixes Bug 6449.
