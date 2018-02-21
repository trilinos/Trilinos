
# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockCrsPerfTest
  SOURCES driver.cpp
  ARGS 
  COMM  mpi
  STANDARD_PASS_OUTPUT
)
