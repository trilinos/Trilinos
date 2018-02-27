
# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockCrsPerfTest
  SOURCES Tpetra_TestBlockCrs.cpp
  ARGS 
  COMM  mpi serial
  STANDARD_PASS_OUTPUT
)
