
# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMAssembly_InsertGlobalIndices_DP
  SOURCES fem_assembly_InsertGlobalIndices_DP.cpp
  ARGS
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMAssembly_LocalElementLoop_DP
  SOURCES fem_assembly_LocalElementLoop_DP.cpp
  ARGS
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMAssembly_TotalElementLoop_DP
  SOURCES fem_assembly_TotalElementLoop_DP.cpp
  ARGS
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMAssembly_TotalElementLoop_SP
  SOURCES fem_assembly_TotalElementLoop_SP.cpp
  ARGS
  COMM  mpi
  STANDARD_PASS_OUTPUT
)


