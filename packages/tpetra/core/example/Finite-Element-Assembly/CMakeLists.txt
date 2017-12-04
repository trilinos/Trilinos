
# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FiniteElementAssembly_type1
  SOURCES fem_assembly_type1.cpp
  ARGS 
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FiniteElementAssembly_type2
  SOURCES fem_assembly_type2.cpp
  ARGS 
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FiniteElementAssembly_type3
  SOURCES fem_assembly_type3.cpp
  ARGS 
  COMM  mpi
  STANDARD_PASS_OUTPUT
)

