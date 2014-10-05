INCLUDE(TribitsAddExecutableAndTest)

# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson01_mpi_only_through_Epetra
  SOURCES lesson01_mpi_only_through_Epetra.cpp 
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

# This example _requires_ MPI.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson01_mpi_on_its_own
  SOURCES lesson01_mpi_on_its_own.cpp 
  ARGS 
  COMM mpi
  STANDARD_PASS_OUTPUT
)

# This example does not _use_ MPI,
# but it works whether or not MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson01_no_mpi
  SOURCES lesson01_no_mpi.cpp 
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
