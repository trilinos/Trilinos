
TRIBITS_ADD_EXECUTABLE_AND_TEST(  EmptySparseKernelDriver  
  SOURCES EmptySparseKernelDriver.cpp  
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
) 

# NOTE (mfh 09 Jan 2015): Examples here are likely obsolete,
# given the move to the Kokkos refactor version of Tpetra.
#TRIBITS_ADD_TEST_DIRECTORIES(PDP10)
