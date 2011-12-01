INCLUDE(TribitsAddExecutableAndTest)   
INCLUDE(TribitsLibraryMacros)

TRIBITS_ADD_EXECUTABLE_AND_TEST(  DummySparseKernelDriver  
  SOURCES DummySparseKernelDriver.cpp  
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
TRIBITS_ADD_TEST_DIRECTORIES(PDP10)
