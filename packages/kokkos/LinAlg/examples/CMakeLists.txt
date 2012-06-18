INCLUDE(TribitsAddExecutableAndTest)   
INCLUDE(TribitsLibraryMacros)

TRIBITS_ADD_EXECUTABLE_AND_TEST(  EmptySparseKernelDriver  
  SOURCES EmptySparseKernelDriver.cpp  
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 

# mfh 18 June 2012: This is failing to build on some machines.  See e.g., 
#
# http://testing.sandia.gov/cdash/viewBuildError.php?buildid=542885
#
# That probably has to do with the recent change to the Kokkos
# interface to sparse matrix kernels, which probably hasn't yet been
# propagated to this test.  Before reenabling the test, make sure that
# it uses the new sparse kernels interface.
#TRIBITS_ADD_TEST_DIRECTORIES(PDP10)
