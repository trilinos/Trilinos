
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Polling_test  
  SOURCES cxx_main.cpp 
  COMM mpi  
  NUM_MPI_PROCS 2
)
