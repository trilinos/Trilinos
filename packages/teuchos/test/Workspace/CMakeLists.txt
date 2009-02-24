
# Required modules
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Workspace_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  NUM_MPI_PROCS 1
  COMM serial mpi
  )
