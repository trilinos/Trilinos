
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  XML_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )
