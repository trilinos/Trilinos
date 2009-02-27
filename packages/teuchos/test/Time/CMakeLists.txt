
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Time_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )
