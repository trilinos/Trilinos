INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Oski_test 
  SOURCES cxx_main.cpp
  COMM mpi serial
  NUM_MPI_PROCS 1
  )
