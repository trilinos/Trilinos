INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  show_stack
  COMM serial mpi
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )
