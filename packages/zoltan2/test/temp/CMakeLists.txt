INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Zoltan2_paramtest
  SOURCES paramtest.cpp
  NUM_MPI_PROCS 1
  COMM mpi serial
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Zoltan2_ostreamtest
  SOURCES ostreamtest.cpp
  NUM_MPI_PROCS 1
  COMM mpi serial
  )
