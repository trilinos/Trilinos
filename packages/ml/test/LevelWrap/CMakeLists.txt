
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MLP_LW
  SOURCES LevelTest.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
