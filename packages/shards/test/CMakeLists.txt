
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_all_serial
  SOURCES Test_Main.cpp Test_Shards_Array.cpp Test_Shards_CellTopology.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
