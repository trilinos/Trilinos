
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_all_serial
  SOURCES Test_Main.cpp Test_Shards_Array.cpp Test_Shards_CellTopology.cpp
  COMM serial
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_all_mpi
  SOURCES Test_Main.cpp Test_Shards_Array.cpp Test_Shards_CellTopology.cpp
  COMM mpi
  STANDARD_PASS_OUTPUT
  )
