
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  test_all
  SOURCES Test_Main.cpp Test_Shards_Array.cpp Test_Shards_CellTopology.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
