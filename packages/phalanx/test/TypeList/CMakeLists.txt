
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  type_list
  SOURCES TypeListTest.cpp
  TESTONLYLIBS phalanx_unit_test_main
  COMM serial mpi
  NUM_MPI_PROCS 1
)
