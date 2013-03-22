INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_test_LL
  SOURCES cxx_main.cpp 
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  PASS_REGULAR_EXPRESSION "tests PASSED"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_memorytest_LL
  SOURCES memorytest_main.cpp 
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1 #1-20
  PASS_REGULAR_EXPRESSION "tests PASSED"
  )
