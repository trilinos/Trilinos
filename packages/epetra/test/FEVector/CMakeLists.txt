INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  FEVector_test
  SOURCES cxx_main.cpp ExecuteTestProblems.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
  
