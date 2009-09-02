INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Vector_test  
  SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp
  ARGS -v
  COMM mpi serial
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
