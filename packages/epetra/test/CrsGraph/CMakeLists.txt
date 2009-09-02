INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_test 
  SOURCES cxx_main.cpp 
  KEYWORDS unit
  ARGS -v
  COMM mpi serial
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
