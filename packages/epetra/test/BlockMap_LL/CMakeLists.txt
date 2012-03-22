
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockMap_test_LL
  SOURCES cxx_main.cpp checkmap.cpp
  ARGS -v
  COMM mpi serial
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
