
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug_5794_IndexBase_LL
  SOURCES cxx_main.cpp
  ARGS 10
  COMM mpi serial
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
