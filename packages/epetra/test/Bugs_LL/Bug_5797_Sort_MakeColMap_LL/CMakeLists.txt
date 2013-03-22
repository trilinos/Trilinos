
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug_5797_Sort_MakeColMap_LL
  SOURCES cxx_main.cpp
  ARGS 25
  COMM mpi serial
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
