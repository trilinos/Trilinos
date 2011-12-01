INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  verySimple
  SOURCES cxx_main.cpp
  ARGS 100
  COMM serial
#  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
