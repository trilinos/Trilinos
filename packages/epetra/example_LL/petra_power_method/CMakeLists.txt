INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  petra_power_method
  SOURCES cxx_main.cpp
  ARGS 100
  COMM serial
#  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
  
