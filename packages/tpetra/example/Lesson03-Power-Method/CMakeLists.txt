INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson03_power_method
  SOURCES lesson03_power_method.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
