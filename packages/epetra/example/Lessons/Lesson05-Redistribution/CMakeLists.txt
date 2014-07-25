INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson05_redistribution
  SOURCES lesson05_redistribution.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
