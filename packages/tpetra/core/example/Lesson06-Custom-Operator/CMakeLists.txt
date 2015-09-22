
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  lesson06_custom_operator
  SOURCES lesson06_custom_operator.cpp
  ARGS 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
