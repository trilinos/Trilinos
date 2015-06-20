

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  testLapackWrappers
  SOURCES testLapackWrappers.cpp
  ARGS
    "--n=1"
    "--n=2"
    "--n=3"
    "--n=4"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
