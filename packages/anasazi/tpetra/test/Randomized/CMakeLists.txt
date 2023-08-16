
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Randomized_Solver_Test
  SOURCES rand_test.cpp 
  ARGS 
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Randomized_Solver_Driver
  SOURCES rand_driver.cpp 
  ARGS 
    "--verbose"
  COMM serial mpi
  )

