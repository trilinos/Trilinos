
INCLUDE(TribitsAddExecutableAndTest)
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalarTraits_test  
  SOURCES ScalarTraits_test.cpp 
  COMM serial mpi  
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
