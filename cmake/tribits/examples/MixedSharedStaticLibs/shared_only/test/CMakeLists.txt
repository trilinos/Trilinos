TRIBITS_ADD_EXECUTABLE_AND_TEST( test
  SOURCES SharedOnly_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "SharedOnly_test returns shared"
  )
