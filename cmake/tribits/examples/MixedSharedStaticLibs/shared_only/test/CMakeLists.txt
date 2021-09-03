tribits_add_executable_and_test( test
  SOURCES SharedOnly_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "SharedOnly_test returns shared"
  )
