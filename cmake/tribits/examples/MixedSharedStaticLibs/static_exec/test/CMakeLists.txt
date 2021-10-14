tribits_add_executable_and_test( test
  SOURCES StaticOnly_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "StaticOnly_test returns static"
  )
