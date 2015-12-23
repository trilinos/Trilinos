TRIBITS_ADD_EXECUTABLE_AND_TEST( test
  SOURCES StaticOnly_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "StaticOnly_test returns static"
  )
