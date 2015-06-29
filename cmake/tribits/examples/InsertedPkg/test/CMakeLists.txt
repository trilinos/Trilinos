TRIBITS_ADD_EXECUTABLE_AND_TEST( test
  SOURCES InsertedPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "InsertedPkg deps are: no_deps"
  )
