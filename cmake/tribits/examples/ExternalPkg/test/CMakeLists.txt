TRIBITS_ADD_EXECUTABLE_AND_TEST( test
  SOURCES ExternalPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "ExteranlPkg deps are: no_deps"
  )
