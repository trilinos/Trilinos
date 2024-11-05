tribits_add_executable_and_test( test
  SOURCES InsertedPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "InsertedPkg deps are: ${EXPECTED_SIMPLECXX_AND_DEPS}"
  )
