TRIBITS_ADD_EXECUTABLE(test  SOURCES  Addon1_test.cpp)

TRIBITS_ADD_ADVANCED_TEST( test
  OVERALL_NUM_MPI_PROCS  1
  TEST_0  EXEC  test
    PASS_REGULAR_EXPRESSION_ALL 
      "Addon1 label is: Addon1"
      "Addon1 deps are: no_deps"
  )
