tribits_add_executable(test  SOURCES  Addon1_test.cpp)

tribits_add_advanced_test( test
  OVERALL_NUM_MPI_PROCS  1
  TEST_0  EXEC  test
    PASS_REGULAR_EXPRESSION_ALL 
      "Addon1 label is: Addon1"
      "Addon1 deps are: headeronlytpl"
  )
