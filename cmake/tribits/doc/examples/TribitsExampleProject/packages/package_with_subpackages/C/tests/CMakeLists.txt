TRIBITS_ADD_EXECUTABLE( c_test SOURCES c_test.cpp )

TRIBITS_ADD_ADVANCED_TEST( test_of_c
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC c_test
    PASS_REGULAR_EXPRESSION_ALL 
      "A label is: A"
      "B label is: B"
      "C label is: C"
      "C deps are: B A"
  )
