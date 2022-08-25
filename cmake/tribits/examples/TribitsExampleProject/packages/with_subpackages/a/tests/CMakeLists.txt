tribits_add_executable(a_test SOURCES a_test.cpp)

tribits_add_advanced_test( test_of_a
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC a_test
    PASS_REGULAR_EXPRESSION_ALL 
      "A label is: A"
      "A deps are: ${EXPECTED_SIMPLECXX_AND_DEPS}"
      "A special value: ${WithSubpackagesA_SPECIAL_VALUE}"
  )
