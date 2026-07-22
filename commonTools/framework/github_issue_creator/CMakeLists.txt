tribits_add_advanced_test( create_trilinos_github_test_failure_issue_unit_tests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND ${Python3_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/create_trilinos_github_test_failure_issue_unit_tests.py -v
    PASS_REGULAR_EXPRESSION "OK"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  )
