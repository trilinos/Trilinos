TRIBITS_ADD_ADVANCED_TEST( ats2_custom_builds_unit_tests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND env
    ARGS SHUNIT_COLOR=none ${CMAKE_CURRENT_SOURCE_DIR}/custom_builds_unit_tests.sh
    PASS_REGULAR_EXPRESSION_ALL
      "ERROR: A supported compiler was not selected for 'ats2' env - anything-cudas-after"
      "Ran 3 tests"
      "[^]OK[$]"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  )
