#
# Create a test executable b_test that only depends on package's man libs
#

tribits_add_executable( b_test SOURCES b_test.cpp
  ADDED_EXE_TARGET_NAME_OUT b_test_TARGET_NAME )

tribits_verbose_print_var(b_test_TARGET_NAME)

# Example of setting a property on an executable
set_target_properties(${b_test_TARGET_NAME}
  PROPERTIES LINKER_LANGUAGE CXX ) # Don't really need to set this!

#
# Create a test-only b_mixed_lang lib which uses the optional test dependent
# MixedLang package
#

add_subdirectory(testlib)

#
# Test b_test
#

tribits_add_advanced_test( test_of_b
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC b_test
    PASS_REGULAR_EXPRESSION_ALL 
      "B label is: B"
      "B deps are: ${EXPECTED_B_DEPS}"
  ADDED_TEST_NAME_OUT  test_of_b_TEST_NAME
  )

tribits_verbose_print_var(test_of_b_TEST_NAME)

# Example of changing CTest property on advanced test
if (test_of_b_TEST_NAME)
  set_tests_properties( ${test_of_b_TEST_NAME}
    PROPERTIES  TIMEOUT  12.0 )
endif()

#
# Define a test involving the optional dependent MixedLang package
# 

assert_defined(${PACKAGE_NAME}_ENABLE_MixedLang)
if (${PACKAGE_NAME}_ENABLE_MixedLang)

  include(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

  tribits_add_executable_and_test( test_of_b_mixed_lang
    SOURCES  b_mixed_lang_test.cpp
    NUM_MPI_PROCS  1
    PASS_REGULAR_EXPRESSION  "B Mixed Language Hello World"
    TESTONLYLIBS  b_test_utils  b_mixed_lang
    ${TAEAT_EXTRALIB_ARGS} # Used for testing and to demonstrate errors
    )

endif()
