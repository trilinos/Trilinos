#
# Create a test executable b_test that only depends on package's man libs
#

TRIBITS_ADD_EXECUTABLE( b_test SOURCES b_test.cpp
  ADDED_EXE_TARGET_NAME_OUT b_test_TARGET_NAME )

TRIBITS_VERBOSE_PRINT_VAR(b_test_TARGET_NAME)

# Example of setting a property on an executable
SET_TARGET_PROPERTIES(${b_test_TARGET_NAME}
  PROPERTIES LINKER_LANGUAGE CXX ) # Don't really need to set this!

#
# Create a test-only b_mixed_lang lib which uses the optional test dependent
# MixedLang package
#

ADD_SUBDIRECTORY(testlib)

#
# Test b_test
#

SET(EXPECTED_B_DEPS)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_WithSubpackagesA)
IF(${PACKAGE_NAME}_ENABLE_WithSubpackagesA)
  SET(EXPECTED_B_DEPS "${EXPECTED_B_DEPS} A no_deps")
ENDIF()
SET(EXPECTED_B_DEPS "${EXPECTED_B_DEPS} no_deps")

TRIBITS_ADD_ADVANCED_TEST( test_of_b
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC b_test
    PASS_REGULAR_EXPRESSION_ALL 
      "B label is: B"
      "B deps are:${EXPECTED_B_DEPS}"
  ADDED_TEST_NAME_OUT  test_of_b_TEST_NAME
  )

TRIBITS_VERBOSE_PRINT_VAR(test_of_b_TEST_NAME)

# Example of changing CTest property on advanced test
IF (test_of_b_TEST_NAME)
  SET_TESTS_PROPERTIES( ${test_of_b_TEST_NAME}
    PROPERTIES  TIMEOUT  3.0 )
ENDIF()

#
# Define a test involving the optional dependent MixedLang package
# 

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MixedLang)
IF (${PACKAGE_NAME}_ENABLE_MixedLang)

  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

  TRIBITS_ADD_EXECUTABLE_AND_TEST( test_of_b_mixed_lang
    SOURCES  b_mixed_lang_test.cpp
    NUM_MPI_PROCS  1
    PASS_REGULAR_EXPRESSION  "B Mixed Language"
    TESTONLYLIBS  b_mixed_lang
    ${TAEAT_EXTRALIB_ARGS} # Used for testing and to demonstrate errors
    )

ENDIF()
