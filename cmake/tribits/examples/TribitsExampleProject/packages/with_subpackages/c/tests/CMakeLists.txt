tribits_add_test( c_util
  NOEXEPREFIX  NOEXESUFFIX
  DIRECTORY ../
  NAME test_of_c_util
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "Called c_util: B ${EXPECTED_B_DEPS} ${EXPECTED_A_AND_DEPS}"
  ADDED_TESTS_NAMES_OUT  c_util_TEST_NAME
  )

tribits_verbose_print_var(c_util_TEST_NAME)

# Example of how to set properties on added test
if (c_util_TEST_NAME)
  set_tests_properties( ${c_util_TEST_NAME}  PROPERTIES  TIMEOUT  15.0 )
  # NOTE: You don't have to set the TIMEOUT separately.  You can pass it into
  # tribits_add_tet().
endif()
  
tribits_add_executable( c_test SOURCES c_test.cpp )

tribits_add_advanced_test( test_of_c
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC c_test
    PASS_REGULAR_EXPRESSION_ALL 
      "A label is: A"
      "B label is: B"
      "C label is: C"
      "C deps are: B A"
  )


#
# Define a test-only lib and two tests that use it.
#

assert_defined(WithSubpackagesB_ENABLE_MixedLang)
if (WithSubpackagesB_ENABLE_MixedLang AND
   (WithSubpackagesB_ENABLE_TESTS OR WithSubpackages_ENABLE_TESTS)
  )

  # Use a TESTONLY lib b_mixed_lang from WithSubpackagesB to define a TESTONLY
  # here!  It is not recommended to share test libraries between package but
  # this example shows it can be done.

  tribits_add_library( c_b_mixed_lang
    TESTONLY
    SOURCES  c_b_mixed_lang.cpp  
    HEADERS  c_b_mixed_lang.hpp
    DEPLIBS  pws_c  b_mixed_lang
    )
  
  tribits_add_executable_and_test( test_of_c_b_mixed_lang
    SOURCES  test_of_c_b_mixed_lang.cpp
    TESTONLYLIBS  c_b_mixed_lang
    NUM_MPI_PROCS  1
    PASS_REGULAR_EXPRESSION
      "B ${EXPECTED_B_DEPS} ${EXPECTED_A_AND_DEPS} B Mixed Language"
    ADDED_EXE_TARGET_NAME_OUT  test_of_c_b_mixed_lang_TARGET_NAME
    ADDED_TESTS_NAMES_OUT  test_of_c_b_mixed_lang_TEST_NAME
    )

  tribits_verbose_print_var(test_of_c_b_mixed_lang_TARGET_NAME)
  tribits_verbose_print_var(test_of_c_b_mixed_lang_TEST_NAME)

  # Example on how to set exe target properties.
  set_target_properties( ${test_of_c_b_mixed_lang_TARGET_NAME}
    PROPERTIES  LINKER_LANGUAGE  CXX )
  # NOTE: You don't have to set the LINKER_LANGUAGE

  if (test_of_c_b_mixed_lang_TEST_NAME)
    set_tests_properties( ${test_of_c_b_mixed_lang_TEST_NAME}
      PROPERTIES  TIMEOUT  3.0 )
    # NOTE: You don't have to set the TIMEOUT.
  endif()

endif()
