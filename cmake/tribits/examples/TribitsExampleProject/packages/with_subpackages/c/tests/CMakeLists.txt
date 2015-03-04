IF(WithSubpackagesB_ENABLE_ExternalPkg)
  SET(EXPECTED_EXTERNALPKG_DEPS "ExternalPkg no_deps ")
ENDIF()

TRIBITS_ADD_TEST( c_util
  NOEXEPREFIX  NOEXESUFFIX
  DIRECTORY ../
  NAME test_of_c_util
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "Called c_util: A no_deps ${EXPECTED_EXTERNALPKG_DEPS}no_deps no_deps"
  ADDED_TESTS_NAMES_OUT  c_util_TEST_NAME
  )

TRIBITS_VERBOSE_PRINT_VAR(c_util_TEST_NAME)

# Example of how to set properties on added test
IF (c_util_TEST_NAME)
  SET_TESTS_PROPERTIES( ${c_util_TEST_NAME}  PROPERTIES  TIMEOUT  15.0 )
  # NOTE: You don't have to set the TIMEOUT separately.  You can pass it into
  # TRIBITS_ADD_TET().
ENDIF()
  
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


#
# Define a test-only lib and two tests that use it.
#

ASSERT_DEFINED(WithSubpackagesB_ENABLE_MixedLang)
IF (WithSubpackagesB_ENABLE_MixedLang AND
   (WithSubpackagesB_ENABLE_TESTS OR WithSubpackages_ENABLE_TESTS)
  )

  # Use a TESTONLY lib b_mixed_lang from WithSubpackagesB to define a TESTONLY
  # here!  It is not recommended to share test libraries between package but
  # this example shows it can be done.

  TRIBITS_ADD_LIBRARY( c_b_mixed_lang
    TESTONLY
    SOURCES  c_b_mixed_lang.cpp  
    HEADERS  c_b_mixed_lang.hpp
    DEPLIBS  pws_c  b_mixed_lang
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST( test_of_c_b_mixed_lang
    SOURCES  test_of_c_b_mixed_lang.cpp
    TESTONLYLIBS  c_b_mixed_lang
    NUM_MPI_PROCS  1
    PASS_REGULAR_EXPRESSION
      "B A no_deps no_deps B Mixed Language"
    ADDED_EXE_TARGET_NAME_OUT  test_of_c_b_mixed_lang_TARGET_NAME
    ADDED_TESTS_NAMES_OUT  test_of_c_b_mixed_lang_TEST_NAME
    )

  TRIBITS_VERBOSE_PRINT_VAR(test_of_c_b_mixed_lang_TARGET_NAME)
  TRIBITS_VERBOSE_PRINT_VAR(test_of_c_b_mixed_lang_TEST_NAME)

  # Example on how to set exe target properties.
  SET_TARGET_PROPERTIES( ${test_of_c_b_mixed_lang_TARGET_NAME}
    PROPERTIES  LINKER_LANGUAGE  CXX )
  # NOTE: You don't have to set the LINKER_LANGUAGE

  IF (test_of_c_b_mixed_lang_TEST_NAME)
    SET_TESTS_PROPERTIES( ${test_of_c_b_mixed_lang_TEST_NAME}
      PROPERTIES  TIMEOUT  3.0 )
    # NOTE: You don't have to set the TIMEOUT.
  ENDIF()

ENDIF()
