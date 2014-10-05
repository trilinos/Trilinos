#
# Create a test executable b_test that only depends on package's man libs
#

TRIBITS_ADD_EXECUTABLE( b_test SOURCES b_test.cpp )

#
# Create a test-only b_mixed_lang lib which uses the optional test dependent
# MixedLang package
#

ADD_SUBDIRECTORY(testlib)

#
# Test b_test
#

SET(EXPECTED_B_DEPS)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_PackageWithSubpackagesSubpackageA)
IF(${PACKAGE_NAME}_ENABLE_PackageWithSubpackagesSubpackageA)
  SET(EXPECTED_B_DEPS "${EXPECTED_B_DEPS} A no_deps")
ENDIF()
SET(EXPECTED_B_DEPS "${EXPECTED_B_DEPS} no_deps")

TRIBITS_ADD_ADVANCED_TEST( test_of_b
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 EXEC b_test
    PASS_REGULAR_EXPRESSION_ALL 
      "B label is: B"
      "B deps are:${EXPECTED_B_DEPS}"
  )

#
# Define an test involving the optional dependent MixedLang package
# 

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MixedLanguage)
IF (${PACKAGE_NAME}_ENABLE_MixedLanguage)

  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

  TRIBITS_ADD_EXECUTABLE_AND_TEST( test_of_b_mixed_lang
    SOURCES b_mixed_lang_test.cpp
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION 
      "B Mixed Language"
    TESTONLYLIBS b_mixed_lang
    ${TAEAT_EXTRALIB_ARGS} # Used for testing and to demonstrate errors
    )

ENDIF()
