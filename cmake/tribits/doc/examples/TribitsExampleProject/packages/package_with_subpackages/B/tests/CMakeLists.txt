TRIBITS_ADD_EXECUTABLE( b_test SOURCES b_test.cpp )

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MixedLanguage)
IF (${PACKAGE_NAME}_ENABLE_MixedLanguage)

  TRIBITS_ADD_LIBRARY( b_mixed_lang
     HEADERS b_mixed_lang.hpp
     SOURCES b_mixed_lang.cpp
     TESTONLY
     )

  TRIBITS_ADD_EXECUTABLE_AND_TEST( test_of_b_mixed_lang
    SOURCES b_mixed_lang_test.cpp
    TESTONLYLIBS b_mixed_lang
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION 
      "B Mixed Language"
     )

ENDIF()
