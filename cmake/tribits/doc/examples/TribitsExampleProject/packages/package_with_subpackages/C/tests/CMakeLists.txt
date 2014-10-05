
TRIBITS_ADD_TEST( c_util
  NOEXEPREFIX  NOEXESUFFIX
  DIRECTORY ../
  NAME test_of_c_util
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "Called c_util: A no_deps no_deps no_deps"
  )



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

ASSERT_DEFINED(PackageWithSubpackagesSubpackageB_ENABLE_MixedLanguage)
IF (PackageWithSubpackagesSubpackageB_ENABLE_MixedLanguage)

  TRIBITS_ADD_LIBRARY( c_b_mixed_lang
    TESTONLY
    SOURCES  c_b_mixed_lang.cpp  
    HEADERS  c_b_mixed_lang.hpp
    DEPLIBS  pws_c  b_mixed_lang
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST( test_of_c_b_mixed_lang
    SOURCES  test_of_c_b_mixed_lang.cpp
    TESTONLYLIBS  c_b_mixed_lang
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION
      "B A no_deps no_deps B Mixed Language"
    )

ENDIF()
