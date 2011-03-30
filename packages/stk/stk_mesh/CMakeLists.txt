
ADD_SUBDIRECTORY(stk_mesh)


# The use cases and unit test currently depend on boost::program_options through 
# lib stk_util_use_cases.
IF ( TPL_Enable_BoostLib AND ${PACKAGE_NAME}_ENABLE_TESTS)
  PACKAGE_ADD_TEST_DIRECTORIES(unit_tests)
  PACKAGE_ADD_TEST_DIRECTORIES(use_cases)
ENDIF()


