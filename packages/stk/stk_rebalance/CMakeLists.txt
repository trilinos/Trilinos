
ADD_SUBDIRECTORY(stk_rebalance)
ADD_SUBDIRECTORY(stk_rebalance_utils)

PACKAGE_ADD_TEST_DIRECTORIES(unit_tests)

# The use cases currently depend on boost::program_options through 
# lib stk_util_use_cases.
IF ( TPL_Enable_BoostLib )
  PACKAGE_ADD_TEST_DIRECTORIES(use_cases)
ENDIF()

