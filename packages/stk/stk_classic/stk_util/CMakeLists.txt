INCLUDE(TribitsPackageMacros)

ADD_SUBDIRECTORY(stk_util)

IF ( STK_ENABLE_BoostLib )
  TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)
ENDIF()

