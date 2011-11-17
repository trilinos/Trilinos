INCLUDE(TribitsPackageMacros)

ADD_SUBDIRECTORY(stk_util)

IF ( TPL_ENABLE_BoostLib )
  TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)
ENDIF()

