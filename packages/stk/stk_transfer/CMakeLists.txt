INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(Transfer)

#this is apparently a header only library so no need to add the stk_transfer directory

TRIBITS_ADD_TEST_DIRECTORIES(use_cases)


TRIBITS_SUBPACKAGE_POSTPROCESS()
