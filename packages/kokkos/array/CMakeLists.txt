INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)

