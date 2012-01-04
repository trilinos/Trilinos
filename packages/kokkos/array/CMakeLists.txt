INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(usecases)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)

