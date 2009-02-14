
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(ThreadPool)

#
# B) Set up package-specific options
#

CHECK_INCLUDE_FILE_CXX(pthread.h HAVE_PTHREAD)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  PRINT_VAR(HAVE_PTHREAD)
ENDIF()
ASSERT_DEFINED(HAVE_PTHREAD)
IF(HAVE_PTHREAD)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "ADDING THREADPOOL TESTS")
  ENDIF()
  PACKAGE_ADD_TEST_DIRECTORIES(test)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
