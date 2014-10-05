
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(ThreadPool)

#
# B) Set up package-specific options
#

#Pthread is a required dependency so this conditional should always be true
#hence the assert after it. 
SET(HAVE_PTHREAD FALSE)
IF(TPL_ENABLE_Pthread)
  SET(HAVE_PTHREAD TRUE)
ENDIF()
IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  PRINT_VAR(HAVE_PTHREAD)
ENDIF()
ASSERT_DEFINED(HAVE_PTHREAD)

IF(${PACKAGE_NAME}_ENABLE_MPI)
  SET(HAVE_MPI TRUE)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF(HAVE_PTHREAD)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "ADDING THREADPOOL TESTS")
  ENDIF()
  TRIBITS_ADD_TEST_DIRECTORIES(test test/hpccg test/hhpccg)
ENDIF()

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
