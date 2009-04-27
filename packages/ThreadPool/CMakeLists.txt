
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(ThreadPool)

#
# B) Set up package-specific options
#

#Pthread is a required dependency so this conditional should always be true
#hence the assert after it. 
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
  PACKAGE_ADD_TEST_DIRECTORIES(test test/hpccg)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
