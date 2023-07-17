

#
# A) Define the package
#

TRIBITS_PACKAGE(Intrepid2)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_INTREPID2_DEBUG
  "Enable a host of runtime debug checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

IF(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK_OVERRIDE)
  ## override the inf flag when it is needed e.g., refactoring dyn rank view
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK
    HAVE_INTREPID2_DEBUG_INF_CHECK
    "Enable a host of runtime debug checking for infs and nans."
    ${${PROJECT_NAME}_ENABLE_DEBUG})
ELSE()
  ## unless override is enabled, inf check is turned off by default
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK
    HAVE_INTREPID2_DEBUG_INF_CHECK
    "Enable a host of runtime debug checking for infs and nans."
    OFF)
ENDIF()   

#
# C) Add the libraries, tests, and examples
#
ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(unit-test)

# performance tests 
TRIBITS_ADD_EXAMPLE_DIRECTORIES(perf-test)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
