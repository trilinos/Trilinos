INCLUDE(SetDefaultAndFromEnv)

SET(CTEST_NIGHTLY_START_TIME "04:00:00 UTC")
# NOTE: Above only is used by centralized VCS like CVS and SVN and does
# nothing for git and other Distributed VCS.  However, it needs to be set here
# to get around a defect in ctest_start() when passing in the APPEND argument
# which is used in ctest -S scripts to run tests on a different node from
# where the build is done (like in many of the ATDM Trilinos builds). For
# details, see https://gitlab.kitware.com/cmake/cmake/issues/20471.

# Set actual CTest/CDash settings

IF (NOT DEFINED CTEST_DROP_METHOD)
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_METHOD "http")
ENDIF()

IF (CTEST_DROP_METHOD STREQUAL "http" OR CTEST_DROP_METHOD STREQUAL "https")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE "testing.sandia.gov")
  SET_DEFAULT_AND_FROM_ENV(CTEST_PROJECT_NAME "Trilinos")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION "/cdash/submit.php?project=Trilinos")
  SET_DEFAULT_AND_FROM_ENV(CTEST_TRIGGER_SITE "")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_CDASH TRUE)
  # Secondary submit to development CDash site
  SET_DEFAULT_AND_FROM_ENV(TRIBITS_2ND_CTEST_DROP_SITE
    "testing-dev.sandia.gov")
  SET_DEFAULT_AND_FROM_ENV(TRIBITS_2ND_CTEST_DROP_LOCATION
    "/cdash/submit.php?project=Trilinos")
ENDIF()
