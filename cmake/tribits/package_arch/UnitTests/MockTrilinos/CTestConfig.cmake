INCLUDE(SetDefaultAndFromEnv)

# Should never be used to submit anywhere!  This is just for the unit test
# mock project!

# Must match what is in CDash project 'Trilinos'
SET(CTEST_NIGHTLY_START_TIME "00:00:00 UTC")

# Set actual CTest/CDash settings

IF (NOT DEFINED CTEST_DROP_METHOD)
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_METHOD "http")
ENDIF()

IF (CTEST_DROP_METHOD STREQUAL "http")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE "dummy.com")
  SET_DEFAULT_AND_FROM_ENV(CTEST_PROJECT_NAME "MockProjectName")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION "/cdash/submit.php?project=MockProjectName")
  SET_DEFAULT_AND_FROM_ENV(CTEST_TRIGGER_SITE "")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_CDASH TRUE)
ENDIF()
