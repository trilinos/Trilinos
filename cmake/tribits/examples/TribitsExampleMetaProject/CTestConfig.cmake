include(SetDefaultAndFromEnv)

set(CTEST_NIGHTLY_START_TIME "04:00:00 UTC") # 10 PM MDT or 9 PM MST

if (NOT DEFINED CTEST_DROP_METHOD)
  set_default_and_from_env(CTEST_DROP_METHOD "https")
endif()

if (CTEST_DROP_METHOD STREQUAL "http" OR CTEST_DROP_METHOD STREQUAL "https")
  set_default_and_from_env(CTEST_DROP_SITE "testing.sandia.gov/cdash")
  set_default_and_from_env(CTEST_PROJECT_NAME "TribitsExampleMetaProject")
  set_default_and_from_env(CTEST_DROP_LOCATION "/submit.php?project=TribitsExampleMetaProject")
  set_default_and_from_env(CTEST_TRIGGER_SITE "")
  set_default_and_from_env(CTEST_DROP_SITE_CDASH TRUE)
endif()
