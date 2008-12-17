SET(CTEST_PROJECT_NAME "Trilinos")
SET(CTEST_NIGHTLY_START_TIME "02:00:00 EST") # Midnight MST

IF(NOT DEFINED CTEST_DROP_METHOD)
  SET(CTEST_DROP_METHOD "http")
ENDIF()

IF(CTEST_DROP_METHOD STREQUAL "http")
  #SET(CTEST_DROP_SITE "datamining.sandia.gov")
  #SET(CTEST_DROP_LOCATION "/CDash/submit.php?project=Trilinos")
  SET(CTEST_DROP_SITE "trilinos.sandia.gov")
  SET(CTEST_DROP_LOCATION "/cdash/submit.php?project=Trilinos")
  SET(CTEST_TRIGGER_SITE "")
ENDIF()
