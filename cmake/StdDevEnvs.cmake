#
# This file automatically configures a TriBITS project (e.g. Trilinos) to use
# standard development environments.  Just include it using:
#
#  -D <Project>_CONFIGURE_OPTIONS_FILE=<some-path>/StdDevEnvs.cmake
#

IF (NOT $ENV{ATTB_ENV} STREQUAL "")
  SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT TRUE)
  IF ("${${PROJECT_NAME}_USE_ATTB_DEV_ENV}" STREQUAL "")
    MESSAGE("-- NOTE: Env var ATTB_ENV=$ENV{ATTB_ENV} => Including ATTBDevEnv.cmake "
      " by default!  (To skip, set -D${PROJECT_NAME}_USE_ATTB_DEV_ENV=FALSE)")
  ENDIF()
ELSE()
  SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT FALSE)
ENDIF()

SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV ${${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT}
  CACHE BOOL "Set if the ATTBDevEnv.cmake file should be included.")

IF (${PROJECT_NAME}_USE_ATTB_DEV_ENV)
  SET(ATTB_DEV_ENV_FILE ${CMAKE_CURRENT_LIST_DIR}/ctest/drivers/ATTB/ATTBDevEnv.cmake)
  TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${ATTB_DEV_ENV_FILE}")
  INCLUDE(${ATTB_DEV_ENV_FILE})
ENDIF()

# ToDo: Add more standard dev envs like SEMSDevEnv.cmake ...
