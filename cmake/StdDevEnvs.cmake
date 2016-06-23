#
# This file automatically configures a TriBITS project (e.g. Trilinos) to use
# standard development environments.  Just include it using:
#
#  -D <Project>_CONFIGURE_OPTIONS_FILE=<some-path>/StdDevEnvs.cmake
#

IF (NOT "$ENV{TRILINOS_SEMS_DEV_ENV_LOADED}" STREQUAL "")
  SET(${PROJECT_NAME}_USE_SEMS_DEV_ENV_DEFAULT TRUE)
  IF ("${${PROJECT_NAME}_USE_SEMS_DEV_ENV}" STREQUAL "")
    MESSAGE("-- NOTE: Env var"
       " TRILINOS_SEMS_DEV_ENV_LOADED={$ENV{TRILINOS_SEMS_DEV_ENV_LOADED}}"
       " => Including SEMSDevEnv.cmake by default!" )
  ENDIF()
ELSE()
  SET(${PROJECT_NAME}_USE_SEMS_DEV_ENV_DEFAULT FALSE)
ENDIF()

IF (NOT $ENV{ATTB_ENV} STREQUAL "")
  #SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT TRUE)
  #IF ("${${PROJECT_NAME}_USE_ATTB_DEV_ENV}" STREQUAL "")
  #  MESSAGE("-- NOTE: Env var ATTB_ENV=$ENV{ATTB_ENV} => Including ATTBDevEnv.cmake "
  #    " by default!  (To skip, set -D${PROJECT_NAME}_USE_ATTB_DEV_ENV=FALSE)")
  #ENDIF()
  SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT FALSE)
  # ToDo: Turn the default back on again once we get working on more platforms
ELSE()
  SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT FALSE)
ENDIF()

SET(${PROJECT_NAME}_USE_SEMS_DEV_ENV ${${PROJECT_NAME}_USE_SEMS_DEV_ENV_DEFAULT}
  CACHE BOOL "Set if the SEMSDevEnv.cmake file should be included.")

SET(${PROJECT_NAME}_USE_ATTB_DEV_ENV ${${PROJECT_NAME}_USE_ATTB_DEV_ENV_DEFAULT}
  CACHE BOOL "Set if the ATTBDevEnv.cmake file should be included.")

IF (${PROJECT_NAME}_USE_SEMS_DEV_ENV)
  MESSAGE("-- SEMS: Setting compilers and TPL paths for loaded SEMS Dev Env"
    " (To skip, set -D${PROJECT_NAME}_USE_SEMS_DEV_ENV=FALSE) ..." )
   SET(SEMS_DEV_ENV_FILE ${CMAKE_CURRENT_LIST_DIR}/std/sems/SEMSDevEnv.cmake)
  TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${SEMS_DEV_ENV_FILE}")
  INCLUDE(${SEMS_DEV_ENV_FILE})
ELSEIF (${PROJECT_NAME}_USE_ATTB_DEV_ENV)
  SET(ATTB_DEV_ENV_FILE ${CMAKE_CURRENT_LIST_DIR}/ctest/drivers/ATTB/ATTBDevEnv.cmake)
  TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${ATTB_DEV_ENV_FILE}")
  INCLUDE(${ATTB_DEV_ENV_FILE})
ENDIF()

# ToDo: Add more standard dev envs ...
