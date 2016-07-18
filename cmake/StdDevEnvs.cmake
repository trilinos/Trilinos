#
# This file automatically configures a TriBITS project (e.g. Trilinos) to use
# standard development environments.  Just include it using:
#
#  -D <Project>_CONFIGURE_OPTIONS_FILE=<some-path>/StdDevEnvs.cmake
#

SET(${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT "")

IF (
  (NOT "$ENV{TRILINOS_SEMS_DEV_ENV_LOADED}" STREQUAL "")
  AND
  ("${${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT}" STREQUAL "")
  AND
  ("${${PROJECT_NAME}_USE_BUILD_ENV}" STREQUAL "")
  )
  MESSAGE("-- SEMS: Env var"
    " TRILINOS_SEMS_DEV_ENV_LOADED={$ENV{TRILINOS_SEMS_DEV_ENV_LOADED}}"
    " => Allowing load of SEMS Build Env by default!" )
  SET(${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT "SEMS")
ENDIF()

IF (
  (NOT "$ENV{ATTB_ENV}" STREQUAL "")
  AND
  ("${${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT}" STREQUAL "")
  AND
  ("${${PROJECT_NAME}_USE_BUILD_ENV}" STREQUAL "")
  )
  #MESSAGE("-- ATTB: Env var ATTB_ENV=$ENV{ATTB_ENV}"
  #  " => Allowing load of ATTB Build Env by default!")
  #SET(${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT "ATTB")
  # ToDo: Turn the default back on again once we get working on more platforms
ENDIF()

SET(${PROJECT_NAME}_USE_BUILD_ENV_VALID_CHOICES "'', 'SEMS', or 'ATTB'")

SET(${PROJECT_NAME}_USE_BUILD_ENV ${${PROJECT_NAME}_USE_BUILD_ENV_DEFAULT}
  CACHE STRING
   "Load known build env (Choose ${${PROJECT_NAME}_USE_BUILD_ENV_VALID_CHOICES})")

PRINT_VAR(${PROJECT_NAME}_USE_BUILD_ENV)

IF (${PROJECT_NAME}_USE_BUILD_ENV STREQUAL "SEMS")
  MESSAGE("-- SEMS: Loading SEMSDevEnv.cmake to set compilers and TPL paths"
    " (To skip, set -D${PROJECT_NAME}_USE_BUILD_ENV=) ..." )
   SET(SEMS_DEV_ENV_FILE ${CMAKE_CURRENT_LIST_DIR}/std/sems/SEMSDevEnv.cmake)
  TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${SEMS_DEV_ENV_FILE}")
  INCLUDE(${SEMS_DEV_ENV_FILE})
ELSEIF (${PROJECT_NAME}_USE_BUILD_ENV STREQUAL "ATTB")
  MESSAGE("-- ATTB: Loading ATTBDevEnv.cmake to set compilers and TPL paths"
    " (To skip, set -D${PROJECT_NAME}_USE_BUILD_ENV=) ..." )
  SET(ATTB_DEV_ENV_FILE ${CMAKE_CURRENT_LIST_DIR}/ctest/drivers/ATTB/ATTBDevEnv.cmake)
  TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${ATTB_DEV_ENV_FILE}")
  INCLUDE(${ATTB_DEV_ENV_FILE})
ELSEIF (${PROJECT_NAME}_USE_BUILD_ENV STREQUAL "")
  # Don't load any known dev env 
ELSE()
  MESSAGE(FATAL_ERROR "ERROR: The value of"
    " `${PROJECT_NAME}_USE_BUILD_ENV='${${PROJECT_NAME}_USE_BUILD_ENV}'"
    " is invalid. Choose one of the supported values"
    " ${${PROJECT_NAME}_USE_BUILD_ENV_VALID_CHOICES}!")
ENDIF()

# ToDo: Add more standard dev envs ...
