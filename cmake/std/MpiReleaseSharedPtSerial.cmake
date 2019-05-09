# This file can be read in using either:
#
#   -C <abs-path>/<file-name>.cmake
#
# or:
#
#   -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/<file-name>.cmake
#

# Handle this being passed in with -C option instead of
# <Project>_CONFIGURE_OPTIONS_FILE.
IF ("${PROJECT_NAME}" STREQUAL "")
  SET(PROJECT_NAME Trilinos)
  INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../tribits/core/utils/AssertDefined.cmake")
ENDIF()

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SetUtils.cmake")

TRIL_SET_BOOL_CACHE_VAR(TPL_ENABLE_MPI ON)
TRIL_SET_BOOL_CACHE_VAR(CMAKE_BUILD_TYPE RELEASE)
TRIL_SET_BOOL_CACHE_VAR(Trilinos_ENABLE_DEBUG OFF)
TRIL_SET_BOOL_CACHE_VAR(BUILD_SHARED_LIBS ON)
TRIL_SET_BOOL_CACHE_VAR(Trilinos_ENABLE_EXPLICIT_INSTANTIATION ON)
TRIL_SET_BOOL_CACHE_VAR(Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)

# NOTE: The order of these includes matters!

include("${CMAKE_CURRENT_LIST_DIR}/BasicCiTestingSettings.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/sems/SEMSDevEnv.cmake")
