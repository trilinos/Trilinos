#
# Backward compatibility script for Trilinos CTest/CDash drivers
#

get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Set project defaults (only the defaults for variables set in
# TribitsCTestDriverCore.cmake)

SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
#MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

SET(Trilinos_REPOSITORY_LOCATION_DEFAULT
  "https://github.com/trilinos/Trilinos.git")
SET(Trilinos_REPOSITORY_LOCATION_NIGHTLY_DEFAULT
  "software.sandia.gov:/space/git/nightly/Trilinos")
SET(CTEST_DROP_SITE_COVERAGE_DEFAULT
  "testing.sandia.gov")
SET(CTEST_DROP_LOCATION_COVERAGE_DEFAULT
  "/extended/cdash/submit.php?project=Trilinos")

SET(Trilinos_REPOS_URL_BASE  https://github.com/trilinos/
  CACHE STRING "Set in TrilinosCTestDriverCore.cmake")

# Many of the existing scripts use the variable TRILINOS_CMAKE_DIR, so we set
# it here.
SET(TRILINOS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/../)

IF (NOT "$ENV{Trilinos_TRIBITS_DIR_OVERRIDE}" STREQUAL "")
  SET(Trilinos_TRIBITS_DIR "$ENV{Trilinos_TRIBITS_DIR_OVERRIDE}")
ELSE()
  SET(Trilinos_TRIBITS_DIR "${TRIBITS_PROJECT_ROOT}/cmake/tribits")
ENDIF()

INCLUDE("${Trilinos_TRIBITS_DIR}/ctest_driver/TribitsCTestDriverCore.cmake")

macro(TRILINOS_CTEST_DRIVER)
  TRIBITS_CTEST_DRIVER()
endmacro()
