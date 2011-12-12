#
# Backward compatibility script for Trilinos CTest/CDash drivers
#

get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Used by TribitsCTestCoreDriver.cmake
SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
#MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

# Many of the existing scripts use the variable TRILINOS_CMAKE_DIR, so we set
# it here.
SET(TRILINOS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/../)

INCLUDE("${TRIBITS_PROJECT_ROOT}/cmake/tribits/ctest/TribitsCTestDriverCore.cmake")

macro(TRILINOS_CTEST_DRIVER)
  TRIBITS_CTEST_DRIVER()
endmacro()
