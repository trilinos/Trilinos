#
# Backward compatibility script for Trilinos CTest/CDash drivers
#

get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Used by TribitsCTestCoreDriver.cmake
SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
#MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

INCLUDE("${TRIBITS_PROJECT_ROOT}/cmake/tribits/ctest/TribitsCTestDriverCore.cmake")

macro(TRILINOS_CTEST_DRIVER)
  TRIBITS_CTEST_DRIVER()
endmacro()
