# This file is provided solely for backwards compatibility with
# existing driver CTest scripts.

# Existing scripts need a value for TRILINOS_CMAKE_DIR.
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(TRILINOS_CMAKE_DIR "${CMAKE_CURRENT_LIST_DIR}/../")

#
# Include the real TriBITS driver script.
#
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../tribits/ctest/TribitsCTestDriverCore.cmake")
macro(TRILINOS_CTEST_DRIVER)
  TRIBITS_CTEST_DRIVER()
endmacro()
