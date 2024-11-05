#
# Set the locations of things for this project
#


set(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
set(CTEST_SOURCE_NAME "TribitsExampleProject")
include("${TRIBITS_PROJECT_ROOT}/ProjectName.cmake")
if (NOT "$ENV{${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  set(${PROJECT_NAME}_TRIBITS_DIR "$ENV{${PROJECT_NAME}_TRIBITS_DIR}")
endif()
if("${${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  # If not set externally, then assume this is inside of tribits example
  # directory.
  set(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../..")
endif()

#
# Include the TriBITS file to get other modules included
#

include("${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver/TribitsCTestDriverCore.cmake")

function(tribitsexproj_ctest_driver)
  set_default_and_from_env( CTEST_BUILD_FLAGS "-j1 -i" )
  set_default_and_from_env( CTEST_PARALLEL_LEVEL "1" )
  tribits_ctest_driver()
endfunction()
