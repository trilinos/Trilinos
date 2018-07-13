#
# Set the locations of things for this project
#


SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
SET(CTEST_SOURCE_NAME "TribitsExampleProject")
INCLUDE("${TRIBITS_PROJECT_ROOT}/ProjectName.cmake")
IF (NOT "$ENV{${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  SET(${PROJECT_NAME}_TRIBITS_DIR "$ENV{${PROJECT_NAME}_TRIBITS_DIR}")
ENDIF()
IF("${${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  # If not set externally, then assume this is inside of tribits example
  # directory.
  SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../..")
ENDIF()

#
# Include the TriBITS file to get other modules included
#

INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver/TribitsCTestDriverCore.cmake")

FUNCTION(TRIBITSEXPROJ_CTEST_DRIVER)
  SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "-j1 -i" )
  SET_DEFAULT_AND_FROM_ENV( CTEST_PARALLEL_LEVEL "1" )
  TRIBITS_CTEST_DRIVER()
ENDFUNCTION()
