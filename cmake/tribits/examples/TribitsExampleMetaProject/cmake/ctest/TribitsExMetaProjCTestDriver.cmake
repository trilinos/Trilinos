#
# Set the locations of things for this project
#

SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
SET(CTEST_SOURCE_NAME "TribitsExampleMetaProject")
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

FUNCTION(TRIBITSEXMETAPROJ_CTEST_DRIVER)
  SET_DEFAULT_AND_FROM_ENV(TribitsExMetaProj_GIT_URL_REPO_BASE
    https://github.com/tribits/)
  SET_DEFAULT(TribitsExMetaProj_REPOSITORY_LOCATION_DEFAULT
    "${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleMetaProject.git")
  SET_DEFAULT(TribitsExMetaProj_REPOSITORY_LOCATION_NIGHTLY_DEFAULT 
    "${TribitsExMetaProj_REPOSITORY_LOCATION_DEFAULT}")
  TRIBITS_CTEST_DRIVER()
ENDFUNCTION()
