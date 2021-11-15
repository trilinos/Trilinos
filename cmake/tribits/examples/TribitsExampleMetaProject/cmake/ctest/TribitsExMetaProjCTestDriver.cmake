#
# Set the locations of things for this project
#

set(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
set(CTEST_SOURCE_NAME "TribitsExampleMetaProject")
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

function(tribitsexmetaproj_ctest_driver)
  set_default_and_from_env(TribitsExMetaProj_GIT_URL_REPO_BASE
    https://github.com/tribits/)
  set_default(TribitsExMetaProj_REPOSITORY_LOCATION_DEFAULT
    "${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleMetaProject.git")
  set_default(TribitsExMetaProj_REPOSITORY_LOCATION_NIGHTLY_DEFAULT 
    "${TribitsExMetaProj_REPOSITORY_LOCATION_DEFAULT}")
  tribits_ctest_driver()
endfunction()
