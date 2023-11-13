################################################################################
#                                                                              #
#                        TribitsExampleProject2                                #
#                                                                              #
################################################################################

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)
include("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
project(${PROJECT_NAME} LANGUAGES NONE)
set(TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE  TRUE)
set(${PROJECT_NAME}_TRIBITS_DIR
   "${CMAKE_CURRENT_LIST_DIR}/../.."  CACHE  STRING
  "TriBITS base directory (default assumes in TriBITS source tree)")
include("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
tribits_project()
