#
# CMake -P script to return a 10-digit integer for the repo version at a given
# commit reference.
#
# Usage:
#
#   cmake \
#     -D PROJECT_NAME=<Project> \
#     -D <Project>_SOURCE_DIR=<projectDir> \
#     -D COMMIT_REF=<commit_ref> \
#     -P tribits_get_version_date.cmake
#   

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

# A) Validate input

if ("${PROJECT_NAME}" STREQUAL "")
  message(FATAL_ERROR "Error, must set PROJECT_NAME!")
endif()

if ("${${PROJECT_NAME}_SOURCE_DIR}" STREQUAL "")
  message(FATAL_ERROR "Error, must set ${PROJECT_NAME}_SOURCE_DIR!")
endif()

if ("${COMMIT_REF}" STREQUAL "")
  message(FATAL_ERROR "Error, must set COMMIT_REF!")
endif()

# B) Include modules

set(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")
set(CMAKE_MODULE_PATH
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/utils
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch
   )

include(TribitsGetVersionDate)

# C) Find Git

find_package(Git QUIET)

# D) Run functions

tribits_get_raw_git_commit_utc_time("${${PROJECT_NAME}_SOURCE_DIR}"  "${COMMIT_REF}"
  raw_git_commit_utc_time )

tribits_get_version_date_from_raw_git_commit_utc_time("${raw_git_commit_utc_time}"
  version_date_out )

message(${version_date_out})
