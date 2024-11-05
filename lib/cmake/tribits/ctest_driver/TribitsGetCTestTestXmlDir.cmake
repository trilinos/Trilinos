#
# cmake -P script to get the CTest testing XML directory
# <build>/Testing/<buildstarttime> given just the <build> directory path.
#
# Usage:
#
#   cmake \
#     -DPROJECT_NAME=<projectName> \
#     -D${PROJECT_NAME}_TRIBITS_DIR=<tribits-dir> \
#     -DCTEST_BUILD_DIR=<build-dir> \
#     -P <tribits-dir>/ctest_driver/TribitsGetCTestTestXmlDir.cmake
#
# This script reads in the <build-dir>/Testing/TAG to get <bulidstarttime> and
# then prints the directory <build>/Testing/<buildstarttime> to STDOUT.
#

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

if ("${PROJECT_NAME}" STREQUAL "")
  message(FATAL_ERROR "Error, PROJECT_NAME must be set!")
endif()

if ("${${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  message(FATAL_ERROR "Error, ${PROJECT_NAME}_TRIBITS_DIR must be set!")
endif()

if ("${CTEST_BUILD_DIR}" STREQUAL "")
  message(FATAL_ERROR "Error, CTEST_BUILD_DIR must be set!")
endif()

set( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver"
  )

include(TribitsReadTagFile)

set(TAG_FILE "${CTEST_BUILD_DIR}/Testing/TAG")

tribits_read_ctest_tag_file("${TAG_FILE}" buildStartTime  cdashGroup  cdashModel)

message("${CTEST_BUILD_DIR}/Testing/${buildStartTime}")
