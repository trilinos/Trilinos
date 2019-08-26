#
# cmake -P script to get the CTest testing XML directory
# <build>/Testing/<buildstarttime> given just the <bulid> directory path.
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

IF ("${PROJECT_NAME}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error, PROJECT_NAME must be set!")
ENDIF()

IF ("${${PROJECT_NAME}_TRIBITS_DIR}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error, ${PROJECT_NAME}_TRIBITS_DIR must be set!")
ENDIF()

IF ("${CTEST_BUILD_DIR}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error, CTEST_BUILD_DIR must be set!")
ENDIF()

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ctest_driver"
  )

INCLUDE(TribitsReadTagFile)

SET(TAG_FILE "${CTEST_BUILD_DIR}/Testing/TAG")

TRIBITS_READ_CTEST_TAG_FILE("${TAG_FILE}" BUILD_START_TIME  CDASH_TRACK)

MESSAGE("${CTEST_BUILD_DIR}/Testing/${BUILD_START_TIME}")
