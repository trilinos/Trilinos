################################################################################
#                                                                              #
#                               TriBITS                                        #
#                                                                              #
################################################################################


IF (NOT TriBITS_PROCESSED_BASE_PROJECT AND NOT PACKAGE_NAME STREQUAL "TriBITS")

  # This CMakeLists.txt file is being processed as the TriBITS projects's base
  # CMakeLists.txt file!  (See comments at bottom of this file.)

  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
  SET(${PROJECT_NAME}_VERSION 1.0)
  PROJECT(${PROJECT_NAME} NONE)
  
  SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "")
  INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
  
  CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

  # TriBITS is the only package so just enable it and test it!
  SET(${PROJECT_NAME}_ENABLE_TriBITS ON)
  SET(${PROJECT_NAME}_ENABLE_TESTS ON CACHE BOOL "Enable test by default.")

  SET(TriBITS_PROCESSED_BASE_PROJECT ON)

  TRIBITS_PROJECT()

ELSE()

  # This CMakeLists.txt file is being processed as the TriBITS package file.

  TRIBITS_PACKAGE(TriBITS)
  
  TRIBITS_ADD_TEST_DIRECTORIES(
    python
    package_arch
    )
  
  ASSERT_DEFINED(CMAKE_GENERATOR)
  IF (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    TRIBITS_ADD_TEST_DIRECTORIES(
      ctest/UnitTests
      doc/examples/UnitTests
      )
  ENDIF()
  
  TRIBITS_PACKAGE_POSTPROCESS()

ENDIF()

# NOTE: In order to allow the `TriBITS/` directory to be both a TriBITS
# package (for inclusion in other TriBITS projects) and to be a TriBITS
# project itself, you only have to put in a simple if statement in this
# top-level CMakeLists.txt file.  That is all!  In fact, this `TriBiTS/`
# directory shows how the same directory can be used for a TriBITS package, a
# TriBITS repository, and a TriBITS project!
