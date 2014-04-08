################################################################################
#                                                                              #
#                               TriBITS                                        #
#                                                                              #
################################################################################


IF (NOT TriBITS_PROCESSED_BASE_PROJECT AND NOT PACKAGE_NAME STREQUAL "TriBITS")

  # This CMakeLists.txt file is being processed as the TriBITSProj projects's
  # base CMakeLists.txt file!  (See comments at bottom of this file.)

  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
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

# NOTE: In order to allow the `tribits` directory to be both a TriBITS package
# (for inclusion in other TriBITS projects) and to be a TriBITS project
# itself, you only have to put in a simple if statement in this top-level
# CMakeLists.txt file and need to use a different package binary directory in
# PackagesList.cmake.  That is all!  In fact, this `tribits` directory shows
# how the same directory can be used for a TriBITS package, a TriBITS
# repository, and a TriBITS project!  However, you can't name the project the
# same as the package, otherwise, TriBITS woiuld try to try to create two
# targets with the name, such as ${PROJECT_NAME}_libs.  Therefore a different
# name and different binary directory must be used for the TriBITS package.
# But that is all.
