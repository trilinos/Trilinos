################################################################################
#                                                                              #
#                               TriBITS                                        #
#                                                                              #
################################################################################


IF (TriBITS_PROCESSED_BASED_PROJECT)

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

ELSE()

  # This CMakeLists.txt file is being processed as the projects's base
  # CMakeLists.txt file!  (See comments at bottom of file)

  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
  PROJECT(${PROJECT_NAME} NONE)
  
  SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "")
  SET(CMAKE_MODULE_PATH "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch")
  INCLUDE(TribitsProject)
  
  CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

  # TriBITS is the only package so just enabled it!
  SET(${PROJECT_NAME}_ENABLE_TriBITS ON)

  SET(TriBITS_PROCESSED_BASED_PROJECT ON)

  TRIBITS_PROJECT()

ENDIF()

# NOTE: In order to allow the tribits directory to be both a TriBITS package
# (for inclusion in other TriBITS projects) and to be a TriBITS project
# itself, you only have to put in a simple if statement and in
# PackagesList.cmake you just need to use a different binary directory.  That
# is all!  In fact, this tribits directory shows how the same directory can be
# used for a TriBITS package, a TriBITS repository, and a TriBITS project!
# However, you can't name the project the same as the package since this will
# try to create two targets with the name ${PROJECT_NAME}_libs.
