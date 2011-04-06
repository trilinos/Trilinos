
#
# CMake script file to generate an xml dependencies file for PROJECT_NAME
#
# This script is designed to be run as:
#
#   $ cmake -P -DPACKAGE_NAME=<packageName> [options] SOME_BASE_DIR/PackageArchDumpDepsXmlScript.cmake
#

# A) Echo input options (must be specified with -D arguments to CMake command)

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("${PROJECT_NAME}_CMAKE_DIR = ${${PROJECT_NAME}_CMAKE_DIR}")
MESSAGE("${PROJECT_NAME}_DEPS_HOME_DIR = ${${PROJECT_NAME}_DEPS_HOME_DIR}")
MESSAGE("${PROJECT_NAME}_EXTRA_REPOSITORIES = ${${PROJECT_NAME}_EXTRA_REPOSITORIES}")
MESSAGE("${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR = ${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}")
# Get the utility CMake code that can determine the dependencies

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_CMAKE_DIR}/utils"
  "${${PROJECT_NAME}_CMAKE_DIR}/package_arch"
  )

INCLUDE(PackageArchGlobalMacros)

# Generate the d

SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)
SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)

PACKAGE_ARCH_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
