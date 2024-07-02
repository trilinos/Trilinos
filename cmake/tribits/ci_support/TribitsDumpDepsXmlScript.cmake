# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# CMake script file to generate an xml dependencies file for PROJECT_NAME
#
# This script is designed to be run as:
#
#   $ cmake \
#       [-D PROJECT_SOURCE_DIR=<projectSourceDir>] \
#       [-D <projectName>_PRE_REPOSITORIES=<prepo0>,<prepo1>,...] \
#       [-D <projectName>_EXTRA_REPOSITORIES=<erepo0>,<erepo1>,...] \
#       -D <projectName>_DEPS_XML_OUTPUT_FILE=<projectDepsFileOut> \
#       -P <tribitsDir>/ci_support/TribitsDumpDepsXmlScript.cmake
#

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

# A) Echo input options (must be specified with -D arguments to CMake command)

# PROJECT_SOURCE_DIR
message("Input: PROJECT_SOURCE_DIR = '${PROJECT_SOURCE_DIR}'")
if ("${PROJECT_SOURCE_DIR}" STREQUAL "")
  get_filename_component( DEFAULT_PROJECT_SOURCE_DIR
     "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE )
  message("-- DEFAULT_PROJECT_SOURCE_DIR='${DEFAULT_PROJECT_SOURCE_DIR}'")
  if (EXISTS "${DEFAULT_PROJECT_SOURCE_DIR}/ProjectName.cmake")
    message("-- Setting default PROJECT_SOURCE_DIR=${DEFAULT_PROJECT_SOURCE_DIR}")
    set(PROJECT_SOURCE_DIR "${DEFAULT_PROJECT_SOURCE_DIR}")
  else()
    message(FATAL_ERROR
      "ERROR: Cannot determine a default PROJECT_SOURCE_DIR location, please set PROJECT_SOURCE_DIR!") 
  endif()
else()
  set(PROJECT_NAME_FILE "${PROJECT_SOURCE_DIR}/ProjectName.cmake")
  if (NOT EXISTS "${PROJECT_NAME_FILE}")
    message(FATAL_ERROR
      "ERROR: PROJECT_SOURCE_DIR='${PROJECT_SOURCE_DIR}' is not a TriBITS project"
      " base dir since it is missing the file ProjectName.cmake!") 
  endif()
  # Else, this is a valid project source dir location
endif()

# Read in ProjectName.cmake to get PROJECT_NAME var
set(PROJECT_NAME_FILE "${PROJECT_SOURCE_DIR}/ProjectName.cmake")
include("${PROJECT_NAME_FILE}")

# Print the other vars
message("Input: ${PROJECT_NAME}_PRE_REPOSITORIES = '${${PROJECT_NAME}_PRE_REPOSITORIES}'")
message("Input: ${PROJECT_NAME}_EXTRA_REPOSITORIES = '${${PROJECT_NAME}_EXTRA_REPOSITORIES}'")
message("Input: ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE = '${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}'")

if ("${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}" STREQUAL "")
    message(FATAL_ERROR
      "ERROR: ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE cannot be empty."
      "  Please set ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE!") 
endif()

#
# Execute the rest of the script now that everything has been asserted or found
#

# Get the TRIBITS_DIR (we can always find this easy since this scrit is in TriBITS)
get_filename_component( ${PROJECT_NAME}_TRIBITS_DIR  "${CMAKE_CURRENT_LIST_DIR}/.."  ABSOLUTE )
message("-- Setting ${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}")

include("${CMAKE_CURRENT_LIST_DIR}/../core/common/TribitsConstants.cmake")
tribits_asesrt_minimum_cmake_version()
include("${CMAKE_CURRENT_LIST_DIR}/../core/common/TribitsCMakePolicies.cmake"  NO_POLICY_SCOPE)

set( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support"
  )

include(TribitsGlobalMacros)
include(TribitsPrintDependencyInfo)
include(TribitsWriteXmlDependenciesFiles)

# Generate the dependencies file

set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES  OFF)
set(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES  FALSE)
if (NOT ${PROJECT_NAME}_PRE_REPOSITORIES) # Make sure is defined!
  set(${PROJECT_NAME}_PRE_REPOSITORIES "")
endif()
if (NOT ${PROJECT_NAME}_EXTRA_REPOSITORIES) # Make sure is defined!
  set(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
endif()
tribits_read_in_native_repositories()
tribits_combine_native_and_extra_repos()
tribits_read_all_project_deps_files_create_deps_graph()
tribits_print_initial_dependency_info()
tribits_write_xml_dependency_files()
