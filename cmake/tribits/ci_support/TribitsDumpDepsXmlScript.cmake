# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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

cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)

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

set( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  "${${PROJECT_NAME}_TRIBITS_DIR}/ci_support"
  )

include(TribitsConstants)
tribits_asesrt_minimum_cmake_version()
include(TribitsCMakePolicies  NO_POLICY_SCOPE)

include(TribitsGlobalMacros)
include(TribitsPrintDependencyInfo)
include(TribitsWriteXmlDependenciesFiles)

# Generate the dependencies file

set(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
set(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
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
