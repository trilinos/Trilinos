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

# A) Echo input options (must be specified with -D arguments to CMake command)

# PROJECT_SOURCE_DIR
MESSAGE("Input: PROJECT_SOURCE_DIR = '${PROJECT_SOURCE_DIR}'")
IF ("${PROJECT_SOURCE_DIR}" STREQUAL "")
  GET_FILENAME_COMPONENT( DEFAULT_PROJECT_SOURCE_DIR
     "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE )
  MESSAGE("-- DEFAULT_PROJECT_SOURCE_DIR='${DEFAULT_PROJECT_SOURCE_DIR}'")
  IF (EXISTS "${DEFAULT_PROJECT_SOURCE_DIR}/ProjectName.cmake")
    MESSAGE("-- Setting default PROJECT_SOURCE_DIR=${DEFAULT_PROJECT_SOURCE_DIR}")
    SET(PROJECT_SOURCE_DIR "${DEFAULT_PROJECT_SOURCE_DIR}")
  ELSE()
    MESSAGE(FATAL_ERROR
      "ERROR: Cannot determine a default PROJECT_SOURCE_DIR location, please set PROJECT_SOURCE_DIR!") 
  ENDIF()
ELSE()
  SET(PROJECT_NAME_FILE "${PROJECT_SOURCE_DIR}/ProjectName.cmake")
  IF (NOT EXISTS "${PROJECT_NAME_FILE}")
    MESSAGE(FATAL_ERROR
      "ERROR: PROJECT_SOURCE_DIR='${PROJECT_SOURCE_DIR}' is not a TriBITS project"
      " base dir since it is missing the file ProjectName.cmake!") 
  ENDIF()
  # Else, this is a valid project source dir location
ENDIF()

# Read in ProjectName.cmake to get PROJECT_NAME var
SET(PROJECT_NAME_FILE "${PROJECT_SOURCE_DIR}/ProjectName.cmake")
INCLUDE("${PROJECT_NAME_FILE}")

# Print the other vars
MESSAGE("Input: ${PROJECT_NAME}_PRE_REPOSITORIES = '${${PROJECT_NAME}_PRE_REPOSITORIES}'")
MESSAGE("Input: ${PROJECT_NAME}_EXTRA_REPOSITORIES = '${${PROJECT_NAME}_EXTRA_REPOSITORIES}'")
MESSAGE("Input: ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE = '${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}'")

IF ("${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "ERROR: ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE cannot be empty."
      "  Please set ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE!") 
ENDIF()

#
# Execute the rest of the script now that everything has been asserted or found
#

# Get the TRIBITS_DIR (we can always find this easy since this scrit is in TriBITS)
GET_FILENAME_COMPONENT( ${PROJECT_NAME}_TRIBITS_DIR  "${CMAKE_CURRENT_LIST_DIR}/.."  ABSOLUTE )
MESSAGE("-- Setting ${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}")

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  )

INCLUDE(TribitsConstants)
TRIBITS_ASESRT_MINIMUM_CMAKE_VERSION()
INCLUDE(TribitsCMakePolicies)

INCLUDE(TribitsGlobalMacros)

# Generate the dependencies file

SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
IF (NOT ${PROJECT_NAME}_PRE_REPOSITORIES) # Make sure is defined!
  SET(${PROJECT_NAME}_PRE_REPOSITORIES "")
ENDIF()
IF (NOT ${PROJECT_NAME}_EXTRA_REPOSITORIES) # Make sure is defined!
  SET(${PROJECT_NAME}_EXTRA_REPOSITORIES "")
ENDIF()
TRIBITS_READ_IN_NATIVE_REPOSITORIES()
TRIBITS_COMBINE_NATIVE_AND_EXTRA_REPOS()
TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
