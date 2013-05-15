# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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
#   $ cmake -P -DPACKAGE_NAME=<packageName> [options] \
#       SOME_BASE_DIR/TribitsDumpDepsXmlScript.cmake
#

# A) Echo input options (must be specified with -D arguments to CMake command)

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")
MESSAGE("${PROJECT_NAME}_EXTRA_REPOSITORIES = ${${PROJECT_NAME}_EXTRA_REPOSITORIES}")
MESSAGE("${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR = ${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}")
# Get the utility CMake code that can determine the dependencies

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch"
  )

INCLUDE(TribitsGlobalMacros)

# Generate the dependencies file

SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)
SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
TRIBITS_READ_IN_NATIVE_REPOSITORIES()
  SET(${PROJECT_NAME}_ALL_REPOSITORIES ${${PROJECT_NAME}_NATIVE_REPOSITORIES}
    ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
