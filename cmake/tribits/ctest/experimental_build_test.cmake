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
# CTest script that is used to do an experimental build/test right from a
# developer's own build directory.
#
# To run this script:
#
# 1) First configure your directory without enabling any packages (to set the
# platform-specific and site-specific options).  The CMakeCache.txt file that
# is created will then be modified by the script as packages are enabled and
# disabled.  NOTE: For the safest results, start with an empty build directory
# (except for your configure script of course).  This seems to be needed when
# switching back and forth between an 'Experimental' and 'Nightly' build for
# instance.
#
# 2) Run the script (overriding any appropriate options) as:
#
#    env ${PROJECT_NAME}_PACKAGES="<PACKAGES>" \
#      ctest -S ${PROJECT_NAME}_TRIBITS_DIR/ctest/experimental_build_test.cmake -VV
#
# where PACAKGES is the semi-colon-separated list of packages being tested
# (e.g. ${PROJECT_NAME}_PACKAGES="Teuchos;Epetra;NOX") and
# ${PROJECT_NAME}_TRIBITS_DIR points back to your home project directory.  You
# can take off the -VV argument if you don't want this to be too verbose.
#
# There are a number of other options that you can change as
# environment varibles.  See the macros SET_DEFAULT_AND_FROM_ENV(...)
# in the file TribitsCTestDriverCore.cmake.  One option that you
# might want to overridde, for instance is CTEST_BUILD_NAME so that
# you can insert a special name into the dashboard.
#
# When this script finishes running, the last package listed in
# ${PROJECT_NAME}_PACAKGES will be enabled in the CMakeCache.txt file.
#
# NOTE: It is better to use the CMake-built make target 'dashboard' to run
# this script as it takes care of the details of manipulating the cache and
# restoring the package enables when it is done.
#

#
# General setup code:
#
# Do not modify any of this directly, use use environment variables instead!
#

#
# Include some CMake/CTest code files
#

SET( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}"
  "${CTEST_SCRIPT_DIRECTORY}/../utils"
  )

INCLUDE(TribitsCTestDriverCore)
INCLUDE(GetLastDirName)
INCLUDE(SetDefaultAndFromEnv)

#
# Override some configuration variables
#

# All these can be changed by env vars
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_DO_UPDATES FALSE)
SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "-Werror")

# Don't change these in the env!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
SET(CTEST_GENERATE_DEPS_XML_OUTPUT_FILE TRUE)
SET(CTEST_WIPE_CACHE FALSE)

# This script should be in PROJECT_BASE/cmake/tribits/ctest
SET_DEFAULT_AND_FROM_ENV(PROJECT_SOURCE_DIR "${CTEST_SCRIPT_DIRECTORY}/../../..")
SET(CTEST_SOURCE_DIRECTORY "${PROJECT_SOURCE_DIR}")

GET_FILENAME_COMPONENT(PWD . REALPATH)
SET(CTEST_BINARY_DIRECTORY "${PWD}")
SET(CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/do-configure")

GET_LAST_DIR_NAME("${CTEST_BINARY_DIRECTORY}" BUILD_DIR_NAME)

# Can be overridden by the environment
SET( CTEST_BUILD_NAME "${HOST_TYPE}-${BUILD_DIR_NAME}" )
SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES OFF )

#
# Run the build/test/submit driver
#

TRIBITS_CTEST_DRIVER()
