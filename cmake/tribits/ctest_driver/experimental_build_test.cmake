# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
# where PACKAGES is the semi-colon-separated list of packages being tested
# (e.g. ${PROJECT_NAME}_PACKAGES="Teuchos;Epetra;NOX") and
# ${PROJECT_NAME}_TRIBITS_DIR points back to your home project directory.  You
# can take off the -VV argument if you don't want this to be too verbose.
#
# There are a number of other options that you can change as
# environment variables.  See the macros set_default_and_from_env(...)
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

set( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}"
  "${CTEST_SCRIPT_DIRECTORY}/../utils"
  )

include(TribitsCTestDriverCore)
include(GetLastDirName)
include(SetDefaultAndFromEnv)

#
# Override some configuration variables
#

# All these can be changed by env vars
set(CTEST_TEST_TYPE Experimental)
set(CTEST_UPDATE_VERSION_ONLY TRUE)
set(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "-Werror")

# Don't change these in the env!
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
set(CTEST_GENERATE_DEPS_XML_OUTPUT_FILE TRUE)
set(CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE FALSE)
set(CTEST_WIPE_CACHE FALSE)

# This script should be in PROJECT_BASE/cmake/tribits/ctest
set_default_and_from_env(PROJECT_SOURCE_DIR "${CTEST_SCRIPT_DIRECTORY}/../../..")
set(CTEST_SOURCE_DIRECTORY "${PROJECT_SOURCE_DIR}")

get_filename_component(PWD . REALPATH)
set(CTEST_BINARY_DIRECTORY "${PWD}")
set(CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/do-configure")

get_last_dir_name("${CTEST_BINARY_DIRECTORY}" BUILD_DIR_NAME)

# Can be overridden by the environment
set( CTEST_BUILD_NAME "${HOST_TYPE}-${BUILD_DIR_NAME}" )
set( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES OFF )

#
# Run the build/test/submit driver
#

tribits_ctest_driver()
