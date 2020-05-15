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
# This file is included from TribitsProject.cmake which has already set to the
# tribits implementation to use in {PROJECT_NAME}_TRIBITS_DIR.
#


SET(CMAKE_MODULE_PATH
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/cmake
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/utils
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/config_tests
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/modules
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/installation
   )

IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("CMAKE_MODULE_PATH='${CMAKE_MODULE_PATH}'")
ENDIF()

INCLUDE(TribitsConstants)
TRIBITS_ASESRT_MINIMUM_CMAKE_VERSION()
INCLUDE(TribitsCMakePolicies)

INCLUDE(TribitsIncludeDirectories)
INCLUDE(TribitsFindPythonInterp)
INCLUDE(TribitsGlobalMacros)
INCLUDE(TribitsConfigureCTestCustom)

INCLUDE(AdvancedSet)
INCLUDE(AdvancedOption)
INCLUDE(TimingUtils)
INCLUDE(SetDefault)


#
# Defines a TriBITS project (the guts).
#

MACRO(TRIBITS_PROJECT_IMPL)

  #
  # A) Basic top-level TriBITS project stuff
  #

  MESSAGE("")
  MESSAGE("Configuring ${PROJECT_NAME} build directory")
  MESSAGE("")

  # A.1) Set some basic system vars and info you can't change
  TRIBITS_ASSERT_AND_SETUP_PROJECT_AND_STATIC_SYSTEM_VARS()

  # A.2) Read user provided options from specified files.  It is important to
  # process these files *very* early on so that they have the same basic
  # effect of setting these variables directly in the cache.
  TRIBITS_READ_IN_OPTIONS_FROM_FILE()

  # A.3) Get some other basic system info that is useful early on in
  # configuration
  TRIBITS_SETUP_BASIC_SYSTEM_VARS()
  TRIBITS_FIND_PYTHON_INTERP()
  FIND_PACKAGE(Git)
  IF (NOT "${GIT_VERSION_STRING_OVERRIDE}" STREQUAL "")
    SET(GIT_VERSION_STRING ${GIT_VERSION_STRING_OVERRIDE}) # For testing!
  ENDIF()

  #
  # A.4) Read in the Project's version file
  #
  # NOTE: The file Version.cmake must be read *before* the global options are
  # read because the variables defined in Version.cmake provide defaults for
  # many of these options.
  #
  TRIBITS_PROJECT_READ_VERSION_FILE(${PROJECT_SOURCE_DIR})

  # Since the version header file is now configured the root build
  # dir needs to be on the include path
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  #
  # B) Set up user options and global variables that will be used throughout
  #

  MESSAGE("")
  MESSAGE("Setting up major user options ...")
  MESSAGE("")

  TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS()

  # Have to start timing after we read in the major options since that
  # determines the timing option var.
  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(GLOBAL_TIME_START_SECONDS)
  ENDIF()

  TRIBITS_READ_IN_NATIVE_REPOSITORIES()

  TRIBITS_COMBINE_NATIVE_AND_EXTRA_REPOS()

  TRIBITS_PROCESS_EXTRA_REPOS_OPTIONS_FILES()

  INCLUDE(TribitsInstallationTestingMacros)
  TRIBITS_FIND_PROJECT_INSTALL()

  #
  # C) Generate version info file and read in ${PROJECT_NAME} packages and
  # TPLs and process dependencies
  #

  TRIBITS_GENERATE_REPO_VERSION_OUTPUT_AND_FILE_AND_INSTALL()

  # Read in and process all of the project's package, TPL, listss and
  # dependency definition files.  The order these are read in are:
  #
  # * Read in all PackagesList.cmake and TPLsList.cmake files for all
  #   native and extra repos in repo order.
  # * Process each repos's cmake/RepositoryDependenciesSetup.cmake file
  #   in repo order.
  # * Process the project's cmake/cmake/ProjectDependenciesSetup.cmake
  # * Process each package's Dependencies.cmake file.  If a package a subpackages,
  #   The subpackage Dependencies.cmake files are read before setting up the
  #   parent package's dependenices.
  #
  TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

  #
  # D) Apply dependency logic to enable and disable TriBITS packages packages
  # and tests
  #

  TRIBITS_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES()

  #
  # E) Stop after all dependencies handling is finished if asked.
  #

  IF (${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
    MESSAGE("")
    MESSAGE("Shortcircuiting after dependency tracking ...")
    RETURN()
  ENDIF()

  #
  # F) Set up the environment on this computer
  #

  MESSAGE("")
  MESSAGE("Probing the environment ...")
  MESSAGE("")

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    TRIBITS_SETUP_ENV()
  ELSE()
    MESSAGE("-- Skipping env setup due to"
      " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
  ENDIF()

  #
  # G) Go get the information for all enabled TPLS
  #

  MESSAGE("")
  MESSAGE("Getting information for all enabled TPLs ...")
  MESSAGE("")

  TRIBITS_PROCESS_ENABLED_TPLS()

  #
  # H) Set up for testing with CTest and ${PROJECT_NAME} test harness
  #

  MESSAGE("")
  MESSAGE("Setting up testing support ...")
  MESSAGE("")

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    TRIBITS_INCLUDE_CTEST_SUPPORT()
  ELSE()
    MESSAGE("-- Skipping testing support setup due to"
      " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
  ENDIF()

  #
  # I) Add the 'dashboard' target
  #
  # NOTE: Must come after setting up for testing
  #

  SET(TRIBITS_ADD_DASHBOARD_TARGET_MODULE
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CTEST_DRIVER_DIR}/TribitsAddDashboardTarget.cmake
    )
  IF (
    EXISTS ${TRIBITS_ADD_DASHBOARD_TARGET_MODULE}
    AND
    NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY
    )
    INCLUDE(${TRIBITS_ADD_DASHBOARD_TARGET_MODULE})
    TRIBITS_ADD_DASHBOARD_TARGET()
  ENDIF()

  #
  # J) Configure individual packages
  #

  MESSAGE("")
  MESSAGE("Configuring individual enabled ${PROJECT_NAME} packages ...")
  MESSAGE("")

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    TRIBITS_REPOSITORY_CONFIGURE_ALL_VERSION_HEADER_FILES(
      ${${PROJECT_NAME}_ALL_REPOSITORIES})
    TRIBITS_REPOSITORY_CONFIGURE_ALL_VERSION_DATE_FILES(
      ${${PROJECT_NAME}_ALL_REPOSITORIES})
  ENDIF()

  TRIBITS_CONFIGURE_ENABLED_PACKAGES()

  #
  # K) Write dummy makefiles for Ninja
  #

  IF (CMAKE_GENERATOR STREQUAL "Ninja")

    IF (${PROJECT_NAME}_WRITE_NINJA_MAKEFILES)

      MESSAGE("")
      MESSAGE("Generating dummy makefiles in each directory to call Ninja ...")
      MESSAGE("")
  
      IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
        TIMER_GET_RAW_SECONDS(NINJA_MAKEFILES_TIME_START_SECONDS)
      ENDIF()
  
      INCLUDE(GenerateNinjaMakefiles)
      GENERATE_NINJA_MAKEFILES(${CMAKE_SOURCE_DIR})
  
      IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
        TIMER_GET_RAW_SECONDS(NINJA_MAKEFILES_TIME_END_SECONDS)
        TIMER_PRINT_REL_TIME(${NINJA_MAKEFILES_TIME_START_SECONDS}
         ${NINJA_MAKEFILES_TIME_END_SECONDS}
           "Total time generate Ninja makefiles ${PROJECT_NAME}")
      ENDIF()

    ELSE()

      MESSAGE("\nNOTE: *NOT* generating dummy Ninja makefiles (see above note"
        " and check CMake version)")

    ENDIF()


  ENDIF()

  #
  # L) Setup for packaging and distribution
  #

  IF (${PROJECT_NAME}_ENABLE_CPACK_PACKAGING)
    MESSAGE("")
    MESSAGE("Set up for creating a distribution ...")
    MESSAGE("")
    IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
      TRIBITS_SETUP_PACKAGING_AND_DISTRIBUTION()
    ELSE()
      MESSAGE("-- Skipping distribution setup due to"
        " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
    ENDIF()
  ELSE()
    MESSAGE("")
    MESSAGE("Skipping setup for distribution because"
      " ${PROJECT_NAME}_ENABLE_CPACK_PACKAGING=OFF")
    MESSAGE("")
  ENDIF()

  #
  # M) Set up for installation
  #

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    TRIBITS_SETUP_FOR_INSTALLATION()
  ENDIF()

  #
  # N) Show final timing and end
  #

  MESSAGE("")
  MESSAGE("Finished configuring ${PROJECT_NAME}!")
  MESSAGE("")
  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(GLOBAL_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${GLOBAL_TIME_START_SECONDS}  ${GLOBAL_TIME_STOP_SECONDS}
      "Total time to configure ${PROJECT_NAME}")
  ENDIF()

ENDMACRO()


#
# @MACRO: TRIBITS_PROJECT_ENABLE_ALL()
#
# Process a project where you enable all of the packages by default.
#
# Usage::
#
#   TRIBITS_PROJECT_ENABLE_ALL()
#
# This macro just sets the global cache var
# `${PROJECT_NAME}_ENABLE_ALL_PACKAGES`_ to ``ON`` by default then calls
# `TRIBITS_PROJECT()`_.  That is all.  This macro is generally used for
# TriBITS projects that have just a single package or by default just want to
# enable all packages.  This is especially useful when you have a TriBITS
# project with just a single package.
#
MACRO(TRIBITS_PROJECT_ENABLE_ALL)
  SET(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON CACHE BOOL "Enable all by default" )
  TRIBITS_PROJECT_IMPL(${ARGN})
ENDMACRO()
