# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# This file is included from TribitsProject.cmake which has already set to the
# tribits implementation to use in {PROJECT_NAME}_TRIBITS_DIR.
#

set(CMAKE_MODULE_PATH
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/cmake
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/utils
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/common
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/test_support
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/config_tests
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/modules
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/installation
   )

if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  message("CMAKE_MODULE_PATH='${CMAKE_MODULE_PATH}'")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsConstants.cmake")
tribits_asesrt_minimum_cmake_version()
include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsCMakePolicies.cmake"  NO_POLICY_SCOPE)

# TriBITS package_arch includes
include(TribitsIncludeDirectories)
include(TribitsFindPythonInterp)
include(TribitsGlobalMacros)
include(TribitsConfigureCTestCustom)
include(TribitsGenerateResourceSpecFile)
include(TribitsPackageDependencies)
include(TribitsPrintDependencyInfo)
include(TribitsPackagingSupport)
include(TribitsConfigureTiming)

# TriBITS utils includes
include(AdvancedSet)
include(AdvancedOption)
include(TimingUtils)
include(SetDefault)


#
# Defines a TriBITS project (the guts).
#

macro(tribits_project_impl)

  #
  # A) Basic top-level TriBITS project stuff
  #

  message("")
  message("Configuring ${PROJECT_NAME} build directory")
  message("")

  # A.1) Set some basic system vars and info you can't change
  tribits_assert_and_setup_project_and_static_system_vars()

  # A.2) Read user provided options from specified files.  It is important to
  # process these files *very* early on so that they have the same basic
  # effect of setting these variables directly in the cache.
  tribits_read_in_options_from_file()

  # A.3) Get some other basic system info that is useful early on in
  # configuration
  tribits_setup_basic_system_vars()
  tribits_find_python_interp()
  find_package(Git)
  if (NOT "${GIT_VERSION_STRING_OVERRIDE}" STREQUAL "")
    set(GIT_VERSION_STRING ${GIT_VERSION_STRING_OVERRIDE}) # For testing!
  endif()

  #
  # A.4) Read in the Project's version file
  #
  # NOTE: The file Version.cmake must be read *before* the global options are
  # read because the variables defined in Version.cmake provide defaults for
  # many of these options.
  #
  tribits_project_read_version_file(${PROJECT_SOURCE_DIR})

  # Since the version header file is now configured the root build
  # dir needs to be on the include path
  tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  #
  # B) Set up user options and global variables that will be used throughout
  #

  message("")
  message("Setting up major user options ...")
  message("")

  tribits_define_global_options_and_define_extra_repos()

  # Have to start timing after we read in the major options since that
  # determines the timing option var.
  tribits_config_code_start_timer(GLOBAL_TIME_START_SECONDS)

  tribits_read_in_native_repositories()

  tribits_combine_native_and_extra_repos()

  tribits_process_extra_repos_options_files()

  include(TribitsInstallationTestingMacros)
  tribits_find_project_install()

  #
  # C) Generate version info file and read in ${PROJECT_NAME} packages and
  # TPLs and process dependencies
  #

  tribits_generate_repo_version_output_and_file_and_install()

  # Read in and process all of the project's package, TPL, listss and
  # dependency definition files.
  tribits_read_all_project_deps_files_create_deps_graph()
  tribits_print_initial_dependency_info()
  tribits_write_xml_dependency_files_if_supported()

  #
  # D) Apply dependency logic to enable and disable TriBITS packages and tests
  #

  tribits_adjust_and_print_package_dependencies()

  #
  # E) Stop after all dependencies handling is finished if asked
  #

  if (${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
    message("")
    message("Shortcircuiting after dependency tracking ...")
    return()
  endif()

  #
  # F) Set up the environment on this computer
  #

  message("")
  message("Probing the environment ...")
  message("")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    tribits_setup_env()
  else()
    message("-- Skipping env setup due to"
      " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
  endif()

  #
  # G) Go get the information for all enabled TPLS
  #

  tribits_process_enabled_tpls()

  #
  # H) Set up for testing with CTest and ${PROJECT_NAME} test harness
  #

  message("")
  message("Setting up testing support ...")
  message("")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    tribits_include_ctest_support()
  else()
    message("-- Skipping testing support setup due to"
      " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
  endif()

  #
  # I) Add the 'dashboard' target
  #
  # NOTE: Must come after setting up for testing
  #

  set(TRIBITS_ADD_DASHBOARD_TARGET_MODULE
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CTEST_DRIVER_DIR}/TribitsAddDashboardTarget.cmake
    )
  if (
    EXISTS ${TRIBITS_ADD_DASHBOARD_TARGET_MODULE}
    AND
    NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY
    )
    include(${TRIBITS_ADD_DASHBOARD_TARGET_MODULE})
    tribits_add_dashboard_target()
  endif()

  #
  # J) Configure individual packages
  #

  message("")
  message("Configuring individual enabled ${PROJECT_NAME} packages ...")
  message("")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    tribits_repository_configure_all_version_header_files(
      ${${PROJECT_NAME}_ALL_REPOSITORIES})
    tribits_repository_configure_all_version_date_files(
      ${${PROJECT_NAME}_ALL_REPOSITORIES})
  endif()

  tribits_configure_enabled_packages()

  #
  # K) Write dummy makefiles for Ninja
  #

  if (CMAKE_GENERATOR STREQUAL "Ninja")

    if (${PROJECT_NAME}_WRITE_NINJA_MAKEFILES)

      message("")
      message("Generating dummy makefiles in each directory to call Ninja ...")
      message("")
  
      tribits_config_code_start_timer(NINJA_MAKEFILES_TIME_START_SECONDS)
  
      include(GenerateNinjaMakefiles)
      generate_ninja_makefiles(${CMAKE_SOURCE_DIR})
  
      tribits_config_code_stop_timer(NINJA_MAKEFILES_TIME_START_SECONDS
         "Total time generate Ninja makefiles ${PROJECT_NAME}")

    else()

      message("\nNOTE: *NOT* generating dummy Ninja makefiles (see above note"
        " and check CMake version)")

    endif()


  endif()

  #
  # L) Setup for packaging and distribution
  #

  if (${PROJECT_NAME}_ENABLE_CPACK_PACKAGING)
    message("")
    message("Set up for creating a distribution ...")
    message("")
    if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
      tribits_setup_packaging_and_distribution()
    else()
      message("-- Skipping distribution setup due to"
        " ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY=ON")
    endif()
  else()
    message("")
    message("Skipping setup for distribution because"
      " ${PROJECT_NAME}_ENABLE_CPACK_PACKAGING=OFF")
    message("")
  endif()

  #
  # M) Set up for installation
  #

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    tribits_setup_for_installation()
  endif()

  #
  # N) Generate resource spec file if applicable
  #

  tribits_generate_ctest_resource_spec_file_project_logic()

  #
  # O) Show final timing and end
  #

  message("")
  message("Finished configuring ${PROJECT_NAME}!")
  message("")
  tribits_config_code_stop_timer(GLOBAL_TIME_START_SECONDS
    "Total time to configure ${PROJECT_NAME}")

endmacro()


# @MACRO: tribits_project_enable_all()
#
# Process a project where you enable all of the packages by default.
#
# Usage::
#
#   tribits_project_enable_all()
#
# This macro just sets the global cache var
# `${PROJECT_NAME}_ENABLE_ALL_PACKAGES`_ to ``ON`` by default then calls
# `tribits_project()`_.  That is all.  This macro is generally used for
# TriBITS projects that have just a single package or by default just want to
# enable all packages.  This is especially useful when you have a TriBITS
# project with just a single package.
#
macro(tribits_project_enable_all)
  set(${PROJECT_NAME}_ENABLE_ALL_PACKAGES ON CACHE BOOL "Enable all by default" )
  tribits_project_impl(${ARGN})
endmacro()
