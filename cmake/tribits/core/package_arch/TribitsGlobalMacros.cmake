# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Standard TriBITS system includes

include("${CMAKE_CURRENT_LIST_DIR}/../utils/TribitsGitRepoVersionInfo.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsConstants.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/../test_support/TribitsTestCategories.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../test_support/TribitsAddTestHelpers.cmake")

include(TribitsSetupMPI)
include(TribitsGeneralMacros)
include(TribitsVerbosePrintVar)
include(TribitsProcessEnabledTpls)
include(TribitsInstallHeaders)
include(TribitsGetVersionDate)
include(TribitsReportInvalidTribitsUsage)
include(TribitsReadAllProjectDepsFilesCreateDepsGraph)
include(TribitsAdjustPackageEnables)
include(TribitsSetUpEnabledOnlyDependencies)
include(TribitsConfigureTiming)

# Standard TriBITS utilities includes
include(TribitsAddOptionAndDefine)
include(TribitsAddEnumCacheVar)
include(AdvancedOption)
include(AdvancedSet)
include(AppendStringVar)
include(AppendStringVarWithSep)
include(AssertAndTouchDefined)
include(CMakeBuildTypesList)
include(FindListElement)
include(GlobalNullSet)
include(PrintNonemptyVar)
include(PrintVar)
include(RemoveGlobalDuplicates)
include(Split)
include(TimingUtils)
include(SetDefaultAndFromEnv) # Used by some call-back files
include(TribitsFilepathHelpers)
include(TribitsDeprecatedHelpers)

# Standard CMake includes
include(CheckIncludeFileCXX)

# Include here so it does not need to be included in each individual
# FindTPL<TPLNAME>.cmake file over and over.
include(TribitsTplFindIncludeDirsAndLibraries)
include(TribitsTplDeclareLibraries) # Deprecated
# ABOVE: We need to include TribitsTplDeclareLibraries.cmake until all client
# projects stop using it.


# Assert and setup project binary directory and other project variables
#
macro(tribits_assert_and_setup_project_and_static_system_vars)

  string(APPEND IN_SOURCE_ERROR_COMMON_MSG
    "\nYou must now run something like:\n"
    "  $ cd ${CMAKE_CURRENT_SOURCE_DIR}/\n"
    "  $ rm -r CMakeCache.txt CMakeFiles/"
    "\n"
    "Please create a different directory and configure ${PROJECT_NAME}"
    " under that such as:\n"
    "  $ cd ${CMAKE_CURRENT_SOURCE_DIR}/\n"
    "  $ mkdir MY_BUILD\n"
    "  $ cd MY_BUILD\n"
    "  $ cmake [OPTIONS] .."
    )

  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt")
    message(FATAL_ERROR "ERROR! "
      "The file ${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt exists from a"
      " likely prior attempt to do an in-source build."
      "${IN_SOURCE_ERROR_COMMON_MSG}"
      )
  endif()

  if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(FATAL_ERROR "ERROR! "
      "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
      " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
      "\n${PROJECT_NAME} does not support in source builds!\n"
      "NOTE: You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under"
      " the source directory for ${PROJECT_NAME} or you will not be able to configure ${PROJECT_NAME} correctly!"
      "${IN_SOURCE_ERROR_COMMON_MSG}"
      )
  endif()

  string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)
  set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
  set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
  print_var(PROJECT_SOURCE_DIR)
  print_var(PROJECT_BINARY_DIR)
  print_var(${PROJECT_NAME}_TRIBITS_DIR)
  print_var(TriBITS_VERSION_STRING)
  # Above, we put these in the cache so we can grep them out of the cache file

  #
  # Print some basic static info provided by CMake automatically
  #

  print_var(CMAKE_VERSION)
  print_var(CMAKE_GENERATOR)

endmacro()


# Set up some really basic system variables
#
# This macro needs to be called *before* the user *.cmake option files are
# read in so that there is an opportunity to override these.
#
macro(tribits_setup_basic_system_vars)

  # CMAKE_HOST_SYSTEM_NAME is provided by CMake automatically but can actually
  # be overridden in the cache.
  print_var(CMAKE_HOST_SYSTEM_NAME)

  site_name(${PROJECT_NAME}_HOSTNAME)
  mark_as_advanced(${PROJECT_NAME}_HOSTNAME)
  print_var(${PROJECT_NAME}_HOSTNAME)

  # NOTE: CMAKE_HOST_SYSTEM_NAME and ${PROJECT_NAME}_HOSTNAME are used by
  # TRIBITS_ADD[_ADVANCED]_test() to include/exclude tests based in the
  # arguments HOSTS, XHOSTS, HOSTTYPES, AND XHOSTTYPES.

endmacro()


# Define an option to include a file that reads in a bunch of options and read
# those files
#
macro(tribits_read_in_options_from_file)

  set( ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE "" CACHE FILEPATH
    "Name of an optional file that is included first to define any cmake options with set( ... CACHE ...) calls.  NOTE: paths can be separated by commas instead of semicolons but paths cannot contain commas."
    )

  split("${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE}"  ","
    ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE)

  set( ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_ALL
    ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE}
    ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND})

  foreach (CONFIG_OPTS_FILE ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_ALL})
    message("-- " "Reading in configuration options from ${CONFIG_OPTS_FILE} ...")
    tribits_trace_file_processing(PROJECT  INCLUDE  "${CONFIG_OPTS_FILE}")
    include(${CONFIG_OPTS_FILE})
  endforeach()


endmacro()


# Assert ${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR set
# correctly
#
function(assert_project_set_group_and_permissions_on_install_base_dir)

  if (
      (NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "")
       AND
      (NOT "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}" STREQUAL "")
    )
    tribits_dir_is_basedir(
      "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}"
      "${CMAKE_INSTALL_PREFIX}"
      isBaseDir)
    if (NOT isBaseDir)
      message(FATAL_ERROR
        "\n"
        "***\n"
        "*** ERROR in ${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR!\n"
        "***\n"
        "\n"
        "${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR=${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}\n"
        "\n"
        "is not a strict base dir of:\n"
        "\n"
        "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}\n"
        "\n"
        "Either remove ${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR from the cache or set it to be a base dir of CMAKE_INSTALL_PREFIX!\n"
        "\n"
        )
    endif()
  endif()

endfunction()


# Define all of the standard global package architecture options.
#
macro(tribits_define_global_options_and_define_extra_repos)

  set( ${PROJECT_NAME}_ENABLE_ALL_PACKAGES OFF CACHE BOOL
    "Enable all packages PT packages (ST packages as well if ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE is true)." )

  set(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON CACHE BOOL
    "Recursively enable all optional packages for set of enabled packages." )

  set( ${PROJECT_NAME}_INSTALL_EXECUTABLES ON CACHE BOOL
    "Enable the installation of executables provided by the ${PROJECT_NAME} packages." )

  advanced_set(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL
    "Recursively enable all packages that have required or optional dependencies for set of enabled packages." )

  if (${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT STREQUAL "")
    set(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT OFF)
  endif()
  set(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES
    ${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT}
    CACHE BOOL
    "Disable (and printing warning) for enabled packages that have hard-disabled upstream dependencies.  Otherwise, is to raises a fatal configure failure." )

  set_cache_on_off_empty( ${PROJECT_NAME}_ENABLE_TESTS ""
    "Enable tests (execs and ctest add_test()) in all packages  (set to ON, OFF, or leave empty)." )

  set_cache_on_off_empty(${PROJECT_NAME}_ENABLE_EXAMPLES ""
    "Enable examples (exec and ctest add_test()) in all packages  (set to ON, OFF, or leave empty).  If left empty, then this will be set to ON if ${PROJECT_NAME}_ENABLE_TESTS=ON" )

  set(${PROJECT_NAME}_SKIP_CTEST_ADD_TEST OFF CACHE BOOL
    "Skip ctest add_test() for all defined tests (but still build any enabled test or example executable targets)." )

  if (${PROJECT_NAME}_ENABLE_TESTS AND ${PROJECT_NAME}_ENABLE_EXAMPLES STREQUAL "")
    message(STATUS "Setting ${PROJECT_NAME}_ENABLE_EXAMPLES=ON because ${PROJECT_NAME}_ENABLE_TESTS=ON")
    set(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  endif()

  advanced_set( ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
    "Set to empty all package enables (set to OFF at end)." )

  advanced_option(${PROJECT_NAME}_REMOVE_DEFAULT_PACKAGE_DISABLES
    "Removes all default disables from the packages list.  Used for testing etc."
    OFF )

  if ("${${PROJECT_NAME}_ENABLE_C_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_C_DEFAULT ON)
  endif()
  advanced_option(${PROJECT_NAME}_ENABLE_C
    "Enable the C compiler and related code"
    ${${PROJECT_NAME}_ENABLE_C_DEFAULT} )

  if ("${${PROJECT_NAME}_C_Standard_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_C_Standard_DEFAULT c99)
  endif()
  advanced_set(${PROJECT_NAME}_C_Standard
    ${${PROJECT_NAME}_C_Standard_DEFAULT}
    CACHE STRING
    "The standard <cstd> to use in --std=<cstd> for GCC compilers." )

  if ("${${PROJECT_NAME}_ENABLE_CXX_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_CXX_DEFAULT ON)
  endif()
  advanced_option(${PROJECT_NAME}_ENABLE_CXX
    "Enable the C++ compiler and related code"
    ${${PROJECT_NAME}_ENABLE_CXX_DEFAULT} )

  # Hard-code a variable with the same name as a now-removed option that is always enabled.
  # This can be removed after clients have been updated.
  set(${PROJECT_NAME}_ENABLE_CXX11 ON)

  if ("${${PROJECT_NAME}_ENABLE_Fortran_DEFAULT}" STREQUAL "")
    if (WIN32)
      set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
    else()
      set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
    endif()
  endif()

  option(${PROJECT_NAME}_ENABLE_Fortran
    "Enable the Fortran compiler and related code"
    ${${PROJECT_NAME}_ENABLE_Fortran_DEFAULT} )

  advanced_option(${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST
    "Skip the Fortran/C++ compatibility test"
    OFF )

  advanced_set(
    ${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE
    "${${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE_DEFAULT}"
    CACHE BOOL
    "If TRUE, the directory and file permissions on the installed directories and files will be set to world readable.  NOTE: Empty '' (the default) leaves default CMake permissions in place."
    )

  advanced_set(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE ""
    CACHE BOOL
    "If TRUE, the directory and file permissions on the installed directories and files will be set to group-writable.  Setting to TRUE also implies ${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE=TRUE.  NOTE: Empty '' (the default) avoid adding the group write permission."
    )

  if ("${${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE_DEFAULT}" STREQUAL "")
    if (${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE)
      set(${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE_DEFAULT TRUE)
    else()
      set(${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE_DEFAULT
        "${${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE}")
    endif()
  endif()
  advanced_set(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE
    "${${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE_DEFAULT}"
    CACHE BOOL
    "If TRUE, the directory and file permissions on the installed directories and files will be set to group readable.  Setting ${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE=ON implies this is 'ON' as well.  NOTE: Empty '' (the default) leaves default CMake permissions in place."
    )

  advanced_set(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP ""
    CACHE STRING
    "If set, then the installed files and directories from ${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR on down will be given over to this owning group.  The default is empty '' which means the default group will not be changed."
    )

  advanced_set(
    ${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR
    "${CMAKE_INSTALL_PREFIX}"
    CACHE FILEPATH
    "Set the base path for which a recursive chmod and chgrp will be called to set the group and permissions after the install is complete.  The default directory is give by CMAKE_INSTALL_PREFIX."
    )

  assert_project_set_group_and_permissions_on_install_base_dir()

  if ("${${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT TRUE)
  else()
    set(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT FALSE)
  endif()
  advanced_set(
    ${PROJECT_NAME}_SET_INSTALL_RPATH ${${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT}
    CACHE BOOL
    "If TRUE, then set RPATH on installed binaries will set to ${PROJECT_NAME}_INSTALL_LIB_DIR automatically"
    )

  if ("${CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT}" STREQUAL "")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT TRUE)
  else()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT FALSE)
  endif()
  advanced_set(
    CMAKE_INSTALL_RPATH_USE_LINK_PATH ${CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then the RPATH for external shared libs will be embedded in installed libs and execs."
    )

  advanced_set(${PROJECT_NAME}_EXTRA_LINK_FLAGS ""
    CACHE STRING
    "Extra flags added to the end of every linked executable"
    )

  # OpenMP is similar to a TPL in some respects, but requires only compiler
  # flags to enable

  option(${PROJECT_NAME}_ENABLE_OpenMP
    "Build with OpenMP support." OFF)

  if (CMAKE_GENERATOR STREQUAL "Ninja")
    if("${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT}" STREQUAL "")
      set(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT ON)
    endif()
    set(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES
      ${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT} CACHE BOOL
      "Generate dummy makefiles to call ninja in every build subdirectory (requires CMake 3.7.0 or newer)." )
  endif()
  if ("${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES}" STREQUAL "")
    set(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES OFF)
  endif()

  advanced_set(${PROJECT_NAME}_PARALLEL_COMPILE_JOBS_LIMIT "" CACHE STRING
    "If not empty '', gives an integer for the max number of object compile jobs for Ninja builds. (Default empty for no limit)")

  advanced_set(${PROJECT_NAME}_PARALLEL_LINK_JOBS_LIMIT "" CACHE STRING
    "If not empty '', gives an integer for the max number of lib and exe link jobs for Ninja builds. (Default empty for no limit)")
  
  if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    set(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT ON)
  else()
    set(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT OFF)
  endif()
  set(${PROJECT_NAME}_ENABLE_DEBUG ${${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT} CACHE BOOL
    "Enable debug checking for ${PROJECT_NAME} packages.  Off by default unless CMAKE_BUILD_TYPE=\"DEBUG\"." )

  if (${PROJECT_NAME}_ENABLE_DEBUG)
    set(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT ON)
  else()
    set(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT OFF)
  endif()
  set(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG
    ${${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT} CACHE BOOL
    "Enable debug checking of the process which finds errors in the project's CMake files (off by default unless ${PROJECT_NAME}_ENABLE_DEBUG=ON)." )

  if ("${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT "WARNING")
  endif()
  advanced_set(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS
    ${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT}
    CACHE STRING
    "Determines how unparsed arguments for TriBITS functions that use cmake_parase_aruments() internally are handled.  Valid choices are 'WARNING', 'SEND_ERROR', and 'FATAL_ERROR'.  The default is 'SEND_ERROR'."
    )
  if (
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "WARNING")
     AND
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "SEND_ERROR")
     AND
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "FATAL_ERROR")
    )
    message(FATAL_ERROR "Error, the value of"
      " ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS ="
      " '${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS}' is invalid!"
      " Valid values include 'WARNING', 'SEND_ERROR', and 'FATAL_ERROR'"
      )
  endif()

  set(${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR ON
    CACHE BOOL
    "Enable support for Teuchos Time Monitors in all Trilinos packages that support it."
    )

  advanced_set(${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS ON
    CACHE BOOL
    "Show warnings about deprecated code"
    )

  advanced_set(${PROJECT_NAME}_HIDE_DEPRECATED_CODE OFF
    CACHE BOOL
    "Show warnings about deprecated code"
    )

  advanced_set(${PROJECT_NAME}_VERBOSE_CONFIGURE OFF
    CACHE BOOL
    "Make the ${PROJECT_NAME} configure process verbose."
    )

  advanced_option(${PROJECT_NAME}_DUMP_LINK_LIBS
    "Dump the link libraries for every library and executable created."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  advanced_set(${PROJECT_NAME}_TRACE_FILE_PROCESSING
    ${${PROJECT_NAME}_VERBOSE_CONFIGURE}
    CACHE BOOL
    "Print out when all of the various files get processed."
    )

  advanced_set(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION ON
    CACHE BOOL
    "Enable explicit template instantiation in all packages that support it"
    )

  advanced_set(${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE OFF
    CACHE BOOL
    "Auto-generate a resource spec file for use with CTest."
    )

  advanced_set(${PROJECT_NAME}_CUDA_NUM_GPUS 1
    CACHE STRING
    "Number of GPUS to make available in the auto-generated resource spec file."
    )

  advanced_set(${PROJECT_NAME}_CUDA_SLOTS_PER_GPU 1
    CACHE STRING
    "Number of slots per GPU in the auto-generated resource spec file."
    )

  set(CTEST_RESOURCE_SPEC_FILE_DOC_EXTRA "")
  if (${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE)
    set(CTEST_RESOURCE_SPEC_FILE_DEFAULT  ${CMAKE_BINARY_DIR}/ctest_resources.json)
    if ("${CTEST_RESOURCE_SPEC_FILE}" STREQUAL "")
      set(CTEST_RESOURCE_SPEC_FILE_DOC_EXTRA
         "  This file is autogenerated by default since ${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE=${${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE}!" )
    endif()
  else()
    set(CTEST_RESOURCE_SPEC_FILE_DEFAULT "")
  endif()

  advanced_set(CTEST_RESOURCE_SPEC_FILE
    "${CTEST_RESOURCE_SPEC_FILE_DEFAULT}"
    CACHE FILEPATH
    "Resource spec file for CTest.${CTEST_RESOURCE_SPEC_FILE_DOC_EXTRA}"
    )

  if (USE_XSDK_DEFAULTS)
    # Need to set BUILD_SHARED_LIBS default here based on USE_XSDK_DEFAULTS
    # and not in tribits_setup_env() in case there is logic in TriBITS or
    # project-specific files that depends on this var getting set here!
    set(BUILD_SHARED_LIBS_DEFAULT  TRUE)
    if ("${BUILD_SHARED_LIBS}" STREQUAL "")
      message("-- " "XSDK: Setting default BUILD_SHARED_LIBS=TRUE")
    endif()
  else()
    set(BUILD_SHARED_LIBS_DEFAULT  FALSE)
  endif()
  set(BUILD_SHARED_LIBS  ${BUILD_SHARED_LIBS_DEFAULT}
    CACHE  BOOL   "Build shared libraries or not.")

  if ("${${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT  FALSE)
  endif()
  advanced_set(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS
    ${${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT}
    CACHE BOOL
    "If set TRUE, then 'SYSTEM' will be passed into include_directories() for TPL includes.")

  if ("${${PROJECT_NAME}_IMPORTED_NO_SYSTEM_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_IMPORTED_NO_SYSTEM_DEFAULT FALSE)
  endif()
  advanced_set(${PROJECT_NAME}_IMPORTED_NO_SYSTEM
    ${${PROJECT_NAME}_IMPORTED_NO_SYSTEM_DEFAULT}
    CACHE BOOL
    "If set TRUE, then set IMPORTED_NO_SYSTEM property on all exported libraries.")

  if (CMAKE_VERSION VERSION_LESS 3.23 AND ${PROJECT_NAME}_IMPORTED_NO_SYSTEM)
    message(FATAL_ERROR "Error, setting ${PROJECT_NAME}_IMPORTED_NO_SYSTEM='${${PROJECT_NAME}_IMPORTED_NO_SYSTEM}' for CMake version '${CMAKE_VERSION}' < 3.23 is not allowed!")
  endif()

  advanced_set(TPL_FIND_SHARED_LIBS ON CACHE BOOL
    "If ON, then the TPL system will find shared libs if they exist, otherwise will only find static libs." )

  advanced_set(${PROJECT_NAME}_LINK_SEARCH_START_STATIC OFF CACHE BOOL
    "If ON, then the property LINK_SEARCH_START_STATIC will be added to all executables." )

  advanced_set(${PROJECT_NAME}_LIBRARY_NAME_PREFIX ""
    CACHE STRING
    "Prefix for all ${PROJECT_NAME} library names. If set to, for example, 'prefix_',
    libraries will be named and installed as 'prefix_<libname>.*'.  Default is '' (no prefix)."
    )

  if ("${${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT FALSE)
  endif()
  advanced_set( ${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS
    ${${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then all of the TPL libs must be found for every enabled TPL."
    )

  if ("${${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT FALSE)  # Maintain backward compatibility
  endif()
  advanced_set( ${PROJECT_NAME}_USE_GNUINSTALLDIRS
    ${${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then CMake GNUInstallDris modules is used to pick standard install paths by default."
    )

  if ("${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}" STREQUAL "")
    # Assume the TriBITS project wants to install headers and libraries by default
    set(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT ON)
  endif()

  advanced_set(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS
    ${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}
    CACHE BOOL
    "Install libraries and headers (default is ${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}).  NOTE: Shared libraries are always installed since they are needed by executables."
    )

  # Creating <Package>Config.cmake files is currently *very* expensive for large
  # TriBITS projects so we disable this by default for TriBITS.
  if ("${${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT OFF)
  endif()
  advanced_set(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
    ${${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT}
    CACHE BOOL
    "Determines if ${PROJECT_NAME}Config.cmake and <PACKAGE>Config.cmake files are created or not."
    )

  if ("${${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES_DEFAULT OFF)
  endif()
  advanced_set(${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES
    ${${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES_DEFAULT}
    CACHE BOOL
    "Skip installing the file ${PROJECT_NAME}Config.cmake."
    )

  if (NOT ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT)
    # We need to generate the dependency logic for export dependency files if
    # asked.
    if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
      set(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)
    else()
      set(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT OFF)
    endif()
  endif()
  advanced_set(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES
     ${${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT} CACHE BOOL
    "Generate packages dependency data-structures needed for dependency export files." )

  # ${PROJECT_NAME}_ELEVATE_SS_TO_PS is depreciated!
  if (${PROJECT_NAME}_ELEVATE_SS_TO_PS_DEFAULT)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "WARNING: ${PROJECT_NAME}_ELEVATE_SS_TO_PS_DEFAULT is deprecated."
        "  Use ${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT instead!")
    endif()
    set(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT ON)
  endif()

  if ("${${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT OFF)
  endif()
  advanced_set( ${PROJECT_NAME}_ELEVATE_ST_TO_PT
    ${${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT}
    CACHE BOOL
    "Elevate all defined ST packages to PT packages." )

  if ("${${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT OFF)
  endif()
  advanced_set( ${PROJECT_NAME}_ENABLE_CPACK_PACKAGING
     ${${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT}
     CACHE BOOL
    "Enable support for creating a distribution using CPack" )

  if ("${${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT TRUE)
  endif()
  advanced_set( ${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION
    ${${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT}
    CACHE BOOL
    "Excluded disabled packages from the CPack-generated distribution.")

  if ("${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT OFF)
  endif()
  advanced_set(
    ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
    ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT}
    CACHE BOOL
    "Allow Secondary Tested (ST) packages and code to be implicitly enabled." )

  if ("${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT NIGHTLY)
  endif()
  advanced_set(${PROJECT_NAME}_TEST_CATEGORIES
     ${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}
     CACHE STRING
    "List of categories of tests to enable: '${${PROJECT_NAME}_VALID_CATEGORIES_STR}' (default `${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}`)."
    )
  tribits_get_invalid_categories(${PROJECT_NAME}_TEST_CATEGORIES)

  if (NOT GIT_EXECUTABLE)
    set(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT OFF)
  elseif ("${${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT}" STREQUAL "" )
    set(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT OFF)
  endif()
  advanced_set(
    ${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE
    ${${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT}
    CACHE BOOL
    "Generate the ${PROJECT_NAME}RepoVersion.txt file.")

  tribits_advanced_set_cache_var_and_default(${PROJECT_NAME}_SHOW_GIT_COMMIT_PARENTS
    BOOL OFF "Show parents' commit info in the repo version output.")

  if ("${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES_DEFAULT OFF)
  endif()
  advanced_set(
    ${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES
    ${${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES_DEFAULT}
    CACHE BOOL
    "Generate VersionDate.cmake and <RepoName>_version_date.h files.")

  if ("${DART_TESTING_TIMEOUT_DEFAULT}"  STREQUAL "")
    set(DART_TESTING_TIMEOUT_DEFAULT  1500)
  endif()
  advanced_set(
    DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT_DEFAULT}
    CACHE STRING
    "Raw CMake/CTest global default test timeout (default 1500).  (NOTE: Does not impact timeouts of tests that have the TIMEOUT property set on a test-by-test basis.)"
    )
  # NOTE: 1500 is the CMake default set in Modules/CTest.cmake.  We need to
  # set the default here because we need to be able to scale it correctly in
  # case the user does not explicitly set this var in the cache.

  advanced_set(${PROJECT_NAME}_SCALE_TEST_TIMEOUT 1.0 CACHE STRING
    "Scale factor for global DART_TESTING_TIMEOUT and individual test TIMEOUT (default 1.0)."
    )
  # NOTE: This value is 1.0, *NOT* 1!  This is used in tribits_scale_timeout()
  # and there are unit tests that rely on this default!

  advanced_set(${PROJECT_NAME}_REL_CPU_SPEED 1.0 CACHE STRING
    "Relative CPU speed of the computer used to scale performance tests (default 1.0)."
    )

  if ("${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT ON)
  endif()
  advanced_set( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}
    CACHE BOOL
    "Determines if a variety of development mode checks are turned on by default or not."
    )

  if (NOT "${${PROJECT_NAME}_ASSERT_MISSING_PACKAGES}" STREQUAL "")
    tribits_deprecated("Warning, ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES="
      "'${${PROJECT_NAME}_ASSERT_MISSING_PACKAGES}' is set and is no"
      " longer supported!  Please set"
      " ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES instead (see build ref)!" )
    if (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
      set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT  FATAL_ERROR)
    else()
     set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT  IGNORE)
    endif()
  endif()

  if ("${${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT}" STREQUAL "")
    if (${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
      set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT  FATAL_ERROR)
    else()
      set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT  IGNORE)
    endif()
  endif()

  set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST
    "FATAL_ERROR" "SEND_ERROR" )
  set(${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_VALUES_LIST
    ${${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_ERROR_VALUES_LIST}
    "WARNING" "NOTICE" "IGNORE" "OFF" )
  tribits_add_enum_cache_var( ${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES
    DEFAULT_VAL ${${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_DEFAULT}
    DOC_STRING
      "Assert that all external and internal dependencies are defined in the project"
    ALLOWED_STRINGS_LIST
      ${${PROJECT_NAME}_ASSERT_DEFINED_DEPENDENCIES_VALUES_LIST}
    IS_ADVANCED )

  if ("${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_DEFAULT}" STREQUAL "")
    if (${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
      set(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_DEFAULT  FATAL_ERROR)
    else()
      set(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_DEFAULT  IGNORE)
    endif()
  endif()
  set(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_VALUES_LIST
      "FATAL_ERROR" "SEND_ERROR" "WARNING" "IGNORE" "OFF")
  tribits_add_enum_cache_var( ${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE
    DEFAULT_VAL "${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_DEFAULT}"
    DOC_STRING
      "Assert correct usage of TriBITS"
    ALLOWED_STRINGS_LIST
      ${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE_VALUES_LIST}
    IS_ADVANCED )

  advanced_set( ${PROJECT_NAME}_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES
    FALSE  CACHE  BOOL
    "If set to TRUE, a 'NOTE' is printed for each missing package that is ignored." )

  advanced_set( ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C code for supported compilers." )

  advanced_set( ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C++ code for supported compilers." )

  multiline_set( ENABLE_SHADOW_WARNINGS_DOC
    "Turn ON or OFF shadowing warnings for all packages where strong warnings have"
    " not been explicitly disabled.  Setting the empty '' let's each package decide." )
  set_cache_on_off_empty( ${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS ""
    "${ENABLE_SHADOW_WARNINGS_DOC}" )
  mark_as_advanced(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)

  advanced_set( ${PROJECT_NAME}_ENABLE_COVERAGE_TESTING OFF
    CACHE BOOL "Enable support for coverage testing by setting needed compiler/linker options." )

  advanced_set( ${PROJECT_NAME}_ENABLE_CHECKED_STL OFF
    CACHE BOOL "Turn on checked STL checking (e.g. -D_GLIBCXX_DEBUG) or not." )

  advanced_set( ${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS OFF
    CACHE BOOL "Turn on debugging symbols (e.g. -g) or not if not a full debug build." )

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "-Werror")
  else()
    set(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "")
  endif()

  advanced_set( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT}"
    CACHE STRING "Flags for treating warnings as errors (for all compilers, -Werror by default for GNU).  To turn off warnings as errors set to ''")

  advanced_set(${PROJECT_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF CACHE BOOL
    "If test output complaining about circular references is found, then the test will fail." )

  advanced_set(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR ""
    CACHE FILEPATH
    "If set to non-null, this is the default directory where package dependency files will be written.")

  if (${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR)
    set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}")
  else()
    set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT "")
  endif()
  advanced_set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
    "${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "Output XML file containing ${PROJECT_NAME} dependenices used by tools (if not empty)." )

  if(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR AND
    ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE
    )
    set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
  else()
    set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT "")
  endif()
  advanced_set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
    "${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "Output XML file used by CDash in ${PROJECT_NAME}-independent format (if not empty)." )

  if(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR AND
    ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE
    )
    set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME}" )
  else()
    set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT "")
  endif()
  advanced_set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
    "${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "HTML ${PROJECT_NAME} dependenices file that will be written to (if not empty)." )

  #
  # Extra repositories
  #

  assert_defined(${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME)

  set(DEFAULT_EXTRA_REPOS_FILE
    "${PROJECT_SOURCE_DIR}/cmake/${${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME}")

  if (EXISTS ${DEFAULT_EXTRA_REPOS_FILE})
    #message("${DEFAULT_EXTRA_REPOS_FILE} does exist!")
    set(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT ${DEFAULT_EXTRA_REPOS_FILE})
  else()
    #message("${DEFAULT_EXTRA_REPOS_FILE} does *NOT* exist!")
    set(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT)
  endif()

  advanced_set(${PROJECT_NAME}_EXTRAREPOS_FILE
    "${${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT}"
    CACHE STRING
    "File containing the list of extra repositories containing add-on packages to process")
  #print_var(${PROJECT_NAME}_EXTRAREPOS_FILE)

  advanced_set(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
    ""
    CACHE STRING
    "Type of testing to pull in extra repositories (Continuous, or Nightly)" )

  advanced_set(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES
    FALSE CACHE BOOL
   "Set if to ignore missing extra repositories (or fail hard)" )

  # Even if a project does not support an extra repos file, it can always
  # support extra repositories defined by the user by the very nature of
  # Tribits.

  advanced_set(${PROJECT_NAME}_PRE_REPOSITORIES
    ""
    CACHE STRING
    "List of pre-extra repositories that contain extra ${PROJECT_NAME} packages."
    )
  split("${${PROJECT_NAME}_PRE_REPOSITORIES}"  "," ${PROJECT_NAME}_PRE_REPOSITORIES)

  advanced_set(${PROJECT_NAME}_EXTRA_REPOSITORIES
    ""
    CACHE STRING
    "List of post-extra repositories that contain extra ${PROJECT_NAME} packages."
    )
  split("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  "," ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  if ("${${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST}"  STREQUAL  "")
    set(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST  TRUE)
  endif()

  tribits_get_and_process_extra_repositories_lists()

  advanced_set(${PROJECT_NAME}_INSTALLATION_DIR
    ""  CACHE  STRING
    "Location of an installed version of ${PROJECT_NAME} that will be built against during installation testing"
    )

  #
  # More options
  #

  if("${${PROJECT_NAME}_INSTALLATION_DIR}" STREQUAL "")
    set(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT OFF)
  else()
    set(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT ON)
  endif()

  advanced_set(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    ${${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT}
    CACHE STRING
    "Enable testing against an installed version of ${PROJECT_NAME}."
    )

  advanced_option(${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING
    "Short-circuit after dependency handling is complete"
    OFF )

  advanced_option(${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY
    "Only trace dependency handling.  Don't configure to build anything!"
    OFF )

  advanced_set(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
    FALSE CACHE BOOL
   "Set to 'ON' to see configure times (Unix/Linux systems only)" )

  advanced_set(${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
    FALSE CACHE BOOL
   "Set to 'ON' to see configure times for individual packages" )

  if ("${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT}"
    STREQUAL ""
    )
    set(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT OFF)
  endif()
  advanced_set(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME
    ${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT}
    CACHE BOOL
    "Set to 'ON' to see start and end date/time for advanced tests." )

  if ("${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT}"
    STREQUAL ""
    )
    set(${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT OFF)
  endif()
  advanced_set(${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST
    ${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT}
    CACHE BOOL
    "Set to 'ON' to see the machine load for advanced tests." )

  if ("${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_DEFAULT}" STREQUAL "")
    set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_DEFAULT  "DEPRECATION")
  endif()

  tribits_add_enum_cache_var(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE
    DEFAULT_VAL "${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_DEFAULT}"
    DOC_STRING "Mode for dealing with usage of TriBITS deprecated functionality"
    ALLOWED_STRINGS_LIST ${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_ALL_VALID_VALUES}
    IS_ADVANCED
    )

  mark_as_advanced(BUILD_TESTING)
  mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)
  mark_as_advanced(DART_TESTING_TIMEOUT)
  mark_as_advanced(EXECUTABLE_OUTPUT_PATH)
  mark_as_advanced(LIBRARY_OUTPUT_PATH)
  mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
  mark_as_advanced(CMAKE_OSX_SYSROOT)

endmacro()


macro(tribits_setup_installation_paths)

  #
  # A) Determine if we are going to be using default paths from GNUInstallDirs module
  #

  set(TRIBITS_USE_GNUINSTALLDIRS TRUE)

  if (NOT ${PROJECT_NAME}_USE_GNUINSTALLDIRS)
    # For backward compatibility and unit testing
    set(TRIBITS_USE_GNUINSTALLDIRS FALSE)
  endif()

  #
  # B) Pick the defaults for the install dirs
  #

  if (TRIBITS_USE_GNUINSTALLDIRS)
    include(GNUInstallDirs)
    set(${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT ${CMAKE_INSTALL_INCLUDEDIR})
    set(${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT ${CMAKE_INSTALL_LIBDIR})
    set(${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT ${CMAKE_INSTALL_BINDIR})
    set(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT "example")
  else()
    set(${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT "include")
    set(${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT "lib")
    set(${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT "bin")
    set(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT "example")
  endif()

  #
  # C) Set the cache variables for the install dirs
  #

  advanced_set( ${PROJECT_NAME}_INSTALL_INCLUDE_DIR
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT}
    CACHE PATH
    "Location where the headers will be installed.  If given as a STRING type and relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'include'"
    )

  advanced_set( ${PROJECT_NAME}_INSTALL_LIB_DIR
    ${${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT}
    CACHE PATH
    "Location where the libraries will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'lib'"
    )

  advanced_set( ${PROJECT_NAME}_INSTALL_RUNTIME_DIR
    ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT}
    CACHE PATH
    "Location where the runtime DLLs and designated programs will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'bin'"
    )

  advanced_set(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR
    ${${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT}
    CACHE PATH
    "Location where assorted examples will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'example'"
    )

  #
  # D) Setup RPATH handling
  #

  print_var(${PROJECT_NAME}_SET_INSTALL_RPATH)
  print_var(CMAKE_INSTALL_RPATH_USE_LINK_PATH)

  if (${PROJECT_NAME}_SET_INSTALL_RPATH)
    if ("${CMAKE_INSTALL_RPATH}" STREQUAL "")
      message("-- " "Setting default for CMAKE_INSTALL_RPATH pointing to ${PROJECT_NAME}_INSTALL_LIB_DIR")
      assert_defined(CMAKE_INSTALL_PREFIX)
      assert_defined(${PROJECT_NAME}_INSTALL_LIB_DIR)
      if (IS_ABSOLUTE ${${PROJECT_NAME}_INSTALL_LIB_DIR})
        set(CMAKE_INSTALL_RPATH
          "${PROJECT_NAME}_INSTALL_LIB_DIR}" )
      else()
        set(CMAKE_INSTALL_RPATH
          "${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}" )
      endif()
    endif()
    if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
      if ("${CMAKE_MACOSX_RPATH}" STREQUAL "")
        message("-- " "Setting default CMAKE_MACOSX_RPATH=TRUE")
        set(CMAKE_MACOSX_RPATH TRUE)
      endif()
      print_var(CMAKE_MACOSX_RPATH)
    endif()
  endif()
  string(REPLACE ":" ";" CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}")
  print_var(CMAKE_INSTALL_RPATH)

  #
  # E) Set permissions on created installation directories
  #

  if (
    (NOT "${${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE}" STREQUAL "")
    OR
    (NOT "${${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE}" STREQUAL "")
    )

    print_var(${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE)
    print_var(${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE)
  
    # Group permissions
    if (${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE 
      OR ${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE
      )
      set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_GROUP
        GROUP_READ GROUP_EXECUTE)
    else()
      set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_GROUP) # Empty
    endif()
  
    # World permissions
    if (${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE)
      set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_WORLD
        WORLD_READ WORLD_EXECUTE)
    else()
      set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_WORLD) # Empty
    endif()
  
    # Directory permissions
    set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      ${CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_GROUP}
      ${CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_WORLD}
      )
  
    # Print the permissions in a way that allows for strong testing
    string(REPLACE ";" " " CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_W_SPACES
      "${CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS}" )
    message("-- " "CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS = "
     "(${CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS_W_SPACES})")

  endif()

endmacro()


# Read in the Project's native repositories.
#
# On output, the variable ${PRJOECT_NAME}_NATIVE_REPOSITORIES is set.
#
macro(tribits_read_in_native_repositories)
  if (${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE)
    if (IS_ABSOLUTE ${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
      set(NATIVE_REPO_FILE ${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
    else()
      set(NATIVE_REPO_FILE
        ${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
    endif()
  else()
    set(NATIVE_REPO_FILE ${PROJECT_SOURCE_DIR}/cmake/NativeRepositoriesList.cmake)
  endif()
  if (EXISTS ${NATIVE_REPO_FILE})
    tribits_trace_file_processing(PROJECT  INCLUDE  "${NATIVE_REPO_FILE}")
    include(${NATIVE_REPO_FILE})
  else()
    set(${PROJECT_NAME}_NATIVE_REPOSITORIES ".")
  endif()
endmacro()


# Combine native and extra repos lists into a single list.
#
# Combines ${PROJECT_NAME}_PRE_REPOSITORIES
# ${PROJECT_NAME}_NATIVE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES
# into a single list ${PROJECT_NAME}_ALL_REPOSITORIES.
#
macro(tribits_combine_native_and_extra_repos)
  assert_defined(${PROJECT_NAME}_PRE_REPOSITORIES)
  assert_defined(${PROJECT_NAME}_NATIVE_REPOSITORIES)
  assert_defined(${PROJECT_NAME}_EXTRA_REPOSITORIES)
  set( ${PROJECT_NAME}_ALL_REPOSITORIES
    ${${PROJECT_NAME}_PRE_REPOSITORIES}
    ${${PROJECT_NAME}_NATIVE_REPOSITORIES}
    ${${PROJECT_NAME}_EXTRA_REPOSITORIES}
    )
endmacro()


# Process extra repo extra options files
#
macro(tribits_process_extra_repos_options_files)
  # Loop through the Repositories, set their base directories and run their
  # options setup callback functions.
  foreach(REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    tribits_get_repo_name_dir(${REPO}  REPO_NAME  REPO_DIR)
    tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR}  ${REPO_DIR}  ${REPO_NAME}_SOURCE_DIR)
    tribits_set_base_repo_dir(${PROJECT_BINARY_DIR}  ${REPO_DIR}  ${REPO_NAME}_BINARY_DIR)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Processing extra options call-backs for ${REPO}")
      print_var(${REPO_NAME}_SOURCE_DIR)
      print_var(${REPO_NAME}_BINARY_DIR)
    endif()
    tribits_repository_setup_extra_options_runner(${REPO_NAME})
  endforeach()
endmacro()


# Copy an simple text file to the binary dir to be included in the tarball
#
macro(tribits_copy_installer_resource _varname _source _destination)
  set("${_varname}" "${_destination}")
  if (EXISTS "${_destination}")
    file(REMOVE_RECURSE "${_destination}")
  ENDIF ()
  configure_file(
    "${_source}"
    "${_destination}"
    COPYONLY)
endmacro()


# Get the versions of all the git repos
#
function(tribits_generate_repo_version_file_string  projectRepoVersionFileStrOut)

  set(projectRepoVersionFileStr "")

  tribits_generate_single_repo_version_string(
    ${CMAKE_CURRENT_SOURCE_DIR} singleRepoVersionStr
    INCLUDE_COMMIT_PARENTS ${${PROJECT_NAME}_SHOW_GIT_COMMIT_PARENTS})
  string(APPEND projectRepoVersionFileStr
    "*** Base Git Repo: ${PROJECT_NAME}\n"
    "${singleRepoVersionStr}\n" )

  set(extraRepoIdx 0)
  foreach(extraRepo ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

    if (${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
      # Read from an extra repo file with potentially different dir.
      list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${extraRepoIdx}
        extraRepoDir )
    else()
       # Not read from extra repo file so dir is same as name
       set(extraRepoDir ${extraRepo})
    endif()

    tribits_generate_single_repo_version_string(
       "${CMAKE_CURRENT_SOURCE_DIR}/${extraRepoDir}" singleRepoVersionStr
       INCLUDE_COMMIT_PARENTS ${${PROJECT_NAME}_SHOW_GIT_COMMIT_PARENTS})
    string(APPEND projectRepoVersionFileStr
      "*** Git Repo: ${extraRepoDir}\n"
      "${singleRepoVersionStr}\n" )

    math(EXPR extraRepoIdx "${extraRepoIdx}+1")

  endforeach()

  set(${projectRepoVersionFileStrOut} ${projectRepoVersionFileStr} PARENT_SCOPE)

endfunction()


# Generate the project repos version file and print to stdout
#
# This function is designed so that it can be unit tested from inside of a
# cmake -P script.
#
function(tribits_generate_repo_version_output_and_file)
  # Get the repos versions
  tribits_generate_repo_version_file_string(projectRepoVersionFileStr)
  # Print the versions
  message("\n${PROJECT_NAME} repos versions:\n"
    "--------------------------------------------------------------------------------\n"
    "${projectRepoVersionFileStr}"
    " --------------------------------------------------------------------------------\n"
    )
  #) Write out the version file
  file(WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
    "${projectRepoVersionFileStr}")
endfunction()


# Create project dependencies file and create install target for these
#
# NOTE: Before calling this function, the extra repos datastructure must be
# filled out!
#
# NOTE: This function cannot be called in a cmake -P script because it has a
# call to install()!  That is why this function is separated out from
# tribits_generate_repo_version_output_and_file().
#
function(tribits_generate_repo_version_output_and_file_and_install)

  #
  # A) Create the ${PROJECT_NAME}RepoVersion.txt file if requested
  #

  print_var(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE)
  if (${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE)

    # A) Make sure that there is a .git dir in the project before generating
    if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
      set(PROJECT_SOURCE_IS_GIT_REPO TRUE)
    else()
      set(PROJECT_SOURCE_IS_GIT_REPO FALSE)
    endif()
    if (PROJECT_SOURCE_IS_GIT_REPO)
      # Get repo versions, print to stdout and write file
      tribits_generate_repo_version_output_and_file()
      # Add install target for this file
      install(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
        DESTINATION "." )
    else()
      message("\nNOTE: Skipping generation of ${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
        " because project source is not a git repo!")
    endif()

    # B) Install the repo version file if it is in source tree which it will
    # be for a tarball (see tribits_setup_packaging_and_distribution()).
    set(REPO_VERSION_FILE_IN_SOURCE_TREE
      ${CMAKE_CURRENT_SOURCE_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME})
    if (EXISTS ${REPO_VERSION_FILE_IN_SOURCE_TREE})
      install(
        FILES "${REPO_VERSION_FILE_IN_SOURCE_TREE}"
        DESTINATION "." )
    endif()

  endif()

endfunction()


# Adjust package enable logic and print out before and after state
#
# On output sets:
#
# * ${PROJECT_NAME}_NUM_ENABLED_PACKAGES: Number of enabled packages (local variable)
# * ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}: Enable status of PACKAGE_NAME (local variable)
#    ToDo: Fill in others as well!
#
macro(tribits_adjust_and_print_package_dependencies)
  tribits_config_code_start_timer(ADJUST_PACKAGE_DEPS_TIME_START_SECONDS)
  tribits_print_enables_before_adjust_package_enables()
  tribits_adjust_package_enables()
  tribits_print_enables_after_adjust_package_enables()
  tribits_handle_project_extra_link_flags_as_a_tpl()
  tribits_set_up_enabled_only_dependencies()
  tribits_config_code_stop_timer(ADJUST_PACKAGE_DEPS_TIME_START_SECONDS
    "\nTotal time to adjust package and TPL enables")
endmacro()


# Tack on ${PROJECT_NAME}_EXTRA_LINK_LIBS as a TPL that every downstream
# external and internal package depends on
#
macro(tribits_handle_project_extra_link_flags_as_a_tpl)

  if (${PROJECT_NAME}_EXTRA_LINK_FLAGS)

    set(lastLibTplName ${PROJECT_NAME}TribitsLastLib)

    # Define the TPL ${PROJECT_NAME}TribitsLastLib and its find module
    set(${lastLibTplName}_FINDMOD
      "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/FindTPLProjectLastLib.cmake")

    # Tack on ${PROJECT_NAME}TribitsLastLib as a dependency to all packages
    foreach(packageName ${${PROJECT_NAME}_DEFINED_PACKAGES})
      tribits_get_package_enable_status(${packageName}  packageEnable  "")
      list(APPEND ${packageName}_LIB_DEFINED_DEPENDENCIES ${lastLibTplName})
      if (packageEnable)
        list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${lastLibTplName})
      endif()
    endforeach()

    # Prepend ${PROJECT_NAME}TribitsLastLib to the list of packages
    list(PREPEND ${PROJECT_NAME}_DEFINED_TPLS ${lastLibTplName})
    list(PREPEND ${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES ${lastLibTplName})
    list(PREPEND ${PROJECT_NAME}_DEFINED_PACKAGES ${lastLibTplName})
    set(TPL_ENABLE_${lastLibTplName} ON)
    set(${lastLibTplName}_PACKAGE_BUILD_STATUS EXTERNAL)

  endif()

endmacro()


#
# Macros for setting up the standard environment
#


macro(tribits_setup_env)

  tribits_config_code_start_timer(SETUP_ENV_TIME_START_SECONDS)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    set(TRIBITS_SETUP_ENV_DEBUG  TRUE)
  endif()

  # Apply XSDK defaults

  if ("${${PROJECT_NAME}_TRIBITS_XSDK_DIR}" STREQUAL "")
    set(${PROJECT_NAME}_TRIBITS_XSDK_DIR  "${${PROJECT_NAME}_TRIBITS_DIR}/xsdk")
  endif()
  if (EXISTS "${${PROJECT_NAME}_TRIBITS_XSDK_DIR}")
    set(USE_XSDK_DEFAULTS_DEFAULT  FALSE) # Set to TRUE for Trilinos 13.0.0?
    set(XSDK_ENABLE_C  ${${PROJECT_NAME}_ENABLE_C})
    set(XSDK_ENABLE_CXX  ${${PROJECT_NAME}_ENABLE_CXX})
    set(XSDK_ENABLE_Fortran  ${${PROJECT_NAME}_ENABLE_Fortran})
    include("${${PROJECT_NAME}_TRIBITS_XSDK_DIR}/XSDKDefaults.cmake")
    # NOTE: BUILD_SHARED_LIBS was set in
    # tribits_define_global_options_and_define_extra_repos() based on
    # USE_XSDK_DEFAULTS in case there is logic in TriBITS that depends on this
    # var getting set there.
  endif()

  # BUILD_SHARED_LIBS
  print_var(BUILD_SHARED_LIBS)

  # Set to release build by default

  if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message(STATUS "Setting CMAKE_BUILD_TYPE=RELEASE since it was not set ...")
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Type of build to perform (i.e. DEBUG, RELEASE, NONE)" )
  else()
    string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UP)
    list(FIND CMAKE_BUILD_TYPES_LIST ${CMAKE_BUILD_TYPE_UP} BUILD_TYPE_IDX)
    if (BUILD_TYPE_IDX EQUAL -1)
      message(SEND_ERROR "Error, the given CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        " is not in the list of valid values \"${CMAKE_BUILD_TYPES_LIST}\"!")
    endif()
  endif()
  print_var(CMAKE_BUILD_TYPE)

  # Override the silly CMAKE_CONFIGURATION_TYPES variable.  This is needed for
  # MSVS!  Later, we Override CMAKE_CONFIGURATION_TYPES to just one
  # configuration after the compiler checks (see below).
  if (CMAKE_CONFIGURATION_TYPES)
    if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      set(CMAKE_CONFIGURATION_TYPE "Debug")
    elseif(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
      set(CMAKE_CONFIGURATION_TYPE "Release")
    else()
      set(CMAKE_CONFIGURATION_TYPE "Release")
    endif()
  else()
    set(CMAKE_CONFIGURATION_TYPE "")
  endif()
  if (TRIBITS_SETUP_ENV_DEBUG)
    print_var(CMAKE_CONFIGURATION_TYPE)
  endif()

  # Set up MPI if MPI is being used

  if ("${TPL_ENABLE_MPI}" STREQUAL "")
    # If TPL_ENABLE_MPI is undefined or empty because this project does not
    # define an MPI TPL, then explicitly disable it.
    set(TPL_ENABLE_MPI FALSE)
  endif()

  if (TPL_ENABLE_MPI)
    tribits_setup_mpi()
  endif()

  # Enable compilers

  assert_defined(${PROJECT_NAME}_ENABLE_C)
  if (${PROJECT_NAME}_ENABLE_C)
    enable_language(C)
    include(CMakeDetermineCCompiler)
    print_var(CMAKE_C_COMPILER_ID)
    print_var(CMAKE_C_COMPILER_VERSION)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  endif()

  assert_defined(${PROJECT_NAME}_ENABLE_CXX)
  if (${PROJECT_NAME}_ENABLE_CXX)
    enable_language(CXX)
    include(CMakeDetermineCXXCompiler)
    print_var(CMAKE_CXX_COMPILER_ID)
    print_var(CMAKE_CXX_COMPILER_VERSION)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  endif()

  assert_defined(${PROJECT_NAME}_ENABLE_Fortran)
  if (${PROJECT_NAME}_ENABLE_Fortran)
    enable_language(Fortran)
  endif()

  # Do some project-specific tweaks for compiler options, etc.
  set(PROJECT_COMPILER_CONFIG_FILE
    # Can be used for things like Kokkos.
    "${${PROJECT_NAME}_SOURCE_DIR}/cmake/ProjectCompilerPostConfig.cmake"
    CACHE FILEPATH
    "Allow for project-specific compiler settings."
   )
  if (EXISTS "${PROJECT_COMPILER_CONFIG_FILE}")
    tribits_trace_file_processing(PROJECT  INCLUDE  "${PROJECT_COMPILER_CONFIG_FILE}")
    include("${PROJECT_COMPILER_CONFIG_FILE}")
  endif()

  # Set up C++ language standard selection.
  if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
  elseif (NOT CMAKE_CXX_STANDARD MATCHES "^(11|14|17|20)$")
    message(FATAL_ERROR "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} is not 11, 14, 17, or 20.")
  ENDIF ()
  set(${PROJECT_NAME}_CXX_STANDARD_FEATURE cxx_std_${CMAKE_CXX_STANDARD})
  if (NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  ENDIF ()
  if (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
  ENDIF ()

  # Set up for strong compiler warnings and warnings as errors

  include(TribitsSetupBasicCompileLinkFlags)
  tribits_setup_basic_compile_link_flags()

  #
  # The compilers are set, the environment is known to CMake.  Now set the
  # installation paths and options.
  #
  tribits_setup_installation_paths()

  # Set up Windows interface stuff

  if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_DEPRECATE
      -D_CRT_NONSTDC_NO_DEPRECATE  -D_SCL_SECURE_NO_WARNINGS)
    set(WIN_INTERFACE_INCL  ${${PROJECT_NAME}_TRIBITS_DIR}/win_interface/include)
    if (EXISTS "${WIN_INTERFACE_INCL}")
      include_directories("${WIN_INTERFACE_INCL}")
      if (TRIBITS_SETUP_ENV_DEBUG)
        message("-- Adding win_interface/include ...")
      endif()
    endif()
  endif()

  if (WIN32 AND NOT CYGWIN)
    set(NATIVE_MS_WINDOWS TRUE)
  else()
    set(NATIVE_MS_WINDOWS FALSE)
  endif()

  # Probe for non-standard headers

  if (${PROJECT_NAME}_ENABLE_CXX)
    check_include_file_cxx(sys/time.h HAVE_SYS_TIME_H)
    check_include_file_cxx(time.h HAVE_TIME_H)
    check_include_file_cxx(stdint.h HAVE_STDINT_H)
    check_include_file_cxx(inttypes.h HAVE_INTTYPES_H)
  endif()

  set(HAVE_ALGORITHM TRUE)
  set(HAVE_CASSERT TRUE)
  set(HAVE_CCTYPE TRUE)
  set(HAVE_CERRNO TRUE)
  set(HAVE_CLIMITS TRUE)
  set(HAVE_CMATH TRUE)
  set(HAVE_COMPLEX TRUE)
  set(HAVE_CSTDARG TRUE)
  set(HAVE_CSTDIO TRUE)
  set(HAVE_CSTDLIB TRUE)
  set(HAVE_CSTRING TRUE)
  set(HAVE_IOMANIP TRUE)
  set(HAVE_IOSTREAM TRUE)
  set(HAVE_ITERATOR TRUE)
  set(HAVE_LIST TRUE)
  set(HAVE_MAP TRUE)
  set(HAVE_MEMORY TRUE)
  set(HAVE_MUTABLE TRUE)
  set(HAVE_NAMESPACES TRUE)
  set(HAVE_NEW_FOR_SCOPING TRUE)
  set(HAVE_NUMERIC TRUE)
  set(HAVE_NUMERIC_LIMITS TRUE)
  set(HAVE_POW TRUE)
  set(HAVE_SET TRUE)
  set(HAVE_SSTREAM TRUE)
  set(HAVE_FSTREAM TRUE)
  set(HAVE_STDEXCEPT TRUE)
  set(HAVE_STRING TRUE)
  set(HAVE_VECTOR TRUE)

  # 2008/12/20: rabartl: Above: All of these defines should be removed
  # because we decided that we were going to assume that all compilers
  # have these C++98 standard features.  We will deal with cases where
  # this is not true but we should not assume the worst right from the
  # beginning.

  # Find Perl

  find_package(Perl)

  # Do Fortran stuff

  include(TribitsFortranMangling)

  # Get BLAS name mangling
  #
  # ToDo: Make this a project-specific specialization

  include(TribitsBLASMangling)

  # Set up some MPI info

  if (TPL_ENABLE_MPI)
    set(HAVE_MPI TRUE)
  else()
    set(HAVE_MPI FALSE)
  endif()

  # OpenMP isn't really a TPL because support is built into the compiler.
  if(${PROJECT_NAME}_ENABLE_OpenMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
      tribits_set_openmp_flags(CXX)
      tribits_set_openmp_flags(C)
      if(OpenMP_Fortran_FLAGS)
        tribits_set_openmp_flags(Fortran)
      else()
      # Older versions of FindOpenMP.cmake don't find Fortran flags.  Mike H said this is safe.
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
      endif()
    else()
      message(FATAL_ERROR "Could not find OpenMP, try setting OpenMP_C_FLAGS and OpenMP_CXX_FLAGS directly")
    endif(OPENMP_FOUND)
  endif(${PROJECT_NAME}_ENABLE_OpenMP)

  # Check if we need the math library or not and find the right one
  if (NOT NATIVE_MS_WINDOWS)
    include(MathLibraryNeeded)
  endif()

  # Check for isnan and isinf support
  if (${PROJECT_NAME}_ENABLE_CXX)
    include(FiniteValue)
  endif()

  # Check for Doxygen/dot - We can use variables set in this check to
  # enable/disable the graphical dependency graphs in doxygen Doxyfiles.
  include(FindDoxygen)

  # You have to override the configuration types for MSVS after the compiler
  # checks!
  set(CMAKE_CONFIGURATION_TYPES  ${CMAKE_CONFIGURATION_TYPE}
    CACHE STRING
    "Override by TriBITS (see TribitsDevelopersGuilde.*)"
    FORCE)
  if (CMAKE_CONFIGURATION_TYPES)
    print_var(CMAKE_CONFIGURATION_TYPES)
  endif()

  tribits_config_code_stop_timer(SETUP_ENV_TIME_START_SECONDS
    "\nTotal time to probe and setup the environment")

  # Set ninja compile and link parallel job limits

  if (${PROJECT_NAME}_PARALLEL_COMPILE_JOBS_LIMIT)
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS
      compile_job_pool=${${PROJECT_NAME}_PARALLEL_COMPILE_JOBS_LIMIT})
    set(CMAKE_JOB_POOL_COMPILE compile_job_pool)
  endif()
  
  if (${PROJECT_NAME}_PARALLEL_LINK_JOBS_LIMIT)
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS
      link_job_pool=${${PROJECT_NAME}_PARALLEL_LINK_JOBS_LIMIT})
    set(CMAKE_JOB_POOL_LINK link_job_pool)
  endif()

endmacro()


macro(tribits_set_openmp_flags  LANG)
  if (NOT "${OpenMP_${LANG}_FLAGS_OVERRIDE}" STREQUAL "")
    set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} ${OpenMP_${LANG}_FLAGS_OVERRIDE}")
  else()
    set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} ${OpenMP_${LANG}_FLAGS}")
  endif()
endmacro()


# Set mapping of labels to subprojects (i.e. TriBITS packages) for local CTest
# only
#
# NOTE: This macro is only used define mapping of labels to subprojects for
# running ctest locally.  This results in summarizing the tests run for each
# subproject (TriBITS package) if any tests were run.  Therefore, it is
# harmless to define the mapping for every TriBITS package.  Only TriBITS
# packages will be listed in the summary if they had one or more tests run.
#
macro(tribits_set_labels_to_subprojects_mapping)
  set(CTEST_LABELS_FOR_SUBPROJECTS ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
endmacro()


# Macro to turn on CTest support
#
macro(tribits_include_ctest_support)

  set(DART_TESTING_TIMEOUT_IN ${DART_TESTING_TIMEOUT})

  if (DART_TESTING_TIMEOUT_IN)
    tribits_scale_timeout(${DART_TESTING_TIMEOUT} DART_TESTING_TIMEOUT)
    if (NOT DART_TESTING_TIMEOUT STREQUAL DART_TESTING_TIMEOUT_IN)
     message("-- DART_TESTING_TIMEOUT=${DART_TESTING_TIMEOUT_IN} being scaled by ${PROJECT_NAME}_SCALE_TEST_TIMEOUT=${${PROJECT_NAME}_SCALE_TEST_TIMEOUT} to ${DART_TESTING_TIMEOUT}")
    endif()
    # Have to set DART_TESTING_TIMEOUT in cache or CMake will not put in right
    # 'TimeOut' in DartConfiguration.tcl file!
    set(DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT} CACHE STRING "" FORCE)
  endif()

  # Set up CTEst/CDash subprojects
  tribits_set_labels_to_subprojects_mapping()
  # NOTE: We do this after all of the packages have been defined but before
  # the DartConfiguration.tcl file has been created.

  include(CTest)  # Generates file DartConfiguration.tcl with 'TimeOut' set!

  if (DART_TESTING_TIMEOUT_IN)
    # Put DART_TESTING_TIMEOUT back to user input value to avoid scaling this
    # up and up on recofigures!
    set(DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT_IN} CACHE STRING
      "Original value set by user reset by TriBITS after scaling" FORCE)
  endif()

  tribits_configure_ctest_custom(${${PROJECT_NAME}_SOURCE_DIR}
    ${${PROJECT_NAME}_BINARY_DIR})

  tribits_add_test_helpers_init()

endmacro()
# NOTE: The above logic with DART_TESTING_TIMEOUT is a huge hack.  For some
# reason, on the first configure CMake will not put the local value of the
# scaled DART_TESTING_TIMEOUT variable into the DartConfiguration.tcl.
# Instead, it uses the value of DART_TESTING_TIMEOUT that is in the cache.
# But on reconfigures, CMake uses the value of the local variable
# DART_TESTING_TIMEOUT and ignores the value in the cache (very irritating).
# Therefore, to get CMake to put in the right value for 'TimeOut', you have to
# force-set DART_TESTING_TIMEOUT in the cache on the first configure.  But to
# avoid rescaling DART_TESTING_TIMEOUT up and up on reconfigures, you have to
# force-set DART_TESTING_TIMEOUT back to the user's input value.  The only
# disadvantage of this approach (other than it is a hack to get around a CMake
# bug) is that you loose the user's documentation string, in case they set
# that with a set( ... CACHE ...) statement in an input *.cmake file.


# Determines if a package should be processed
#
function(tribits_determine_if_process_package  PACKAGE_NAME
  PROCESS_PACKAGE_OUT  PACKAGE_ENABLE_STR_OUT
  )

  set(PROCESS_PACKAGE FALSE)
  set(PACKAGE_ENABLE_STR "")

  if (${PACKAGE_NAME}_SUBPACKAGES)
    # Process the package if any of the subpackages are enable
    foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
      set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
      if (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        set(PROCESS_PACKAGE TRUE)
        append_string_var_with_sep(PACKAGE_ENABLE_STR ", " ${TRIBITS_SUBPACKAGE})
      endif()
    endforeach()
  else()
    # If the package itself is enabled, of course process it
    if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
      set(PROCESS_PACKAGE TRUE)
      append_string_var_with_sep(PACKAGE_ENABLE_STR ", " "Libs")
    endif()
  endif()

  # If subpackages or package is enabled, then check tests/examples
  if (PROCESS_PACKAGE)
    if (${PACKAGE_NAME}_ENABLE_TESTS)
      append_string_var_with_sep(PACKAGE_ENABLE_STR ", " "Tests")
    endif()
    if (${PACKAGE_NAME}_ENABLE_EXAMPLES)
      append_string_var_with_sep(PACKAGE_ENABLE_STR ", " "Examples")
    endif()
  endif()

  set(${PROCESS_PACKAGE_OUT} ${PROCESS_PACKAGE} PARENT_SCOPE)
  set(${PACKAGE_ENABLE_STR_OUT} ${PACKAGE_ENABLE_STR} PARENT_SCOPE)

endfunction()


# Reads in the project's version file into the current scope
#
macro(tribits_project_read_version_file  PROJECT_SOURCE_DIR_IN)
  set(PROJECT_VERSION_FILE ${PROJECT_SOURCE_DIR_IN}/Version.cmake)
  if (EXISTS ${PROJECT_VERSION_FILE})
    # Set REPOSITORY_NAME in case Version.cmake is written generically!
    set(REPOSITORY_NAME ${PROJECT_NAME})
    tribits_trace_file_processing(PROJECT  INCLUDE  "${PROJECT_VERSION_FILE}")
    include(${PROJECT_VERSION_FILE})
  endif()
endmacro()


# Read in and the Repository's specific Version.cmake file and then configure
# its ${REPO_NAME}_version.h file.
#
# The file ${REPO_NAME}_version.h is only configured if the repository contains
# the files Version.cmake and Copyright.txt
#
# NOTE: This is done as a function so that the read-in version variables don't
# bleed into the outer scope.
#
function(tribits_repository_configure_version_header_file
  REPOSITORY_NAME  REPOSITORY_DIR  ADD_INSTALL_TARGET
  OUTPUT_VERSION_HEADER_FILE
  )

  if (TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE_DEBUG_DUMP)
    message("TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE: "
      "'${REPOSITORY_NAME}'  '${REPOSITORY_DIR}"
      "  '${OUTPUT_VERSION_HEADER_FILE}'")
  endif()

  string(TOUPPER ${REPOSITORY_NAME} REPOSITORY_NAME_UC)

  tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR} ${REPOSITORY_DIR}
    REPOSITORY_ABS_DIR)

  set(REPOSITORY_VERSION_FILE ${REPOSITORY_ABS_DIR}/Version.cmake)
  set(REPOSITORY_COPYRIGHT_FILE ${REPOSITORY_ABS_DIR}/Copyright.txt)

  if (EXISTS ${REPOSITORY_VERSION_FILE} AND EXISTS ${REPOSITORY_COPYRIGHT_FILE})

    # Read the copyright header info
    tribits_trace_file_processing(REPOSITORY  READ  "${REPOSITORY_COPYRIGHT_FILE}")
    file(READ "${REPOSITORY_COPYRIGHT_FILE}" REPOSITORY_COPYRIGHT_HEADER)

    # Read the version variables and translate into standard form
    tribits_trace_file_processing(REPOSITORY  INCLUDE  "${REPOSITORY_VERSION_FILE}")
    include(${REPOSITORY_VERSION_FILE})
    set(REPOSITORY_MAJOR_VERSION ${${REPOSITORY_NAME}_MAJOR_VERSION})
    set(REPOSITORY_MAJOR_MINOR_VERSION ${${REPOSITORY_NAME}_MAJOR_MINOR_VERSION})
    set(REPOSITORY_VERSION_STRING ${${REPOSITORY_NAME}_VERSION_STRING})

    # Configure the file with everything set
    if (TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE_DEBUG_DUMP)
      message("-- Writing the file ${OUTPUT_VERSION_HEADER_FILE} ...")
    endif()
    tribits_trace_file_processing(REPOSITORY  CONFIGURE  "${OUTPUT_VERSION_HEADER_FILE}")
    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_PACKAGE_ARCH_DIR}/Tribits_version.h.in
      ${OUTPUT_VERSION_HEADER_FILE})

    if (ADD_INSTALL_TARGET)
      # Install version header file
      tribits_install_headers(HEADERS  ${OUTPUT_VERSION_HEADER_FILE})
    endif()

  endif()

endfunction()


# Configure each of the Repositories' version header file
#
function(tribits_repository_configure_all_version_header_files)
  #print_var(ARGN)
  foreach(REPO ${ARGN})
    tribits_get_repo_name_dir(${REPO}  REPO_NAME  REPO_DIR)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Considering configuring version file for '${REPO_NAME}'")
    endif()
    tribits_repository_configure_version_header_file( ${REPO_NAME}  ${REPO_DIR}  TRUE
      "${${PROJECT_NAME}_BINARY_DIR}/${REPO_DIR}/${REPO_NAME}_version.h")
  endforeach()
endfunction()


# Generate the VersionDate.cmake and ${REPO_NAME}_version_date.h files for a
# TriBITS Repository
#
# NOTE: This is done as a function so that the read-in version variables don't
# bleed into the outer scope.
#
function(tribits_repository_configure_version_date_files
  REPOSITORY_NAME  REPOSITORY_DIR  ADD_INSTALL_TARGET
  )

  if (TRIBITS_REPOSITORY_CONFIGURE_VERSION_DATE_FILES_DEBUG_DUMP)
    message("TRIBITS_REPOSITORY_CONFIGURE_VERSION_DATE_FILES: "
      "'${REPOSITORY_NAME}'  '${REPOSITORY_DIR}" )
  endif()

  string(TOUPPER ${REPOSITORY_NAME} REPOSITORY_NAME_UC)

  tribits_set_base_repo_dir(${PROJECT_SOURCE_DIR} ${REPOSITORY_DIR}
    REPO_SOURCE_ABS_DIR)

  tribits_set_base_repo_dir(${PROJECT_BINARY_DIR} ${REPOSITORY_DIR}
    REPO_BINARY_ABS_DIR)

  set(REPO_GIT_VERSION_DATE)
  if (NOT IS_DIRECTORY "${REPO_SOURCE_ABS_DIR}/.git")
    message("-- NOTE: Can't fill in version date files for ${REPOSITORY_NAME} since"
      " ${REPO_SOURCE_ABS_DIR}/.git/ does not exist!")
  elseif (GIT_VERSION_STRING VERSION_LESS "2.10.0")
    message("-- NOTE: Can't fill in version date files for ${REPOSITORY_NAME} since"
      " GIT_VERSION_STRING=${GIT_VERSION_STRING} < 2.10.0")
  else()
    # Generate the version date integer
    tribits_get_raw_git_commit_utc_time("${REPO_SOURCE_ABS_DIR}" "HEAD"
      REPO_GIT_COMMIT_UTC_TIME)
    tribits_get_version_date_from_raw_git_commit_utc_time("${REPO_GIT_COMMIT_UTC_TIME}"
      REPO_GIT_VERSION_DATE )
  endif()

  if (REPO_GIT_VERSION_DATE)
    # Configure the VersionDate.cmake file in the repo binary dir and include it
    set(REPO_VERSION_DATE_CMAKE_FILE "${REPO_BINARY_ABS_DIR}/VersionDate.cmake")
    configure_file(
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_PACKAGE_ARCH_DIR}/VersionDate.cmake.in"
      "${REPO_VERSION_DATE_CMAKE_FILE}" )
    tribits_trace_file_processing(REPOSITORY  INCLUDE  "${REPO_VERSION_DATE_CMAKE_FILE}")
    include(${REPO_VERSION_DATE_CMAKE_FILE})
  endif()

  # Configure the <RepoName>_version_date.h file in the repo binary dir
  if (REPO_GIT_VERSION_DATE)
    set(REPOSITORY_VERSION_DATE_MACRO_DEF
      "#define ${REPOSITORY_NAME_UC}_VERSION_DATE ${REPO_GIT_VERSION_DATE}" )
  else()
    set(REPOSITORY_VERSION_DATE_MACRO_DEF
      "#undef ${REPOSITORY_NAME_UC}_VERSION_DATE" )
  endif()
  set(REPO_VERSION_DATE_HEADER_FILE "${REPO_BINARY_ABS_DIR}/${REPOSITORY_NAME}_version_date.h")
  tribits_trace_file_processing(REPOSITORY  CONFIGURE  "${REPO_VERSION_DATE_HEADER_FILE}")
  configure_file(
    "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_PACKAGE_ARCH_DIR}/Tribits_version_date.h.in"
    "${REPO_VERSION_DATE_HEADER_FILE}" )

  # Add the install target for <RepoName>_version_date.h
  if (ADD_INSTALL_TARGET)
    tribits_install_headers(HEADERS  ${REPO_VERSION_DATE_HEADER_FILE})
  endif()

endfunction()


# Configure each of the Repositories' version date files
#
function(tribits_repository_configure_all_version_date_files)
  #print_var(ARGN)
  if (${PROJECT_NAME}_GENERATE_VERSION_DATE_FILES)
    foreach(REPO ${ARGN})
      tribits_get_repo_name_dir(${REPO}  REPO_NAME  REPO_DIR)
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("Considering configuring version date files for '${REPO_NAME}'")
      endif()
      tribits_repository_configure_version_date_files(${REPO_NAME} ${REPO_DIR} TRUE)
    endforeach()
  endif()
endfunction()


# Configure the enabled packages
#
# This macro actually calls add_subdirectory(<packageDir>) on the enabled
# TriBITS packages.
#
macro(tribits_configure_enabled_packages)

  tribits_config_code_start_timer(CONFIGURE_PACKAGES_TIME_START_SECONDS)

  #
  # A) Global variable initialization
  #

  global_null_set(${PROJECT_NAME}_LIBRARIES "")
  global_null_set(${PROJECT_NAME}_ETI_PACKAGES "")

  #
  # B) Define the source and binary directories for all of the packages that
  # have been enabled.  These are used to allow packages to refer to each
  # other even downstream packages (which is pretty messed up really).
  #

  tribits_filter_package_list_from_var(${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
    INTERNAL  ON  NONEMPTY  ${PROJECT_NAME}_enabledInternalTopLevelPackages)

  foreach(TRIBITS_PACKAGE  IN LISTS  ${PROJECT_NAME}_enabledInternalTopLevelPackages)

    # Get all the package sources independent of whether they are enabled or not.
    # There are some messed up packages that grab parts out of unrelated
    # downstream packages that might not even be enabled.  To support this,
    # allow this.

    tribits_determine_if_process_package(${TRIBITS_PACKAGE}
       PROCESS_PACKAGE  PACKAGE_ENABLE_STR)

    if (PROCESS_PACKAGE)

      if (${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR)
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          print_var(${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR)
        endif()
        if(IS_ABSOLUTE ${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
          set(${TRIBITS_PACKAGE}_BINARY_DIR ${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
        else()
          set(${TRIBITS_PACKAGE}_BINARY_DIR
            ${CMAKE_CURRENT_BINARY_DIR}/${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
        endif()
      else()
        set(${TRIBITS_PACKAGE}_BINARY_DIR
          ${CMAKE_CURRENT_BINARY_DIR}/${${TRIBITS_PACKAGE}_REL_SOURCE_DIR})
      endif()
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        print_var(${TRIBITS_PACKAGE}_BINARY_DIR)
      endif()

    endif()

  endforeach()

  #
  # C) Loop over all of the packages and process their CMakeLists.txt files if
  # they are enabled or if any of their subpackages are enabled.
  #

  # Include these here so they don't need to be included in each package's
  # CMakeLists.txt files.
  include(TribitsPackageMacros)
  include(TribitsSubPackageMacros)
  include(AddSubdirectories)

  set(CONFIGURED_A_PACKAGE FALSE)
  set(ENABLED_PACKAGE_LIBS_TARGETS)

  # Tell packages that are also repos they are being processed as a package.
  set(TRIBITS_PROCESSING_PACKAGE TRUE)

  foreach(TRIBITS_PACKAGE  IN LISTS  ${PROJECT_NAME}_enabledInternalTopLevelPackages)

    tribits_determine_if_process_package(${TRIBITS_PACKAGE}
      PROCESS_PACKAGE  PACKAGE_ENABLE_STR)

    if (PROCESS_PACKAGE)

      message("Processing enabled top-level package: ${TRIBITS_PACKAGE} (${PACKAGE_ENABLE_STR})")

      if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

        tribits_package_config_code_start_timer(PROCESS_THIS_PACKAGE_TIME_START_SECONDS)

        set(PACKAGE_NAME ${TRIBITS_PACKAGE}) # Used in CMake code in downstream package
        set(PARENT_PACKAGE_NAME ${TRIBITS_PACKAGE})
        string(TOUPPER "${PARENT_PACKAGE_NAME}" PARENT_PACKAGE_NAME_UC)

        if (NOT EXISTS ${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt)
          message(FATAL_ERROR
            "Error, the file ${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt does not exist!")
        endif()

        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          print_var(${TRIBITS_PACKAGE}_SOURCE_DIR)
          print_var(${TRIBITS_PACKAGE}_BINARY_DIR)
        endif()

        set(TRIBITS_PACKAGE_CMAKELIST_FILE
          "${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt")
        tribits_trace_file_processing(PACKAGE  ADD_SUBDIR
          "${TRIBITS_PACKAGE_CMAKELIST_FILE}")
        if (NOT ${TRIBITS_PACKAGE}_SOURCE_DIR STREQUAL ${PROJECT_NAME}_SOURCE_DIR)
          add_subdirectory(${${TRIBITS_PACKAGE}_SOURCE_DIR}
            ${${TRIBITS_PACKAGE}_BINARY_DIR})
        else()
          include("${TRIBITS_PACKAGE_CMAKELIST_FILE}")
        endif()
        if ((NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS) AND
	    (NOT TARGET ${PACKAGE_NAME}::all_libs)
          )
          tribits_report_invalid_tribits_usage(
            "ERROR: Forgot to call tribits_package_postprocess() in"
            " ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
        endif()

        list(APPEND ENABLED_PACKAGE_LIBS_TARGETS ${TRIBITS_PACKAGE}::all_libs)
        list(APPEND ${PROJECT_NAME}_LIBRARIES ${${TRIBITS_PACKAGE}_LIBRARIES})

        tribits_package_config_code_stop_timer(PROCESS_THIS_PACKAGE_TIME_START_SECONDS
          "-- Total time to configure top-level package ${TRIBITS_PACKAGE}")

      endif()

      set(CONFIGURED_A_PACKAGE TRUE)

    endif()

  endforeach()

  #
  # D) Loop backwards over ETI packages if ETI is enabled
  #

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    # Do this regardless of whether project level ETI is enabled
    if("${${PROJECT_NAME}_ETI_PACKAGES}" STREQUAL "")
      message("\nNo ETI support requested by packages.\n")
    else()
      #if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("\nProcessing explicit instantiation support for enabled packages ...\n")
      #endif()
      set(REVERSE_ETI_LIST ${${PROJECT_NAME}_ETI_PACKAGES})
      list(REVERSE REVERSE_ETI_LIST)
      foreach(PACKAGE_NAME ${REVERSE_ETI_LIST})
        message("Processing ETI support: ${PACKAGE_NAME}")
        tribits_package_config_code_start_timer(PROCESS_ETI_START_SECONDS)
        set(ETIFILE
          ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/ExplicitInstantiationSupport.cmake)
        if(NOT EXISTS "${ETIFILE}")
          message(FATAL_ERROR
            "Could not find ${PACKAGE_NAME} ETI support file ${ETIFILE}")
        endif()
        tribits_trace_file_processing(PACKAGE  INCLUDE  "${ETIFILE}")
        include("${ETIFILE}")
        tribits_package_config_code_stop_timer(PROCESS_ETI_START_SECONDS
          "-- Time to process ETI support for package ${PACKAGE_NAME}")
      endforeach()
    endif()

  endif()

  #
  # E) Check if no packages are enabled and if that is allowed
  #

  advanced_set( ${PROJECT_NAME}_ALLOW_NO_PACKAGES ON
    CACHE BOOL "Allow configuration to finish even if no packages are enabled")

  if (NOT CONFIGURED_A_PACKAGE)
    if (${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      set(MSG_TYPE WARNING)
    else()
      set(MSG_TYPE ERROR)
    endif()
    message(
      "\n***"
      "\n*** ${MSG_TYPE}:  There were no packages configured so no libraries"
        " or tests/examples will be built!"
      "\n***\n"
      )
    if (NOT ${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      message(SEND_ERROR "Stopping configure!")
    endif()
  else()
    assert_and_touch_defined(${PROJECT_NAME}_ALLOW_NO_PACKAGES)
  endif()

  #
  # F) Process the global variables and other cleanup
  #

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    remove_global_duplicates(${PROJECT_NAME}_LIBRARIES)

    # Add global 'libs' target
    if(ENABLED_PACKAGE_LIBS_TARGETS)
      list(REVERSE ENABLED_PACKAGE_LIBS_TARGETS)
      # Make it so when no packages are enabled it is not a cmake error
      if (NOT TARGET ${PROJECT_NAME}_libs)
        add_custom_target(${PROJECT_NAME}_libs)
        add_dependencies(${PROJECT_NAME}_libs ${ENABLED_PACKAGE_LIBS_TARGETS})
      endif()
      add_custom_target(libs)
      add_dependencies(libs ${ENABLED_PACKAGE_LIBS_TARGETS})
    endif()

    # Add empty <PackageName>_libs targets for top-level packages if asked
    if (${PROJECT_NAME}_DEFINE_MISSING_PACKAGE_LIBS_TARGETS)
      foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
        if (NOT TARGET ${TRIBITS_PACKAGE}_libs)
          add_custom_target(${TRIBITS_PACKAGE}_libs
            COMMENT "Dummy target for ${TRIBITS_PACKAGE}_libs that builds nothing!")
        endif()
      endforeach()
    endif()
    # NOTE: For motivation for above, see the comment about the setting of
    # ${PROJECT_NAME}_DEFINE_MISSING_PACKAGE_LIBS_TARGETS=ON in
    # package-by-package mode in tribits_ctest_driver().  This option is
    # purposefully not documented and not defined as a cache variable since it
    # is an internal TriBITS implementation detail.

  endif()

  tribits_config_code_stop_timer(CONFIGURE_PACKAGES_TIME_START_SECONDS
    "\nTotal time to configure enabled packages")

endmacro()


# Create custom 'install_package_by_package' target
#
function(tribits_add_install_package_by_package_target)

  set(TRIBITS_ENABLED_PACKAGES_BINARY_DIRS)
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
    list(APPEND TRIBITS_ENABLED_PACKAGES_BINARY_DIRS "${${TRIBITS_PACKAGE}_BINARY_DIR}")
  endforeach()

  set(tribits_install_src
    "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")

  configure_file(
    "${tribits_install_src}/cmake_pbp_install.cmake.in"
    cmake_pbp_install.cmake
    @ONLY )

  advanced_set(${PROJECT_NAME}_INSTALL_PBP_RUNNER "" CACHE FILEPATH
    "Program used to run cmake -P cmake_pbp_install.cmake to change user for 'install_package_by_package' target")

  add_custom_target(install_package_by_package
   ${${PROJECT_NAME}_INSTALL_PBP_RUNNER}
    ${CMAKE_COMMAND} -P cmake_pbp_install.cmake
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

endfunction()


macro(tribits_setup_for_installation)

  # Set up to install <Package>Config.cmake, <Project>Config.cmake, and export
  # makefiles.
  add_subdirectory(
    "${${PROJECT_NAME}_TRIBITS_DIR}/core/installation/add_project_install_commands"
    add_project_install_commands)

  # Set up for fixing group and permissions after the install
  add_subdirectory(
    "${${PROJECT_NAME}_TRIBITS_DIR}/core/installation/add_install_group_and_perms_fixups"
    add_install_group_and_perms_fixups)

  # Create custom 'install_package_by_package' target
  tribits_add_install_package_by_package_target()

endmacro()

#  LocalWords:
#  LocalWords: Sandia SANDIA Redistributions
#  LocalWords: tribits TriBITS TRIBITS
#  LocalWords: cmake CMake CMAKE CMakeCache CMakeFiles
#  LocalWords: ctest CPACK
#  LocalWords: foreach endforeach endif endmacro
#  LocalWords: BOOL
#  LocalWords: libs LIBS config PackageName SUBPACKAGES nonenabled
