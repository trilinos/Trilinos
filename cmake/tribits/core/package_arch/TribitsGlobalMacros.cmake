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

# Standard TriBITS system includes
INCLUDE(TribitsConstants)
INCLUDE(TribitsProcessExtraRepositoriesList)
INCLUDE(TribitsProcessPackagesAndDirsLists)
INCLUDE(TribitsProcessTplsLists)
INCLUDE(TribitsAdjustPackageEnables)
INCLUDE(TribitsSetupMPI)
INCLUDE(TribitsTestCategories)
INCLUDE(TribitsGeneralMacros)
INCLUDE(TribitsAddTestHelpers)
INCLUDE(TribitsVerbosePrintVar)
INCLUDE(TribitsProcessEnabledTpl)
INCLUDE(TribitsInstallHeaders)

# Standard TriBITS utilities includes
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(AdvancedOption)
INCLUDE(AdvancedSet)
INCLUDE(AppendStringVar)
INCLUDE(AppendStringVarWithSep)
INCLUDE(AssertAndTouchDefined)
INCLUDE(CMakeBuildTypesList)
INCLUDE(FindListElement)
INCLUDE(GlobalNullSet)
INCLUDE(PrintNonemptyVar)
INCLUDE(PrintVar)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(Split)
INCLUDE(TimingUtils)
INCLUDE(SetDefaultAndFromEnv) # Used by some call-back files

# Standard CMake includes
INCLUDE(CheckIncludeFileCXX)

# Include here so it does not need to be included in each individual
# FindTPL<TPLNAME>.cmake file over and over.
INCLUDE(TribitsTplFindIncludeDirsAndLibraries)
INCLUDE(TribitsTplDeclareLibraries) # Deprecated
# ABOVE: We need to include TribitsTplDeclareLibraries.cmake until all client
# projects stop using it.


#
# Assert and setup project binary directory and other project variables.
#
MACRO(TRIBITS_ASSERT_AND_SETUP_PROJECT_AND_STATIC_SYSTEM_VARS)

  APPEND_STRING_VAR(IN_SOURCE_ERROR_COMMON_MSG
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

  IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt")
    MESSAGE(FATAL_ERROR "ERROR! "
      "The file ${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt exists from a"
      " likely prior attempt to do an in-source build."
      "${IN_SOURCE_ERROR_COMMON_MSG}"
      )
  ENDIF()

  IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    MESSAGE(FATAL_ERROR "ERROR! "
      "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
      " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
      "\n${PROJECT_NAME} does not support in source builds!\n"
      "NOTE: You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under"
      " the source directory for ${PROJECT_NAME} or you will not be able to configure ${PROJECT_NAME} correctly!"
      "${IN_SOURCE_ERROR_COMMON_MSG}"
      )
  ENDIF()

  STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)
  SET(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
  SET(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
  PRINT_VAR(PROJECT_SOURCE_DIR)
  PRINT_VAR(PROJECT_BINARY_DIR)
  PRINT_VAR(${PROJECT_NAME}_TRIBITS_DIR)
  PRINT_VAR(TriBITS_VERSION_STRING)
  # Above, we put these in the cache so we can grep them out of the cache file

  #
  # Print some basic static info provided by CMake automatically
  #

  PRINT_VAR(CMAKE_VERSION)
  PRINT_VAR(CMAKE_GENERATOR)

ENDMACRO()


#
# Set up some really basic system variables.
#
# This macro needs to be called *before* the user *.cmake option files are
# read in so that there is an opportunity to override these.
#
MACRO(TRIBITS_SETUP_BASIC_SYSTEM_VARS)

  # CMAKE_HOST_SYSTEM_NAME is provided by CMake automatically but can actually
  # be overridden in the cache.
  PRINT_VAR(CMAKE_HOST_SYSTEM_NAME)

  SITE_NAME(${PROJECT_NAME}_HOSTNAME)
  MARK_AS_ADVANCED(${PROJECT_NAME}_HOSTNAME)
  PRINT_VAR(${PROJECT_NAME}_HOSTNAME)

  # NOTE: CMAKE_HOST_SYSTEM_NAME and ${PROJECT_NAME}_HOSTNAME are used by
  # TRIBITS_ADD[_ADVANCED]_TEST() to include/exclude tests based in the
  # arguments HOSTS, XHOSTS, HOSTTYPES, AND XHOSTTYPES.

ENDMACRO()


#
# Define an option to include a file that reads in a bunch of options
#
#

MACRO(TRIBITS_READ_IN_OPTIONS_FROM_FILE)

  SET( ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE "" CACHE FILEPATH
    "Name of an optional file that is included first to define any cmake options with SET( ... CACHE ...) calls.  NOTE: paths can be separated by commas instead of semicolons but paths cannot contain commas."
    )

  SPLIT("${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE}"  ","
    ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE)

  SET( ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_ALL
    ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE}
    ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_APPEND})

  FOREACH (CONFIG_OPTS_FILE ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE_ALL})
    MESSAGE("-- " "Reading in configuration options from ${CONFIG_OPTS_FILE} ...")
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${CONFIG_OPTS_FILE}")
    INCLUDE(${CONFIG_OPTS_FILE})
  ENDFOREACH()


ENDMACRO()


#
# Define all of the standard global package architecture options.
#

MACRO(TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS)

  SET( ${PROJECT_NAME}_ENABLE_ALL_PACKAGES OFF CACHE BOOL
    "Enable all packages PT packages (ST packages as well if ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE is true)." )

  SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON CACHE BOOL
    "Recursively enable all optional packages for set of enabled packages." )

  SET( ${PROJECT_NAME}_INSTALL_EXECUTABLES ON CACHE BOOL
    "Enable the installation of executables provided by the ${PROJECT_NAME} packages." )

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL
    "Recursively enable all packages that have required or optional dependencies for set of enabled packages." )

  IF (${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT STREQUAL "")
    SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT OFF)
  ENDIF()
  SET(${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES
    ${${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES_DEFAULT}
    CACHE BOOL
    "Disable (and printing warning) for enabled packages that have hard-disabled upstream dependencies.  Otherwise, is to raises a fatal configure failure." )

  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_TESTS ""
    "Enable tests in all packages  (set to ON, OFF, or leave empty)." )

  SET_CACHE_ON_OFF_EMPTY(${PROJECT_NAME}_ENABLE_EXAMPLES ""
    "Enable examples in all packages  (set to ON, OFF, or leave empty).  If left empty, then this will be set to ON if ${PROJECT_NAME}_ENABLE_TESTS=ON" )

  IF (${PROJECT_NAME}_ENABLE_TESTS AND ${PROJECT_NAME}_ENABLE_EXAMPLES STREQUAL "")
    MESSAGE(STATUS "Setting ${PROJECT_NAME}_ENABLE_EXAMPLES=ON because ${PROJECT_NAME}_ENABLE_TESTS=ON")
    SET(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
    "Set to empty all package enables (set to OFF at end)." )

  ADVANCED_OPTION(${PROJECT_NAME}_REMOVE_DEFAULT_PACKAGE_DISABLES
    "Removes all default disables from the packages list.  Used for testing etc."
    OFF )

  IF ("${${PROJECT_NAME}_ENABLE_C_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_C_DEFAULT ON)
  ENDIF()
  ADVANCED_OPTION(${PROJECT_NAME}_ENABLE_C
    "Enable the C compiler and related code"
    ${${PROJECT_NAME}_ENABLE_C_DEFAULT} )

  IF ("${${PROJECT_NAME}_C_Standard_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_C_Standard_DEFAULT c99)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_C_Standard
    ${${PROJECT_NAME}_C_Standard_DEFAULT}
    CACHE STRING
    "The standard <cstd> to use in --std=<cstd> for GCC compilers." )

  IF ("${${PROJECT_NAME}_ENABLE_CXX_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_CXX_DEFAULT ON)
  ENDIF()
  ADVANCED_OPTION(${PROJECT_NAME}_ENABLE_CXX
    "Enable the C++ compiler and related code"
    ${${PROJECT_NAME}_ENABLE_CXX_DEFAULT} )

  IF ("${${PROJECT_NAME}_ENABLE_CXX11_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_CXX11_DEFAULT OFF)
  ENDIF()
  ADVANCED_OPTION(${PROJECT_NAME}_ENABLE_CXX11
    "Enable the C++11 compiler options and related code (see ${PROJECT_NAME}_CXX11_FLAGS)"
    ${${PROJECT_NAME}_ENABLE_CXX11_DEFAULT} )

  IF ("${${PROJECT_NAME}_ENABLE_Fortran_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
  ENDIF()

  OPTION(${PROJECT_NAME}_ENABLE_Fortran
    "Enable the Fortran compiler and related code"
    ${${PROJECT_NAME}_ENABLE_Fortran_DEFAULT} )

  ADVANCED_OPTION(${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST
    "Skip the Fortran/C++ compatibility test"
    OFF )

  IF ("${${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT TRUE)
  ELSE()
    SET(${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT FALSE)
  ENDIF()
  ADVANCED_SET(
    ${PROJECT_NAME}_SET_INSTALL_RPATH ${${PROJECT_NAME}_SET_INSTALL_RPATH_DEFAULT}
    CACHE BOOL
    "If TRUE, then set RPATH on installed binaries will set to ${PROJECT_NAME}_INSTALL_LIB_DIR automatically"
    )

  IF ("${CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT}" STREQUAL "")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT TRUE)
  ELSE()
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT FALSE)
  ENDIF()
  ADVANCED_SET(
    CMAKE_INSTALL_RPATH_USE_LINK_PATH ${CMAKE_INSTALL_RPATH_USE_LINK_PATH_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then the RPATH for external shared libs will be embedded in installed libs and execs."
    )

  ADVANCED_SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ""
    CACHE STRING
    "Extra flags added to the end of every linked executable"
    )

  # OpenMP is similar to a TPL in some respects, but requires only compiler
  # flags to enable

  OPTION(${PROJECT_NAME}_ENABLE_OpenMP
    "Build with OpenMP support." OFF)

  IF (CMAKE_GENERATOR STREQUAL "Ninja")
    IF("${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT}" STREQUAL "")
      SET(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT ON)
    ENDIF()
    SET(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES
      ${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES_DEFAULT} CACHE BOOL
      "Generate dummy makefiles to call ninja in every bulid subdirectory (requires CMake 3.7.0 or newer)." )
  ENDIF()
  IF ("${${PROJECT_NAME}_WRITE_NINJA_MAKEFILES}" STREQUAL "")
    SET(${PROJECT_NAME}_WRITE_NINJA_MAKEFILES OFF)
  ENDIF()
  
  IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    SET(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT ON)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT OFF)
  ENDIF()
  SET(${PROJECT_NAME}_ENABLE_DEBUG ${${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT} CACHE BOOL
    "Enable debug checking for ${PROJECT_NAME} packages.  Off by default unless CMAKE_BUILD_TYPE=\"DEBUG\"." )

  IF (${PROJECT_NAME}_ENABLE_DEBUG)
    SET(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT ON)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT OFF)
  ENDIF()
  SET(${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG
    ${${PROJECT_NAME}_ENABLE_CONFIGURE_DEBUG_DEFAULT} CACHE BOOL
    "Enable debug checking of the process which finds errors in the project's CMake files (off by default unless ${PROJECT_NAME}_ENABLE_DEBUG=ON)." )

  IF ("${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT "WARNING")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS
    ${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS_DEFAULT}
    CACHE STRING
    "Determins how unparsed arguments for TriBITS functions that use CMAKE_PARASE_ARUMENTS() internally are handled.  Valid choices are 'WARNING', 'SEND_ERROR', and 'FATAL_ERROR'.  The default is 'SEND_ERROR'."
    )
  IF (
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "WARNING")
     AND
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "SEND_ERROR")
     AND
    (NOT ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS STREQUAL "FATAL_ERROR")
    )
    MESSAGE(FATAL_ERROR "Error, the value of"
      " ${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS ="
      " '${${PROJECT_NAME}_CHECK_FOR_UNPARSED_ARGUMENTS}' is invalid!"
      " Valid valules include 'WANRING', 'SEND_ERROR', and 'FATAL_ERROR'"
      )
  ENDIF()

  SET(${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR ON
    CACHE BOOL
    "Enable support for Teuchos Time Monitors in all Trilinos packages that support it."
    )

  ADVANCED_SET(${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS ON
    CACHE BOOL
    "Show warnings about deprecated code"
    )

  ADVANCED_SET(${PROJECT_NAME}_HIDE_DEPRECATED_CODE OFF
    CACHE BOOL
    "Show warnings about deprecated code"
    )

  ADVANCED_SET(${PROJECT_NAME}_VERBOSE_CONFIGURE OFF
    CACHE BOOL
    "Make the ${PROJECT_NAME} configure process verbose."
    )

  IF ("${${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT  ${${PROJECT_NAME}_VERBOSE_CONFIGURE})
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_TRACE_ADD_TEST ${${PROJECT_NAME}_TRACE_ADD_TEST_DEFAULT}
    CACHE BOOL
    "Show a configure time trace of every test added or not added any why (one line)." )

  ADVANCED_OPTION(${PROJECT_NAME}_DUMP_LINK_LIBS
    "Dump the link libraries for every library and executable created."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  ADVANCED_SET(${PROJECT_NAME}_TRACE_FILE_PROCESSING
    ${${PROJECT_NAME}_VERBOSE_CONFIGURE}
    CACHE BOOL
    "Print out when all of the various files get processed."
    )

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION OFF
    CACHE BOOL
    "Enable explicit template instantiation in all packages that support it"
    )

  IF (USE_XSDK_DEFAULTS)
    # Need to set BUILD_SHARED_LIBS default here based on USE_XSDK_DEFAULTS
    # and not in TRIBITS_SETUP_ENV() in case there is logic in TriBITS or
    # project-specific files that depends on this var getting set here!
    SET(BUILD_SHARED_LIBS_DEFAULT  TRUE)
    IF ("${BUILD_SHARED_LIBS}" STREQUAL "")
      MESSAGE("-- " "XSDK: Setting default BUILD_SHARED_LIBS=TRUE")
    ENDIF()
  ELSE()
    SET(BUILD_SHARED_LIBS_DEFAULT  FALSE)
  ENDIF()
  SET(BUILD_SHARED_LIBS  ${BUILD_SHARED_LIBS_DEFAULT}
    CACHE  BOOL   "Build shared libraries or not.")

  IF ("${${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT  FALSE)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS
    ${${PROJECT_NAME}_TPL_SYSTEM_INCLUDE_DIRS_DEFAULT}
    CACHE BOOL
    "If set TRUE, then 'SYSTEM' will be passed into INCLUDE_DIRECTORIES() for TPL includes.")

  ADVANCED_SET(TPL_FIND_SHARED_LIBS ON CACHE BOOL
    "If ON, then the TPL system will find shared libs if they exist, otherwise will only find static libs." )

  ADVANCED_SET(${PROJECT_NAME}_LINK_SEARCH_START_STATIC OFF CACHE BOOL
    "If ON, then the property LINK_SEARCH_START_STATIC will be added to all executables." )

  ADVANCED_SET(${PROJECT_NAME}_LIBRARY_NAME_PREFIX ""
    CACHE STRING
    "Prefix for all ${PROJECT_NAME} library names. If set to, for example, 'prefix_',
    libraries will be named and installed as 'prefix_<libname>.*'.  Default is '' (no prefix)."
    )

  IF ("${${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT FALSE)
  ENDIF()
  ADVANCED_SET( ${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS
    ${${PROJECT_NAME}_MUST_FIND_ALL_TPL_LIBS_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then all of the TPL libs must be found for every enabled TPL."
    )

  IF ("${${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT FALSE)  # Maintain backward compatibility
  ENDIF()
  ADVANCED_SET( ${PROJECT_NAME}_USE_GNUINSTALLDIRS
    ${${PROJECT_NAME}_USE_GNUINSTALLDIRS_DEFAULT}
    CACHE BOOL
    "If set to TRUE, then CMake GNUInstallDris modules is used to pick standard install paths by default."
    )

  IF ("${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}" STREQUAL "")
    # Assume the TriBITS project wants to install headers and libraries by default
    SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT ON)
  ENDIF()

  ADVANCED_SET(${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS
    ${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}
    CACHE BOOL
    "Install libraries and headers (default is ${${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS_DEFAULT}).  NOTE: Shared libraries are always installed since they are needed by executables."
    )

  IF ("${${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT}" STREQUAL "")
    IF(WIN32 AND NOT CYGWIN)
      SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT OFF)
    ELSE()
      SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT ON)
    ENDIF()
  ENDIF()

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES
    ${${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT}
    CACHE BOOL
    "Determines if export makefiles will be created and installed."
    )

  # Creating <Package>Config.cmake files is currently *very* expensive for large
  # TriBITS projects so we disable this by default for TriBITS.
  IF ("${${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT OFF)
  ENDIF()

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
    ${${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT}
    CACHE BOOL
    "Determines if ${PROJECT_NAME}Config.cmake and <PACKAGE>Config.cmake files are created or not."
    )

  IF (NOT ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT)
    # We need to generate the dependency logic for export dependency files if
    # asked.
    IF (${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES OR
      ${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
      )
      SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)
    ELSE()
      SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT OFF)
    ENDIF()
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES
     ${${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT} CACHE BOOL
    "Generate packages dependency data-structures needed for dependency export files." )

  # ${PROJECT_NAME}_ELEVATE_SS_TO_PS is depreciated!
  IF (${PROJECT_NAME}_ELEVATE_SS_TO_PS_DEFAULT)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("-- " "WARNING: ${PROJECT_NAME}_ELEVATE_SS_TO_PS_DEFAULT is deprecated."
        "  Use ${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT instead!")
    ENDIF()
    SET(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT ON)
  ENDIF()

  IF ("${${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET( ${PROJECT_NAME}_ELEVATE_ST_TO_PT
    ${${PROJECT_NAME}_ELEVATE_ST_TO_PT_DEFAULT}
    CACHE BOOL
    "Elevate all defined ST SE packages to PT packages." )

  IF ("${${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_CPACK_PACKAGING
     ${${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT}
     CACHE BOOL
    "Enable support for creating a distribution using CPack" )

  IF ("${${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT TRUE)
  ENDIF()
  ADVANCED_SET( ${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION
    ${${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT}
    CACHE BOOL
    "Excluded disabled packages from the CPack-generated distribution.")

  IF ("${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET(
    ${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE
    ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE_DEFAULT}
    CACHE BOOL
    "Allow Secondary Tested (ST) packages and code to be implicitly enabled." )

  IF ("${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT NIGHTLY)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_TEST_CATEGORIES
     ${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}
     CACHE STRING
    "List of categories of tests to enable: '${${PROJECT_NAME}_VALID_CATEGORIES_STR}' (default `${${PROJECT_NAME}_TEST_CATEGORIES_DEFAULT}`)."
    )
  TRIBITS_GET_INVALID_CATEGORIES(${PROJECT_NAME}_TEST_CATEGORIES)

  IF ("${${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT}" STREQUAL "" )
    SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET(
    ${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE
    ${${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT}
    CACHE BOOL
    "Generate a <ProjectName>RepoVersion.txt file.")

  IF ("${DART_TESTING_TIMEOUT_DEFAULT}"  STREQUAL "")
    SET(DART_TESTING_TIMEOUT_DEFAULT  1500)
  ENDIF()
  ADVANCED_SET(
    DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT_DEFAULT}
    CACHE STRING
    "Raw CMake/CTest global default test timeout (default 1500).  (NOTE: Does not impact timeouts of tests that have the TIMEOUT property set on a test-by-test basis.)"
    )
  # NOTE: 1500 is the CMake default set in Modules/CTest.cmake.  We need to
  # set the default here because we need to be able to scale it correctly in
  # case the user does not explicilty set this var in the cache.

  ADVANCED_SET(${PROJECT_NAME}_SCALE_TEST_TIMEOUT 1.0 CACHE STRING
    "Scale factor for global DART_TESTING_TIMEOUT and individual test TIMEOUT (default 1.0)."
    )
  # NOTE: This value is 1.0, *NOT* 1!  This is used in TRIBITS_SCALE_TIMEOUT()
  # and there are unit tests that rely on this default!

  ADVANCED_SET(${PROJECT_NAME}_REL_CPU_SPEED 1.0 CACHE STRING
    "Relative CPU speed of the computer used to scale performance tests (default 1.0)."
    )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}
    CACHE BOOL
    "Determines if a variety of development mode checks are turned on by default or not." )

  ADVANCED_SET( ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL
    "Determines if asserts are performed on missing packages or not." )

  ADVANCED_SET( ${PROJECT_NAME}_WARN_ABOUT_MISSING_EXTERNAL_PACKAGES
    FALSE  CACHE  BOOL
    "If set to TRUE, a 'NOTE' is printed for each missing package that is ignored." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C code for supported compilers." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C++ code for supported compilers." )

  MULTILINE_SET( ENABLE_SHADOW_WARNINGS_DOC
    "Turn ON or OFF shadowing warnings for all packages where strong warnings have"
    " not been explicitly disabled.  Setting the empty '' let's each package decide." )
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS ""
    "${ENABLE_SHADOW_WARNINGS_DOC}" )
  MARK_AS_ADVANCED(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_COVERAGE_TESTING OFF
    CACHE BOOL "Enable support for coverage testing by setting needed compiler/linker options." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_CHECKED_STL OFF
    CACHE BOOL "Turn on checked STL checking (e.g. -D_GLIBCXX_DEBUG) or not." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS OFF
    CACHE BOOL "Turn on debugging symbols (e.g. -g) or not if not a full debug build." )

  IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "-Werror")
  ELSE()
    SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "")
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT}"
    CACHE STRING "Flags for treating warnings as errors (for all compilers, -Werror by default for GNU).  To turn off warnings as errors set to ''")

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF CACHE BOOL
    "If test output complaining about circular references is found, then the test will fail." )

  ADVANCED_SET(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR ""
    CACHE FILEPATH
    "If set to non-null, this is the default directory where package dependency files will be written.")

  IF (${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR)
    SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}")
  ELSE()
    SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT "")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
    "${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "Output XML file containing ${PROJECT_NAME} dependenices used by tools (if not empty)." )

  IF(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR AND
    ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE
    )
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
  ELSE()
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT "")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
    "${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "Output XML file used by CDash in ${PROJECT_NAME}-independent format (if not empty)." )

  IF(${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR AND
    ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE
    )
    SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT
      "${${PROJECT_NAME}_DEPS_DEFAULT_OUTPUT_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME}" )
  ELSE()
    SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT "")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
    "${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "HTML ${PROJECT_NAME} dependenices file that will be written to (if not empty)." )

  #
  # Extra repositories
  #

  ASSERT_DEFINED(${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME)

  SET(DEFAULT_EXTRA_REPOS_FILE
    "${PROJECT_SOURCE_DIR}/cmake/${${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME}")

  IF (EXISTS ${DEFAULT_EXTRA_REPOS_FILE})
    #MESSAGE("${DEFAULT_EXTRA_REPOS_FILE} does exist!")
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT ${DEFAULT_EXTRA_REPOS_FILE})
  ELSE()
    #MESSAGE("${DEFAULT_EXTRA_REPOS_FILE} does *NOT* exist!")
    SET(${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT)
  ENDIF()

  ADVANCED_SET(${PROJECT_NAME}_EXTRAREPOS_FILE
    "${${PROJECT_NAME}_EXTRAREPOS_FILE_DEFAULT}"
    CACHE FILENAME
    "File containing the list of extra repositories containing add-on packages to process")
  #PRINT_VAR(${PROJECT_NAME}_EXTRAREPOS_FILE)

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
    ""
    CACHE STRING
    "Type of testing to pull in extra repositories (Continuous, or Nightly)" )

  ADVANCED_SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES
    FALSE CACHE BOOL
   "Set if to ignore missing extra repositories (or fail hard)" )

  # Even if a project does not support an extra repos file, it can always
  # support extra repositories defined by the user by the very nature of
  # Tribits.

  ADVANCED_SET(${PROJECT_NAME}_PRE_REPOSITORIES
    ""
    CACHE STRING
    "List of pre-extra repositories that contain extra ${PROJECT_NAME} packages."
    )
  SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  "," ${PROJECT_NAME}_PRE_REPOSITORIES)

  ADVANCED_SET(${PROJECT_NAME}_EXTRA_REPOSITORIES
    ""
    CACHE STRING
    "List of post-extra repositories that contain extra ${PROJECT_NAME} packages."
    )
  SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  "," ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  IF ("${${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST}"  STREQUAL  "")
    SET(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST  TRUE)
  ENDIF()

  TRIBITS_GET_AND_PROCESS_EXTRA_REPOSITORIES_LISTS()

  ADVANCED_SET(${PROJECT_NAME}_INSTALLATION_DIR
    ""  CACHE  STRING
    "Location of an installed version of ${PROJECT_NAME} that will be built against during installation testing"
    )

  #
  # More options
  #

  IF("${${PROJECT_NAME}_INSTALLATION_DIR}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT OFF)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT ON)
  ENDIF()

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    ${${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT}
    CACHE STRING
    "Enable testing against an installed version of ${PROJECT_NAME}."
    )

  ADVANCED_OPTION(${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING
    "Short-circuit after dependency handling is complete"
    OFF )

  ADVANCED_OPTION(${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY
    "Only trace dependency handling.  Don't configure to build anything!"
    OFF )

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
    FALSE CACHE BOOL
   "Set to 'ON' to see configure times (Unix/Linux systems only)" )

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
    FALSE CACHE BOOL
   "Set to 'ON' to see configure times for individual packages" )

  IF ("${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT}"
    STREQUAL ""
    )
    SET(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME
    ${${PROJECT_NAME}_SHOW_TEST_START_END_DATE_TIME_DEFAULT}
    CACHE BOOL
    "Set to 'ON' to see start and end date/time for advanced tests." )

  IF ("${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT}"
    STREQUAL ""
    )
    SET(${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT OFF)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST
    ${${PROJECT_NAME}_SHOW_MACHINE_LOAD_IN_TEST_DEFAULT}
    CACHE BOOL
    "Set to 'ON' to see the machine load for advanced tests." )

  MARK_AS_ADVANCED(BUILD_TESTING)
  MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
  MARK_AS_ADVANCED(DART_TESTING_TIMEOUT)
  MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
  MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
  MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)

ENDMACRO()


MACRO(TRIBITS_SETUP_INSTALLATION_PATHS)

  #
  # A) Determine if we are going to be using default paths from GNUInstallDirs module
  #

  SET(TRIBITS_USE_GNUINSTALLDIRS TRUE)

  IF (NOT ${PROJECT_NAME}_USE_GNUINSTALLDIRS)
    # For backward compatibility and unit testing
    SET(TRIBITS_USE_GNUINSTALLDIRS FALSE)
  ENDIF()

  #
  # B) Pick the defaults for the install dirs
  #

  IF (TRIBITS_USE_GNUINSTALLDIRS)
    INCLUDE(GNUInstallDirs)
    SET(${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT ${CMAKE_INSTALL_INCLUDEDIR})
    SET(${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT ${CMAKE_INSTALL_LIBDIR})
    SET(${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT ${CMAKE_INSTALL_BINDIR})
    SET(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT "example")
  ELSE()
    SET(${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT "include")
    SET(${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT "lib")
    SET(${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT "bin")
    SET(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT "example")
  ENDIF()

  #
  # C) Set the cache variables for the install dirs
  #

  ADVANCED_SET( ${PROJECT_NAME}_INSTALL_INCLUDE_DIR
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR_DEFAULT}
    CACHE PATH
    "Location where the headers will be installed.  If given as a STRING type and relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'include'"
    )

  ADVANCED_SET( ${PROJECT_NAME}_INSTALL_LIB_DIR
    ${${PROJECT_NAME}_INSTALL_LIB_DIR_DEFAULT}
    CACHE PATH
    "Location where the libraries will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'lib'"
    )

  ADVANCED_SET( ${PROJECT_NAME}_INSTALL_RUNTIME_DIR
    ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR_DEFAULT}
    CACHE PATH
    "Location where the runtime DLLs and designated programs will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'bin'"
    )

  ADVANCED_SET(${PROJECT_NAME}_INSTALL_EXAMPLE_DIR
    ${${PROJECT_NAME}_INSTALL_EXAMPLE_DIR_DEFAULT}
    CACHE PATH
    "Location where assorted examples will be installed.  If given as a STRING type relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'example'"
    )

  #
  # D) Setup RPATH handling
  #

  PRINT_VAR(${PROJECT_NAME}_SET_INSTALL_RPATH)
  PRINT_VAR(CMAKE_INSTALL_RPATH_USE_LINK_PATH)

  IF (${PROJECT_NAME}_SET_INSTALL_RPATH)
    IF ("${CMAKE_INSTALL_RPATH}" STREQUAL "")
      MESSAGE("-- " "Setting default for CMAKE_INSTALL_RPATH pointing to ${PROJECT_NAME}_INSTALL_LIB_DIR")
      ASSERT_DEFINED(CMAKE_INSTALL_PREFIX)
      ASSERT_DEFINED(${PROJECT_NAME}_INSTALL_LIB_DIR)
      IF (IS_ABSOLUTE ${${PROJECT_NAME}_INSTALL_LIB_DIR})
        SET(CMAKE_INSTALL_RPATH
          "${PROJECT_NAME}_INSTALL_LIB_DIR}" )
      ELSE()
        SET(CMAKE_INSTALL_RPATH
          "${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}" )
      ENDIF()
    ENDIF()
    IF (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
      IF ("${CMAKE_MACOSX_RPATH}" STREQUAL "")
        MESSAGE("-- " "Setting default CMAKE_MACOSX_RPATH=TRUE")
        SET(CMAKE_MACOSX_RPATH TRUE)
      ENDIF()
      PRINT_VAR(CMAKE_MACOSX_RPATH)
    ENDIF()
  ENDIF()
  STRING(REPLACE ":" ";" CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}")
  PRINT_VAR(CMAKE_INSTALL_RPATH)

ENDMACRO()


#
# Repository specializaiton call-back functions
#
# NOTE: The Tribits system promises to only include these call-back files once
# (in order) and to only the call call-back macros they provide once (in
# order).
#


MACRO(CREATE_EMPTY_TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)
  MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)
  ENDMACRO()
ENDMACRO()


MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS_RUNNER  REPO_NAME)
  SET(CALLBACK_SETUP_EXTRA_OPTIONS_FILE
    "${${REPO_NAME}_SOURCE_DIR}/cmake/CallbackSetupExtraOptions.cmake")
  #PRINT_VAR(CALLBACK_SETUP_EXTRA_OPTIONS_FILE)
  IF (EXISTS ${CALLBACK_SETUP_EXTRA_OPTIONS_FILE})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Processing call-back file and macros in"
        " '${CALLBACK_SETUP_EXTRA_OPTIONS_FILE}'")
    ENDIF()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    CREATE_EMPTY_TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()
    # Include the file which will define the callback macros
    SET(REPOSITORY_NAME ${REPO_NAME})
    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${CALLBACK_SETUP_EXTRA_OPTIONS_FILE}")
    INCLUDE(${CALLBACK_SETUP_EXTRA_OPTIONS_FILE})
    # Call the callback macros to inject repository-specific behavir
    TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()
    # Set back the callback macros to empty to ensure that nonone calls them
    CREATE_EMPTY_TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS()
  ENDIF()
ENDMACRO()


MACRO(CREATE_EMPTY_TRIBITS_REPOSITORY_DEFINE_PACKAGING)
  MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)
  ENDMACRO()
ENDMACRO()


MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING_RUNNER  REPO_NAME)
  SET(CALLBACK_DEFINE_PACKAGING_FILE
    "${${REPO_NAME}_SOURCE_DIR}/cmake/CallbackDefineRepositoryPackaging.cmake")
  #PRINT_VAR(CALLBACK_DEFINE_PACKAGING_FILE)
  IF (EXISTS ${CALLBACK_DEFINE_PACKAGING_FILE})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Processing call-back file and macros in"
        " '${CALLBACK_DEFINE_PACKAGING_FILE}'")
    ENDIF()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    CREATE_EMPTY_TRIBITS_REPOSITORY_DEFINE_PACKAGING()
    # Include the file which will define the callback macros
    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${CALLBACK_DEFINE_PACKAGING_FILE}")
    INCLUDE(${CALLBACK_DEFINE_PACKAGING_FILE})
    # Call the callback macros to inject repository-specific behavir
    TRIBITS_REPOSITORY_DEFINE_PACKAGING()
    # Set back the callback macros to empty to ensure that nonone calls them
    CREATE_EMPTY_TRIBITS_REPOSITORY_DEFINE_PACKAGING()
  ENDIF()
ENDMACRO()


MACRO(CREATE_EMPTY_TRIBITS_PROJECT_DEFINE_PACKAGING)
  MACRO(TRIBITS_PROJECT_DEFINE_PACKAGING)
  ENDMACRO()
ENDMACRO()


MACRO(TRIBITS_PROJECT_DEFINE_PACKAGING_RUNNER)
  SET(CALLBACK_DEFINE_PACKAGING_FILE
    "${PROJECT_SOURCE_DIR}/cmake/CallbackDefineProjectPackaging.cmake")
  #PRINT_VAR(CALLBACK_DEFINE_PACKAGING_FILE)
  IF (EXISTS ${CALLBACK_DEFINE_PACKAGING_FILE})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Processing call-back file and macros in"
        " '${CALLBACK_DEFINE_PACKAGING_FILE}'")
    ENDIF()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    CREATE_EMPTY_TRIBITS_PROJECT_DEFINE_PACKAGING()
    # Include the file which will define the callback macros
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE
      "${CALLBACK_DEFINE_PACKAGING_FILE}")
    INCLUDE(${CALLBACK_DEFINE_PACKAGING_FILE})
    # Call the callback macros to inject project-specific behavir
    TRIBITS_PROJECT_DEFINE_PACKAGING()
    # Set back the callback macros to empty to ensure that nonone calls them
    CREATE_EMPTY_TRIBITS_PROJECT_DEFINE_PACKAGING()
  ENDIF()
ENDMACRO()


#
# Read in the Project's native repositories.,
#
# On output, the variable ${PRJOECT_NAME}_NATIVE_REPOSITORIES is set.
#
MACRO(TRIBITS_READ_IN_NATIVE_REPOSITORIES)
  IF (${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE)
    IF (IS_ABSOLUTE ${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
      SET(NATIVE_REPO_FILE ${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
    ELSE()
      SET(NATIVE_REPO_FILE
        ${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_NATIVE_REPO_FILE_OVERRRIDE})
    ENDIF()
  ELSE()
    SET(NATIVE_REPO_FILE ${PROJECT_SOURCE_DIR}/cmake/NativeRepositoriesList.cmake)
  ENDIF()
  IF (EXISTS ${NATIVE_REPO_FILE})
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${NATIVE_REPO_FILE}")
    INCLUDE(${NATIVE_REPO_FILE})
  ELSE()
    SET(${PROJECT_NAME}_NATIVE_REPOSITORIES ".")
  ENDIF()
ENDMACRO()


#
# Combine native and extra repos lists into a single list.
#
# Combines ${PROJECT_NAME}_PRE_REPOSITORIES
# ${PROJECT_NAME}_NATIVE_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES
# into a single list ${PROJECT_NAME}_ALL_REPOSITORIES.
#
MACRO(TRIBITS_COMBINE_NATIVE_AND_EXTRA_REPOS)
  ASSERT_DEFINED(${PROJECT_NAME}_PRE_REPOSITORIES)
  ASSERT_DEFINED(${PROJECT_NAME}_NATIVE_REPOSITORIES)
  ASSERT_DEFINED(${PROJECT_NAME}_EXTRA_REPOSITORIES)
  SET( ${PROJECT_NAME}_ALL_REPOSITORIES
    ${${PROJECT_NAME}_PRE_REPOSITORIES}
    ${${PROJECT_NAME}_NATIVE_REPOSITORIES}
    ${${PROJECT_NAME}_EXTRA_REPOSITORIES}
    )
ENDMACRO()


#
# Process extra repo extra options files
#
MACRO(TRIBITS_PROCESS_EXTRA_REPOS_OPTIONS_FILES)
  # Loop through the Repositories, set their base directories and run their
  # options setup callback functions.
  FOREACH(REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    TRIBITS_GET_REPO_NAME_DIR(${REPO}  REPO_NAME  REPO_DIR)
    TRIBITS_SET_BASE_REPO_DIR(${PROJECT_SOURCE_DIR}  ${REPO_DIR}  ${REPO_NAME}_SOURCE_DIR)
    TRIBITS_SET_BASE_REPO_DIR(${PROJECT_BINARY_DIR}  ${REPO_DIR}  ${REPO_NAME}_BINARY_DIR)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Processing extra options call-backs for ${REPO}")
      PRINT_VAR(${REPO_NAME}_SOURCE_DIR)
      PRINT_VAR(${REPO_NAME}_BINARY_DIR)
    ENDIF()
    TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS_RUNNER(${REPO_NAME})
  ENDFOREACH()
ENDMACRO()


#
# Copy an simple text file to the binary dir to be included in the tarball
#

MACRO(TRIBITS_COPY_INSTALLER_RESOURCE _varname _source _destination)
  SET("${_varname}" "${_destination}")
  IF (EXISTS "${_destination}")
    FILE(REMOVE_RECURSE "${_destination}")
  ENDIF ()
  CONFIGURE_FILE(
    "${_source}"
    "${_destination}"
    COPYONLY)
ENDMACRO()

#
# Run the git log command to get the verison info for a git rep
#

FUNCTION(TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING  GIT_REPO_DIR
   SINGLE_REPO_VERSION_STRING_OUT
  )

  IF (NOT GIT_EXEC)
    MESSAGE(SEND_ERROR "ERROR, the program '${GIT_NAME}' could not be found!"
      "  We can not generate the repo version file!")
  ENDIF()

  # A) Get the basic version info.

  EXECUTE_PROCESS(
    COMMAND ${GIT_EXEC} log -1 --pretty=format:"%h [%ad] <%ae>"
    WORKING_DIRECTORY ${GIT_REPO_DIR}
    RESULT_VARIABLE GIT_RETURN
    OUTPUT_VARIABLE GIT_OUTPUT
    )
  # NOTE: Above we have to add quotes '"' or CMake will not accept the
  # command.  However, git will put those quotes in the output so we have to
  # strip them out later :-(

  IF (NOT GIT_RETURN STREQUAL 0)
    MESSAGE(FATAL_ERROR "ERROR, ${GIT_EXEC} command returned ${GIT_RETURN}!=0"
      " for extra repo ${GIT_REPO_DIR}!")
    SET(GIT_VERSION_INFO "Error, could not get version info!")
  ELSE()
    # Strip the quotes off :-(
    STRING(LENGTH "${GIT_OUTPUT}" GIT_OUTPUT_LEN)
    MATH(EXPR OUTPUT_NUM_CHARS_TO_KEEP "${GIT_OUTPUT_LEN}-2")
    STRING(SUBSTRING "${GIT_OUTPUT}" 1 ${OUTPUT_NUM_CHARS_TO_KEEP}
      GIT_VERSION_INFO)
  ENDIF()

  # B) Get the first 80 chars of the summary message for more info

  EXECUTE_PROCESS(
    COMMAND ${GIT_EXEC} log -1 --pretty=format:"%s"
    WORKING_DIRECTORY ${GIT_REPO_DIR}
    RESULT_VARIABLE GIT_RETURN
    OUTPUT_VARIABLE GIT_OUTPUT
    )

  IF (NOT GIT_RETURN STREQUAL 0)
    MESSAGE(FATAL_ERROR "ERROR, ${GIT_EXEC} command returned ${GIT_RETURN}!=0"
      " for extra repo ${GIT_REPO_DIR}!")
    SET(GIT_VERSION_SUMMARY "Error, could not get version summary!")
  ELSE()
    # Strip ouf quotes and quote the 80 char string
    SET(MAX_SUMMARY_LEN 80)
    MATH(EXPR MAX_SUMMARY_LEN_PLUS_2 "${MAX_SUMMARY_LEN}+2")
    STRING(LENGTH "${GIT_OUTPUT}" GIT_OUTPUT_LEN)
    MATH(EXPR OUTPUT_NUM_CHARS_TO_KEEP "${GIT_OUTPUT_LEN}-2")
    STRING(SUBSTRING "${GIT_OUTPUT}" 1 ${OUTPUT_NUM_CHARS_TO_KEEP}
      GIT_OUTPUT_STRIPPED)
    IF (GIT_OUTPUT_LEN GREATER ${MAX_SUMMARY_LEN_PLUS_2})
      STRING(SUBSTRING "${GIT_OUTPUT_STRIPPED}" 0 ${MAX_SUMMARY_LEN}
         GIT_SUMMARY_STR)
    ELSE()
      SET(GIT_SUMMARY_STR "${GIT_OUTPUT_STRIPPED}")
    ENDIF()
  ENDIF()

  SET(${SINGLE_REPO_VERSION_STRING_OUT}
    "${GIT_VERSION_INFO}\n${GIT_SUMMARY_STR}" PARENT_SCOPE)

ENDFUNCTION()


#
# Get the versions of all the git repos
#

FUNCTION(TRIBITS_GENERATE_REPO_VERSION_FILE_STRING  PROJECT_REPO_VERSION_FILE_STRING_OUT)

  SET(REPO_VERSION_FILE_STR "")

  TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(
     ${CMAKE_CURRENT_SOURCE_DIR}
     SINGLE_REPO_VERSION)
  APPEND_STRING_VAR(REPO_VERSION_FILE_STR
    "*** Base Git Repo: ${PROJECT_NAME}\n"
    "${SINGLE_REPO_VERSION}\n" )

  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRA_REPO ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

    #PRINT_VAR(EXTRA_REPO)
    #PRINT_VAR(EXTRAREPO_IDX)
    #PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)

    IF (${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS)
      # Read from an extra repo file with potentially different dir.
      LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
        EXTRAREPO_DIR )
    ELSE()
       # Not read from extra repo file so dir is same as name
       SET(EXTRAREPO_DIR ${EXTRA_REPO})
    ENDIF()
    #PRINT_VAR(EXTRAREPO_DIR)

    TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(
       "${CMAKE_CURRENT_SOURCE_DIR}/${EXTRAREPO_DIR}"
       SINGLE_REPO_VERSION)
    APPEND_STRING_VAR(REPO_VERSION_FILE_STR
      "*** Git Repo: ${EXTRAREPO_DIR}\n"
      "${SINGLE_REPO_VERSION}\n" )

    #PRINT_VAR(REPO_VERSION_FILE_STR)

    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  ENDFOREACH()

  SET(${PROJECT_REPO_VERSION_FILE_STRING_OUT} ${REPO_VERSION_FILE_STR} PARENT_SCOPE)

ENDFUNCTION()


#
# Generate the project repos version file and print to stdout
#
# This function is designed so that it can be unit tested from inside of a
# cmake -P script.
#

FUNCTION(TRIBITS_GENERATE_REPO_VERSION_OUTPUT_AND_FILE)
  # Get the repos versions
  TRIBITS_GENERATE_REPO_VERSION_FILE_STRING(PROJECT_REPO_VERSION_FILE_STRING)
  # Print the versions
  MESSAGE("\n${PROJECT_NAME} repos versions:\n"
    "--------------------------------------------------------------------------------\n"
    "${PROJECT_REPO_VERSION_FILE_STRING}"
    " --------------------------------------------------------------------------------\n"
    )
  #) Write out the version file
  FILE(WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
    "${PROJECT_REPO_VERSION_FILE_STRING}")
ENDFUNCTION()


#
# Create project dependencies file and create install target
#
# NOTE: Before calling this function, the extra repos datastructure must be
# filled out!
#
# NOTE: This function can not be called in a cmake -P script because it has a
# call to INSTALL()!  That is why this function is seprated out from
# TRIBITS_GENERATE_REPO_VERSION_OUTPUT_AND_FILE().
#

FUNCTION(TRIBITS_GENERATE_REPO_VERSION_OUTPUT_AND_FILE_AND_INSTALL)

  #
  # A) Create the ${PROJECT_NAME}RepoVersion.txt file if requested
  #

  IF (${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE)

    # A) Make sure that there is a .git dir in the project before generating
    IF (EXISTS "${PROJECT_SOURCE_DIR}/.git")
      SET(PROJECT_SOURCE_IS_GIT_REPO TRUE)
    ELSE()
      SET(PROJECT_SOURCE_IS_GIT_REPO FALSE)
    ENDIF()
    IF (PROJECT_SOURCE_IS_GIT_REPO)
      # Find git first here so we  don't have to find it in called function so
      # it can be unit tested.
      FIND_PROGRAM(GIT_EXEC ${GIT_NAME})
      # Get repo versions, print to stdout and write file
      TRIBITS_GENERATE_REPO_VERSION_OUTPUT_AND_FILE()
      # Add install target for this file
      INSTALL(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
        DESTINATION "." )
    ELSE()
      MESSAGE("\nNOTE: Skipping generation of ${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}"
        " because project source is not a git repo!")
    ENDIF()

    # B) Install the repo version file if it is in source tree which it will
    # be for a tarball (see TRIBITS_SETUP_PACKAGING_AND_DISTRIBUTION()).
    SET(REPO_VERSION_FILE_IN_SOURCE_TREE
      ${CMAKE_CURRENT_SOURCE_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME})
    IF (EXISTS ${REPO_VERSION_FILE_IN_SOURCE_TREE})
      INSTALL(
        FILES "${REPO_VERSION_FILE_IN_SOURCE_TREE}"
        DESTINATION "." )
    ENDIF()

  ENDIF()


ENDFUNCTION()


#
# Sets ${PROJECT_NAME}_EXTRA_REPOSITORIES from
# ${PROJECT_NAME}_EXTRA_REPOSITORIES and ${PROJECT_NAME}_EXTRA_REPOSITORIES if
# it is not alrady set.  Also, it replaces ',' with ';' in the latter.
#
# This function is needed in use cases where extra repos are used where the
# extra repos are not read in through an ExtraRepositoriesList.cmake file and
# instead are directly passed in by the user.
#
MACRO(TRIBITS_SET_ALL_EXTRA_REPOSITORIES)
  IF ("${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES}"   STREQUAL  "")
    # Allow list to be seprated by ',' instead of just by ';'.  This is needed
    # by the unit test driver code
    SPLIT("${${PROJECT_NAME}_PRE_REPOSITORIES}"  ","
      ${PROJECT_NAME}_PRE_REPOSITORIES)
    SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","
      ${PROJECT_NAME}_EXTRA_REPOSITORIES)
    SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
      ${${PROJECT_NAME}_PRE_REPOSITORIES}  ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
  ENDIF()
ENDMACRO()


#
# Macro that processes the list of package and TPLs for the set of 'PRE' or
# 'POST' extra repos.
#
MACRO(TRIBITS_READ_EXTRA_REPOSITORIES_LISTS)

  LIST(LENGTH  ${PROJECT_NAME}_PRE_REPOSITORIES  PRE_EXTRAREPOS_LEN)
  LIST(LENGTH  ${PROJECT_NAME}_EXTRA_REPOSITORIES  POST_EXTRAREPOS_LEN)
  MATH(EXPR  ALL_EXTRAREPOS_LEN  "${PRE_EXTRAREPOS_LEN} + ${POST_EXTRAREPOS_LEN}")

  # See if processing 'PRE' or 'POST' extra repos
  IF (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "PRE")
    SET(EXTRAREPO_IDX_START  0)
    SET(EXTRAREPO_IDX_END  ${PRE_EXTRAREPOS_LEN})
  ELSEIF (READ_PRE_OR_POST_EXRAREPOS  STREQUAL  "POST")
    SET(EXTRAREPO_IDX_START  ${PRE_EXTRAREPOS_LEN})
    SET(EXTRAREPO_IDX_END  ${ALL_EXTRAREPOS_LEN})
  ELSE()
    MESSAGE(FATAL_ERROR "Invalid value for READ_PRE_OR_POST_EXRAREPOS='${READ_PRE_OR_POST_EXRAREPOS}' ")
  ENDIF()
  # NOTE: For some reason, we can't pass this argument to the function and
  # have it read.  Instead, we have to pass it a local variable.  I will never
  # understand CMake.

  SET(EXTRAREPO_IDX  ${EXTRAREPO_IDX_START})
  WHILE(EXTRAREPO_IDX  LESS  EXTRAREPO_IDX_END)

    LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES  ${EXTRAREPO_IDX}  EXTRA_REPO )

    #PRINT_VAR(EXTRA_REPO)
    #PRINT_VAR(EXTRAREPO_IDX)
    #PRINT_VAR(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    SET(${EXTRA_REPO}_SOURCE_DIR "${PROJECT_SOURCE_DIR}/${EXTRA_REPO}")
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${EXTRA_REPO}_SOURCE_DIR)
    ENDIF()
    # ToDo: TriBITS:73: Get ${EXTRA_REPO}_SOURCE_DIR from
    # ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIR when it exists.

    SET(EXTRAREPO_PACKSTAT "")
    IF (${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS)
      LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
        EXTRAREPO_PACKSTAT )
    ENDIF()

    IF (EXTRAREPO_PACKSTAT STREQUAL NOPACKAGES)

      MESSAGE("")
      MESSAGE("Skipping reading packages and TPLs for ${READ_PRE_OR_POST_EXRAREPOS} extra repo ${EXTRA_REPO} because marked NOPACKAGES ... ")
      MESSAGE("")
      # ToDo: TriBITS:73: Don't print the above message by default.  It is
      # just clutter.

    ELSE()

      # Read in the add-on packages from the extra repo

      #PRINT_VAR(${EXTRA_REPO}_PACKAGES_LIST_FILE)
      IF (${EXTRA_REPO}_PACKAGES_LIST_FILE)
        SET(EXTRAREPO_PACKAGES_FILE
          "${PROJECT_SOURCE_DIR}/${${EXTRA_REPO}_PACKAGES_LIST_FILE}")
      ELSE()
        SET(EXTRAREPO_PACKAGES_FILE
          "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME}")
      ENDIF()

      MESSAGE("")
      MESSAGE("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra packages from ${EXTRAREPO_PACKAGES_FILE} ... ")
      MESSAGE("")

      IF (NOT EXISTS "${EXTRAREPO_PACKAGES_FILE}")
        IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          MESSAGE(
            "\n***"
            "\n*** NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}' on request!"
            "\n***\n")
            # ToDo: TriBITS:73: Shorten above message to just one line
        ELSE()
          MESSAGE( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' packages list file '${EXTRAREPO_PACKAGES_FILE}'!")
          # ToDo: TriBITS:73: Change to FATAL_ERROR to abort early
        ENDIF()
      ELSE()
        SET(REPOSITORY_NAME  ${EXTRA_REPO})
        TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE  "${EXTRAREPO_PACKAGES_FILE}")
        INCLUDE("${EXTRAREPO_PACKAGES_FILE}")
        SET(APPEND_TO_PACKAGES_LIST  TRUE)
        TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${EXTRA_REPO} ${EXTRA_REPO})
      ENDIF()

      # Read in the add-on TPLs from the extra repo

      SET(${EXTRA_REPO}_TPLS_FILE
        "${${EXTRA_REPO}_SOURCE_DIR}/${${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME}")

      MESSAGE("")
      MESSAGE("Reading list of ${READ_PRE_OR_POST_EXRAREPOS} extra TPLs from ${${EXTRA_REPO}_TPLS_FILE} ... ")
      MESSAGE("")

      IF (NOT EXISTS "${${EXTRA_REPO}_TPLS_FILE}")
        IF (${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
          MESSAGE(
            "\n***"
            "\n*** NOTE: Ignoring missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}' on request!"
            "\n***\n")
          # ToDo: TriBITS:73: Shorten above warning to just one line
        ELSE()
          MESSAGE( SEND_ERROR
            "ERROR: Skipping missing ${READ_PRE_OR_POST_EXRAREPOS} extra repo '${EXTRA_REPO}' TPLs list file '${${EXTRA_REPO}_TPLS_FILE}'!")
        ENDIF()
      ELSE()
        TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE  "${${EXTRA_REPO}_TPLS_FILE}")
        INCLUDE("${${EXTRA_REPO}_TPLS_FILE}")
        SET(APPEND_TO_TPLS_LIST  TRUE)
        TRIBITS_PROCESS_TPLS_LISTS(${EXTRA_REPO}  ${EXTRA_REPO})
      ENDIF()

    ENDIF()

    MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

  ENDWHILE()

ENDMACRO()

#
# Read in ${PROJECT_NAME} packages and TPLs, process dependencies, write XML
# files
#
MACRO(TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML)

  TRIBITS_SET_ALL_EXTRA_REPOSITORIES()

  # Set to empty
  SET(${PROJECT_NAME}_PACKAGES)
  SET(${PROJECT_NAME}_PACKAGE_DIRS)
  SET(${PROJECT_NAME}_TPLS)

  #
  # A) Read list of packages and TPLs from 'PRE' extra repos
  #

  SET(READ_PRE_OR_POST_EXRAREPOS  PRE)
  TRIBITS_READ_EXTRA_REPOSITORIES_LISTS()

  #
  # B) Read list of packages and TPLs from native repos
  #

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SET_UP_DEPENDENCIES_TIME_START_SECONDS)
  ENDIF()

  FOREACH(NATIVE_REPO ${${PROJECT_NAME}_NATIVE_REPOSITORIES})

    TRIBITS_GET_REPO_NAME_DIR(${NATIVE_REPO}  NATIVE_REPO_NAME  NATIVE_REPO_DIR)
    #PRINT_VAR(NATIVE_REPO_NAME)
    #PRINT_VAR(NATIVE_REPO_DIR)

    # Need to make sure this gets set because logic in Dependencies.cmake files
    # looks for the presents of this variable.
    TRIBITS_SET_BASE_REPO_DIR(${PROJECT_SOURCE_DIR} ${NATIVE_REPO_DIR}
      ${NATIVE_REPO_NAME}_SOURCE_DIR)
    #PRINT_VAR(${NATIVE_REPO_NAME}_SOURCE_DIR)

    #
    # B.1) Define the lists of all ${NATIVE_REPO_NAME} native packages and TPLs
    #

    # B.1.a) Read the core ${NATIVE_REPO_NAME} packages
    IF (${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE)
      IF (IS_ABSOLUTE "${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
        MESSAGE(FATAL_ERROR
          "ToDo: Implement abs path for ${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE")
      ELSE()
        SET(${NATIVE_REPO_NAME}_PACKAGES_FILE
          "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${NATIVE_REPO_NAME}_PACKAGES_FILE_OVERRIDE}")
      ENDIF()
    ELSE()
      SET(${NATIVE_REPO_NAME}_PACKAGES_FILE
        "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_PACKAGES_FILE_NAME}")
    ENDIF()

    MESSAGE("")
    MESSAGE("Reading list of native packages from ${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    MESSAGE("")

    IF (NATIVE_REPO STREQUAL ".")
      SET(REPOSITORY_NAME ${PROJECT_NAME})
    ELSE()
      SET(REPOSITORY_NAME ${NATIVE_REPO_NAME})
    ENDIF()
    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_PACKAGES_FILE}")
    INCLUDE(${${NATIVE_REPO_NAME}_PACKAGES_FILE})

    TRIBITS_PROCESS_PACKAGES_AND_DIRS_LISTS(${NATIVE_REPO_NAME} ${NATIVE_REPO_DIR})

    # B.1.b) Read the core TPLs dependencies

    SET(${NATIVE_REPO_NAME}_TPLS_FILE
      "${${NATIVE_REPO_NAME}_SOURCE_DIR}/${${PROJECT_NAME}_TPLS_FILE_NAME}")

    MESSAGE("")
    MESSAGE("Reading list of native TPLs from ${${NATIVE_REPO_NAME}_TPLS_FILE}")
    MESSAGE("")

    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE
      "${${NATIVE_REPO_NAME}_TPLS_FILE}")
    INCLUDE(${${NATIVE_REPO_NAME}_TPLS_FILE})
    TRIBITS_PROCESS_TPLS_LISTS(${NATIVE_REPO_NAME}  ${NATIVE_REPO_DIR})

  ENDFOREACH()

  #
  # C) Read list of packages and TPLs from 'POST' extra repos
  #

  SET(READ_PRE_OR_POST_EXRAREPOS  POST)
  TRIBITS_READ_EXTRA_REPOSITORIES_LISTS()

  #
  # D) Process lists of packages, TPLs, etc.
  #

  #
  # D.1) Package dependencies for all of the packages for all of the defined
  # packages (not just the core packages)
  #

  TRIBITS_READ_ALL_PACKAGE_DEPENDENCIES()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SET_UP_DEPENDENCIES_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${SET_UP_DEPENDENCIES_TIME_START_SECONDS}
      ${SET_UP_DEPENDENCIES_TIME_STOP_SECONDS}
      "\nTotal time to read in and process all package dependencies")
  ENDIF()

  #
  # D.2) Write out the XML dependency files for the full list of dependencies!
  #

  SET(TRIBITS_DUMP_XML_DEPS_MODULE
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CI_SUPPORT_DIR}/TribitsDumpXmlDependenciesFiles.cmake
    )
  IF (EXISTS ${TRIBITS_DUMP_XML_DEPS_MODULE})
    INCLUDE(${TRIBITS_DUMP_XML_DEPS_MODULE})
    TRIBITS_WRITE_XML_DEPENDENCY_FILES()
  ENDIF()

ENDMACRO()


#
# Print out a list with white-space separators with an initial doc string
#
FUNCTION(TRIBITS_PRINT_PREFIX_STRING_AND_LIST  DOCSTRING   LIST_TO_PRINT)
  STRING(REPLACE ";" " " LIST_TO_PRINT_STR "${LIST_TO_PRINT}")
  LIST(LENGTH  LIST_TO_PRINT  NUM_ELEMENTS)
  IF (NUM_ELEMENTS GREATER "0")
    MESSAGE("${DOCSTRING}:  ${LIST_TO_PRINT_STR} ${NUM_ELEMENTS}")
  ELSE()
    MESSAGE("${DOCSTRING}:  ${NUM_ELEMENTS}")
  ENDIF()
ENDFUNCTION()


#
# Function that prints the current set of enabled/disabled packages given
# input list of packages.
#
FUNCTION(TRIBITS_PRINT_ENABLED_PACKAGES_LIST_FROM_VAR  PACKAGES_LIST_VAR
  DOCSTRING  ENABLED_FLAG  INCLUDE_EMPTY
  )
  IF (ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_ENABLED_LIST(${PACKAGES_LIST_VAR}  ${PROJECT_NAME}
      ENABLED_PACKAGES  NUM_ENABLED)
  ELSEIF (ENABLED_FLAG AND INCLUDE_EMPTY)
    TRIBITS_GET_NONDISABLED_LIST(${PACKAGES_LIST_VAR}  ${PROJECT_NAME}
      ENABLED_PACKAGES  NUM_ENABLED)
  ELSEIF (NOT ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_DISABLED_LIST(${PACKAGES_LIST_VAR}  ${PROJECT_NAME}
      ENABLED_PACKAGES  NUM_ENABLED)
  ELSE() # NOT ENABLED_FLAG AND INCLUDE_EMPTY
    TRIBITS_GET_NONENABLED_LIST(${PACKAGES_LIST_VAR}  ${PROJECT_NAME}
      ENABLED_PACKAGES  NUM_ENABLED)
  ENDIF()
  TRIBITS_PRINT_PREFIX_STRING_AND_LIST("${DOCSTRING}"  "${ENABLED_PACKAGES}")
ENDFUNCTION()


#
# Function that prints the current set of enabled/disabled packages
#
FUNCTION(TRIBITS_PRINT_ENABLED_PACKAGE_LIST  DOCSTRING  ENABLED_FLAG  INCLUDE_EMPTY)
  TRIBITS_PRINT_ENABLED_PACKAGES_LIST_FROM_VAR( ${PROJECT_NAME}_PACKAGES
    "${DOCSTRING}" ${ENABLED_FLAG} ${INCLUDE_EMPTY} )
ENDFUNCTION()


#
# Function that prints the current set of enabled/disabled SE packages
#
FUNCTION(TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST  DOCSTRING  ENABLED_FLAG  INCLUDE_EMPTY)
  IF (ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
      ENABLED_SE_PACKAGES  NUM_ENABLED)
  ELSEIF (ENABLED_FLAG AND INCLUDE_EMPTY)
    TRIBITS_GET_NONDISABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
      ENABLED_SE_PACKAGES  NUM_ENABLED)
  ELSEIF (NOT ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_DISABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
      ENABLED_SE_PACKAGES  NUM_ENABLED)
  ELSE() # NOT ENABLED_FLAG AND INCLUDE_EMPTY
    TRIBITS_GET_NONENABLED_LIST( ${PROJECT_NAME}_SE_PACKAGES  ${PROJECT_NAME}
      ENABLED_SE_PACKAGES  NUM_ENABLED)
  ENDIF()
  TRIBITS_PRINT_PREFIX_STRING_AND_LIST("${DOCSTRING}"  "${ENABLED_SE_PACKAGES}")
ENDFUNCTION()


#
# Function that prints the current set of enabled/disabled TPLs
#
FUNCTION(TRIBITS_PRINT_ENABLED_TPL_LIST  DOCSTRING  ENABLED_FLAG  INCLUDE_EMPTY)
  IF (ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_ENABLED_LIST( ${PROJECT_NAME}_TPLS  TPL
      ENABLED_TPLS  NUM_ENABLED)
  ELSEIF (ENABLED_FLAG AND INCLUDE_EMPTY)
    TRIBITS_GET_NONDISABLED_LIST( ${PROJECT_NAME}_TPLS  TPL
      ENABLED_TPLS  NUM_ENABLED)
  ELSEIF (NOT ENABLED_FLAG AND NOT INCLUDE_EMPTY)
    TRIBITS_GET_DISABLED_LIST( ${PROJECT_NAME}_TPLS  TPL
      ENABLED_TPLS  NUM_ENABLED)
  ELSE() # NOT ENABLED_FLAG AND INCLUDE_EMPTY
    TRIBITS_GET_NONENABLED_LIST( ${PROJECT_NAME}_TPLS  TPL
       ENABLED_TPLS  NUM_ENABLED)
  ENDIF()
  TRIBITS_PRINT_PREFIX_STRING_AND_LIST("${DOCSTRING}"  "${ENABLED_TPLS}")
ENDFUNCTION()


#
# Adjust package enable logic and print out before and after state
#
# On output sets:
#
#    ${PROJECT_NAME}_NUM_ENABLED_PACKAGES: Number of enabled packages (local variable)
#    ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}: Enable status of PACKAGE_NAME (local variable)
#    ToDo: Fill in others as well!
#
MACRO(TRIBITS_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(ADJUST_PACKAGE_DEPS_TIME_START_SECONDS)
  ENDIF()

  TRIBITS_PRINT_ENABLED_PACKAGE_LIST(
    "\nExplicitly enabled packages on input (by user)" ON FALSE)
  TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST(
    "\nExplicitly enabled SE packages on input (by user)" ON FALSE)
  TRIBITS_PRINT_ENABLED_PACKAGE_LIST(
    "\nExplicitly disabled packages on input (by user or by default)" OFF FALSE)
  TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST(
    "\nExplicitly disabled SE packages on input (by user or by default)" OFF FALSE)
  TRIBITS_PRINT_ENABLED_TPL_LIST(
    "\nExplicitly enabled TPLs on input (by user)" ON FALSE)
  TRIBITS_PRINT_ENABLED_TPL_LIST(
    "\nExplicitly disabled TPLs on input (by user or by default)" OFF FALSE)

  TRIBITS_ADJUST_PACKAGE_ENABLES()

  TRIBITS_PRINT_PREFIX_STRING_AND_LIST(
    "\nFinal set of enabled packages" "${${PROJECT_NAME}_ENABLED_PACKAGES}")
  TRIBITS_PRINT_PREFIX_STRING_AND_LIST(
    "\nFinal set of enabled SE packages" "${${PROJECT_NAME}_ENABLED_SE_PACKAGES}")
  TRIBITS_PRINT_ENABLED_PACKAGE_LIST(
    "\nFinal set of non-enabled packages" OFF TRUE)
  TRIBITS_PRINT_ENABLED_SE_PACKAGE_LIST(
    "\nFinal set of non-enabled SE packages" OFF TRUE)
  TRIBITS_PRINT_ENABLED_TPL_LIST(
    "\nFinal set of enabled TPLs" ON FALSE)
  TRIBITS_PRINT_ENABLED_TPL_LIST(
    "\nFinal set of non-enabled TPLs" OFF TRUE)

  TRIBITS_SET_UP_ENABLED_ONLY_DEPENDENCIES()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(ADJUST_PACKAGE_DEPS_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${ADJUST_PACKAGE_DEPS_TIME_START_SECONDS}
      ${ADJUST_PACKAGE_DEPS_TIME_STOP_SECONDS}
      "\nTotal time to adjust package and TPL enables")
  ENDIF()

ENDMACRO()


#
# Macro that gathers information from enabled TPLs
#

MACRO(TRIBITS_PROCESS_ENABLED_TPLS)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(CONFIGURE_TPLS_TIME_START_SECONDS)
  ENDIF()

  FOREACH(TPL_NAME ${${PROJECT_NAME}_TPLS})
    IF (TPL_ENABLE_${TPL_NAME})
      TRIBITS_PROCESS_ENABLED_TPL(${TPL_NAME})
    ENDIF()
  ENDFOREACH()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(CONFIGURE_TPLS_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${CONFIGURE_TPLS_TIME_START_SECONDS}
      ${CONFIGURE_TPLS_TIME_STOP_SECONDS}
      "\nTotal time to configure enabled TPLs")
  ENDIF()

ENDMACRO()


#
# Macros for setting up the standard environment
#


MACRO(TRIBITS_SETUP_ENV)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SETUP_ENV_TIME_START_SECONDS)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    SET(TRIBITS_SETUP_ENV_DEBUG  TRUE)
  ENDIF()

  # Apply XSDK defaults

  IF ("${${PROJECT_NAME}_TRIBITS_XSDK_DIR}" STREQUAL "")
    SET(${PROJECT_NAME}_TRIBITS_XSDK_DIR  "${${PROJECT_NAME}_TRIBITS_DIR}/xsdk")
  ENDIF()
  IF (EXISTS "${${PROJECT_NAME}_TRIBITS_XSDK_DIR}")
    SET(USE_XSDK_DEFAULTS_DEFAULT  FALSE) # Set to TRUE for Trilinos 13.0.0?
    SET(XSDK_ENABLE_C  ${${PROJECT_NAME}_ENABLE_C})
    SET(XSDK_ENABLE_CXX  ${${PROJECT_NAME}_ENABLE_CXX})
    SET(XSDK_ENABLE_Fortran  ${${PROJECT_NAME}_ENABLE_Fortran})
    INCLUDE("${${PROJECT_NAME}_TRIBITS_XSDK_DIR}/XSDKDefaults.cmake")
    # NOTE: BUILD_SHARED_LIBS was set in
    # TRIBITS_DEFINE_GLOBAL_OPTIONS_AND_DEFINE_EXTRA_REPOS() based on
    # USE_XSDK_DEFAULTS in case there is logic in TriBITS that depends on this
    # var getting set there.
  ENDIF()

  # BUILD_SHARED_LIBS
  PRINT_VAR(BUILD_SHARED_LIBS)

  # Set to release build by default

  IF ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    MESSAGE(STATUS "Setting CMAKE_BUILD_TYPE=RELEASE since it was not set ...")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Type of build to perform (i.e. DEBUG, RELEASE, NONE)" )
  ELSE()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UP)
    LIST(FIND CMAKE_BUILD_TYPES_LIST ${CMAKE_BUILD_TYPE_UP} BUILD_TYPE_IDX)
    IF (BUILD_TYPE_IDX EQUAL -1)
      MESSAGE(SEND_ERROR "Error, the given CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        " is not in the list of valid values \"${CMAKE_BUILD_TYPES_LIST}\"!")
    ENDIF()
  ENDIF()
  PRINT_VAR(CMAKE_BUILD_TYPE)

  # Override the silly CMAKE_CONFIGURATION_TYPES variable.  This is needed for
  # MSVS!  Later, we Override CMAKE_CONFIGURATION_TYPES to just one
  # configuration after the compiler checks (see below).
  IF (CMAKE_CONFIGURATION_TYPES)
    IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      SET(CMAKE_CONFIGURATION_TYPE "Debug")
    ELSEIF(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
      SET(CMAKE_CONFIGURATION_TYPE "Release")
    ELSE()
      SET(CMAKE_CONFIGURATION_TYPE "Release")
    ENDIF()
  ELSE()
    SET(CMAKE_CONFIGURATION_TYPE "")
  ENDIF()
  IF (TRIBITS_SETUP_ENV_DEBUG)
    PRINT_VAR(CMAKE_CONFIGURATION_TYPE)
  ENDIF()

  # Set up MPI if MPI is being used

  IF ("${TPL_ENABLE_MPI}" STREQUAL "")
    # If TPL_ENABLE_MPI is undefined or empty because this project does not
    # define an MPI TPL, then explicitly disable it.
    SET(TPL_ENABLE_MPI FALSE)
  ENDIF()

  IF (TPL_ENABLE_MPI)
    TRIBITS_SETUP_MPI()
  ENDIF()

  # Enable compilers

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C)
  IF (${PROJECT_NAME}_ENABLE_C)
    ENABLE_LANGUAGE(C)
    INCLUDE(CMakeDetermineCCompiler)
    PRINT_VAR(CMAKE_C_COMPILER_ID)
    PRINT_VAR(CMAKE_C_COMPILER_VERSION)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX)
  IF (${PROJECT_NAME}_ENABLE_CXX)
    ENABLE_LANGUAGE(CXX)
    INCLUDE(CMakeDetermineCXXCompiler)
    PRINT_VAR(CMAKE_CXX_COMPILER_ID)
    PRINT_VAR(CMAKE_CXX_COMPILER_VERSION)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_Fortran)
    ENABLE_LANGUAGE(Fortran)
  ENDIF()

  # Do some project-specific tweaks for compiler options, etc.
  SET(PROJECT_COMPILER_CONFIG_FILE
    # Can be used for things like Kokkos.
    "${${PROJECT_NAME}_SOURCE_DIR}/cmake/ProjectCompilerPostConfig.cmake"
    CACHE FILEPATH
    "Allow for project-specific compiler settings."
   )
  IF (EXISTS "${PROJECT_COMPILER_CONFIG_FILE}")
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${PROJECT_COMPILER_CONFIG_FILE}")
    INCLUDE("${PROJECT_COMPILER_CONFIG_FILE}")
  ENDIF()

  # Set up for strong compiler warnings and warnings as errors

  INCLUDE(TribitsSetupBasicCompileLinkFlags)
  TRIBITS_SETUP_BASIC_COMPILE_LINK_FLAGS()

  #
  # The compilers are set, the environment is known to CMake.  Now set the
  # installation paths and options.
  #
  TRIBITS_SETUP_INSTALLATION_PATHS()

  # Set up Windows interface stuff

  IF (MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE
      -D_CRT_NONSTDC_NO_DEPRECATE  -D_SCL_SECURE_NO_WARNINGS)
    SET(WIN_INTERFACE_INCL  ${${PROJECT_NAME}_TRIBITS_DIR}/win_interface/include)
    IF (EXISTS "${WIN_INTERFACE_INCL}")
      INCLUDE_DIRECTORIES("${WIN_INTERFACE_INCL}")
      IF (TRIBITS_SETUP_ENV_DEBUG)
        MESSAGE("-- Adding win_interface/include ...")
      ENDIF()
    ENDIF()
  ENDIF()

  IF (WIN32 AND NOT CYGWIN)
    SET(NATIVE_MS_WINDOWS TRUE)
  ELSE()
    SET(NATIVE_MS_WINDOWS FALSE)
  ENDIF()

  # Probe for non-standard headers

  IF (${PROJECT_NAME}_ENABLE_CXX)
    CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
    CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)
    CHECK_INCLUDE_FILE_CXX(stdint.h HAVE_STDINT_H)
    CHECK_INCLUDE_FILE_CXX(inttypes.h HAVE_INTTYPES_H)
  ENDIF()

  SET(HAVE_ALGORITHM TRUE)
  SET(HAVE_CASSERT TRUE)
  SET(HAVE_CCTYPE TRUE)
  SET(HAVE_CERRNO TRUE)
  SET(HAVE_CLIMITS TRUE)
  SET(HAVE_CMATH TRUE)
  SET(HAVE_COMPLEX TRUE)
  SET(HAVE_CSTDARG TRUE)
  SET(HAVE_CSTDIO TRUE)
  SET(HAVE_CSTDLIB TRUE)
  SET(HAVE_CSTRING TRUE)
  SET(HAVE_IOMANIP TRUE)
  SET(HAVE_IOSTREAM TRUE)
  SET(HAVE_ITERATOR TRUE)
  SET(HAVE_LIST TRUE)
  SET(HAVE_MAP TRUE)
  SET(HAVE_MEMORY TRUE)
  SET(HAVE_MUTABLE TRUE)
  SET(HAVE_NAMESPACES TRUE)
  SET(HAVE_NEW_FOR_SCOPING TRUE)
  SET(HAVE_NUMERIC TRUE)
  SET(HAVE_NUMERIC_LIMITS TRUE)
  SET(HAVE_POW TRUE)
  SET(HAVE_SET TRUE)
  SET(HAVE_SSTREAM TRUE)
  SET(HAVE_FSTREAM TRUE)
  SET(HAVE_STDEXCEPT TRUE)
  SET(HAVE_STRING TRUE)
  SET(HAVE_VECTOR TRUE)

  # 2008/12/20: rabartl: Above: All of these defines should be removed
  # because we decided that we were going to assume that all compilers
  # have these C++98 standard features.  We will deal with cases where
  # this is not true but we should not assume the worst right from the
  # beginning.

  # Find Perl

  FIND_PACKAGE(Perl)

  # Do Fortran stuff

  INCLUDE(TribitsFortranMangling)

  # Get BLAS name mangling
  #
  # ToDo: Make this a project-specific specialization

  INCLUDE(TribitsBLASMangling)

  # Determine C++-0x supported features

  IF (${PROJECT_NAME}_ENABLE_CXX AND ${PROJECT_NAME}_ENABLE_CXX11)
    INCLUDE(TribitsCXX11Support)
    TRIBITS_FIND_CXX11_FLAGS() # Aborts if can't find C++11 flags!
    TRIBITS_CHECK_CXX11_SUPPORT(CXX11_WORKS)  # Double check that C++11 flags!
    IF (CXX11_WORKS)
      MESSAGE("-- ${PROJECT_NAME}_ENABLE_CXX11=${${PROJECT_NAME}_ENABLE_CXX11}")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${PROJECT_NAME}_CXX11_FLAGS}")
        IF (TRIBITS_SETUP_ENV_DEBUG OR TRIBITS_ENABLE_CXX11_DEBUG_DUMP)
          PRINT_VAR(CMAKE_CXX_FLAGS)
        ENDIF()
    ELSE()
      MESSAGE(FATAL_ERROR
        "Error, C++11 support does not appear to be supported"
        " with this C++ compiler and/or with the C++11 flags"
        " ${PROJECT_NAME}_CXX11_FLAGS='${${PROJECT_NAME}_CXX11_FLAGS}'!"
        " If the flags ${PROJECT_NAME}_CXX11_FLAGS='${${PROJECT_NAME}_CXX11_FLAGS}'"
        " where set manually, then try clearing the CMake cache and configure"
        " without setting "
        " ${PROJECT_NAME}_CXX11_FLAGS and let the configure process try to"
        " find flags that work automatically.  However, if these compile-time"
        " tests still fail, consider selecting a different C++ compiler"
        " (and compatible compilers for other languages) that supports C++11."
        " Or, if C++11 support in this project is not needed or desired, then set"
        " -D${PROJECT_NAME}_ENABLE_CXX11=OFF.")
    ENDIF()
  ENDIF()

  # Set up some MPI info

  IF (TPL_ENABLE_MPI)
    SET(HAVE_MPI TRUE)
  ELSE()
    SET(HAVE_MPI FALSE)
  ENDIF()

  # OpenMP isn't really a TPL because support is built into the compiler.
  IF(${PROJECT_NAME}_ENABLE_OpenMP)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
      TRIBITS_SET_OPENMP_FLAGS(CXX)
      TRIBITS_SET_OPENMP_FLAGS(C)
      IF(OpenMP_Fortran_FLAGS)
        TRIBITS_SET_OPENMP_FLAGS(Fortran)
      ELSE()
      # Older versions of FindOpenMP.cmake don't find Fortran flags.  Mike H said this is safe.
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
      ENDIF()
    ELSE()
      MESSAGE(FATAL_ERROR "Could not find OpenMP, try setting OpenMP_C_FLAGS and OpenMP_CXX_FLAGS directly")
    ENDIF(OPENMP_FOUND)
  ENDIF(${PROJECT_NAME}_ENABLE_OpenMP)

  # Check if we need the math library or not and find the right one
  IF (NOT NATIVE_MS_WINDOWS)
    INCLUDE(MathLibraryNeeded)
  ENDIF()

  # Check for isnan and isinf support
  IF (${PROJECT_NAME}_ENABLE_CXX)
    INCLUDE(FiniteValue)
  ENDIF()

  # Check for Doxygen/dot - We can use variables set in this check to
  # enable/disable the grapical dependency graphs in doxygen Doxyfiles.
  INCLUDE(FindDoxygen)

  # Set the hack library to get link options on

  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    IF (TRIBITS_SETUP_ENV_DEBUG)
      MESSAGE(STATUS "Creating dummy last_lib for appending the link flags: "
        "${${PROJECT_NAME}_EXTRA_LINK_FLAGS}")
    ENDIF()
    IF (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c)
      FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c
        "typedef int last_lib_dummy_t;\n")
    ENDIF()
    ADD_LIBRARY(last_lib STATIC ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c)
    TARGET_LINK_LIBRARIES(last_lib ${${PROJECT_NAME}_EXTRA_LINK_FLAGS})
  ENDIF()

  # You have to override the configuration types for MSVS after the compiler
  # checks!
  SET(CMAKE_CONFIGURATION_TYPES  ${CMAKE_CONFIGURATION_TYPE}
    CACHE STRING
    "Override by TriBITS (see TribitsDevelopersGuilde.*)"
    FORCE)
  IF (CMAKE_CONFIGURATION_TYPES)
    PRINT_VAR(CMAKE_CONFIGURATION_TYPES)
  ENDIF()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(SETUP_ENV_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${SETUP_ENV_TIME_START_SECONDS}
      ${SETUP_ENV_TIME_STOP_SECONDS}
      "\nTotal time to probe and setup the environment")
  ENDIF()

ENDMACRO()


MACRO(TRIBITS_SET_OPENMP_FLAGS  LANG)
  IF (NOT "${OpenMP_${LANG}_FLAGS_OVERRIDE}" STREQUAL "")
    SET(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} ${OpenMP_${LANG}_FLAGS_OVERRIDE}")
  ELSE()
    SET(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} ${OpenMP_${LANG}_FLAGS}")
  ENDIF()
ENDMACRO()


#
# Set mapping of labels to subprojects (i.e. TriBITS packages) for local CTest
# only.
#
# NOTE: This macro is only used define mapping of labels to subprojects for
# running ctest locally.  This results in summarizing the tests run for each
# subproject (TriBITS package) if any tests were run.  Therefore, it is
# harmless to define the mapping for every TriBITS package.  Only TriBITS
# packages will be listed in the summary if they had one or more tests run.
#

MACRO(TRIBITS_SET_LABELS_TO_SUBPROJECTS_MAPPING)
  SET(CTEST_LABELS_FOR_SUBPROJECTS ${${PROJECT_NAME}_PACKAGES})
ENDMACRO()


#
# Macro to turn on CTest support
#

MACRO(TRIBITS_INCLUDE_CTEST_SUPPORT)

  SET(DART_TESTING_TIMEOUT_IN ${DART_TESTING_TIMEOUT})

  IF (DART_TESTING_TIMEOUT_IN)
    TRIBITS_SCALE_TIMEOUT(${DART_TESTING_TIMEOUT} DART_TESTING_TIMEOUT)
    IF (NOT DART_TESTING_TIMEOUT STREQUAL DART_TESTING_TIMEOUT_IN)
     MESSAGE("-- DART_TESTING_TIMEOUT=${DART_TESTING_TIMEOUT_IN} being scaled by ${PROJECT_NAME}_SCALE_TEST_TIMEOUT=${${PROJECT_NAME}_SCALE_TEST_TIMEOUT} to ${DART_TESTING_TIMEOUT}")
    ENDIF()
    # Have to set DART_TESTING_TIMEOUT in cache or CMake will not put in right
    # 'TimeOut' in DartConfiguration.tcl file!
    SET(DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT} CACHE STRING "" FORCE)
  ENDIF()

  # Set up CTEst/CDash subprojects
  TRIBITS_SET_LABELS_TO_SUBPROJECTS_MAPPING()
  # NOTE: We do this after all of the packages have been defined but before
  # the DartConfiguration.tcl file has been created.

  INCLUDE(CTest)  # Generates file DartConfiguration.tcl with 'TimeOut' set!

  IF (DART_TESTING_TIMEOUT_IN)
    # Put DART_TESTING_TIMEOUT back to user input value to avoid scaling this
    # up and up on recofigures!
    SET(DART_TESTING_TIMEOUT ${DART_TESTING_TIMEOUT_IN} CACHE STRING
      "Original value set by user reset by TriBITS after scaling" FORCE)
  ENDIF()

  TRIBITS_CONFIGURE_CTEST_CUSTOM(${${PROJECT_NAME}_SOURCE_DIR}
    ${${PROJECT_NAME}_BINARY_DIR})

ENDMACRO()
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
# that with a SET( ... CACHE ...) statement in an input *.cmake file.


#
# Function that determines if a package should be processed
#

FUNCTION(TRIBITS_DETERMINE_IF_PROCESS_PACKAGE  PACKAGE_NAME
  PROCESS_PACKAGE_OUT  PACKAGE_ENABLE_STR_OUT
  )

  SET(PROCESS_PACKAGE FALSE)
  SET(PACKAGE_ENABLE_STR "")

  IF (${PACKAGE_NAME}_SUBPACKAGES)
    # Process the package if any of the the subpackages are enable
    FOREACH(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
      SET(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
      IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        SET(PROCESS_PACKAGE TRUE)
        APPEND_STRING_VAR_WITH_SEP(PACKAGE_ENABLE_STR ", " ${TRIBITS_SUBPACKAGE})
      ENDIF()
    ENDFOREACH()
  ELSE()
    # If the package itself is enabled, of course process it
    IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
      SET(PROCESS_PACKAGE TRUE)
      APPEND_STRING_VAR_WITH_SEP(PACKAGE_ENABLE_STR ", " "Libs")
    ENDIF()
  ENDIF()

  # If subpackages or package is enabled, then check tests/examples
  IF (PROCESS_PACKAGE)
    IF (${PACKAGE_NAME}_ENABLE_TESTS)
      APPEND_STRING_VAR_WITH_SEP(PACKAGE_ENABLE_STR ", " "Tests")
    ENDIF()
    IF (${PACKAGE_NAME}_ENABLE_EXAMPLES)
      APPEND_STRING_VAR_WITH_SEP(PACKAGE_ENABLE_STR ", " "Examples")
    ENDIF()
  ENDIF()

  SET(${PROCESS_PACKAGE_OUT} ${PROCESS_PACKAGE} PARENT_SCOPE)
  SET(${PACKAGE_ENABLE_STR_OUT} ${PACKAGE_ENABLE_STR} PARENT_SCOPE)

ENDFUNCTION()


#
# Macro that reads in the project's version file into the current scope
#

MACRO(TRIBITS_PROJECT_READ_VERSION_FILE  PROJECT_SOURCE_DIR_IN)
  SET(PROJECT_VERSION_FILE ${PROJECT_SOURCE_DIR_IN}/Version.cmake)
  IF (EXISTS ${PROJECT_VERSION_FILE})
    # Set REPOSITORY_NAME in case Version.cmake is written generically!
    SET(REPOSITORY_NAME ${PROJECT_NAME})
    TRIBITS_TRACE_FILE_PROCESSING(PROJECT  INCLUDE  "${PROJECT_VERSION_FILE}")
    INCLUDE(${PROJECT_VERSION_FILE})
  ENDIF()
ENDMACRO()


#
# Function that reads in and the Repository's specific Version.cmake file and
# then configures its ${REPO_NAME}_version.h file.
#
# The file ${REPO_NAME}_version.h is only configured if the repository contains
# the files Version.cmake and Copyright.txt
#
# NOTE: This is done as a function so that the read-in version variables don't
# bleed into the outer scope.
#
FUNCTION(TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE
  REPOSITORY_NAME  REPOSITORY_DIR  ADD_INSTALL_TARGET
  OUTPUT_VERSION_HEADER_FILE
  )

  IF (TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE_DEBUG_DUMP)
    MESSAGE("TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE: "
      "'${REPOSITORY_NAME}'  '${REPOSITORY_DIR}"
      "  '${OUTPUT_VERSION_HEADER_FILE}'")
  ENDIF()

  STRING(TOUPPER ${REPOSITORY_NAME} REPOSITORY_NAME_UC)

  TRIBITS_SET_BASE_REPO_DIR(${PROJECT_SOURCE_DIR} ${REPOSITORY_DIR}
    REPOSITORY_ABS_DIR)

  SET(REPOSITORY_VERSION_FILE ${REPOSITORY_ABS_DIR}/Version.cmake)
  SET(REPOSITORY_COPYRIGHT_FILE ${REPOSITORY_ABS_DIR}/Copyright.txt)

  IF (EXISTS ${REPOSITORY_VERSION_FILE} AND EXISTS ${REPOSITORY_COPYRIGHT_FILE})

    # Read the copyright header info
    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  READ  "${REPOSITORY_COPYRIGHT_FILE}")
    FILE(READ "${REPOSITORY_COPYRIGHT_FILE}" REPOSITORY_COPYRIGHT_HEADER)

    # Read the version variables and translate into standard form
    TRIBITS_TRACE_FILE_PROCESSING(REPOSITORY  INCLUDE  "${REPOSITORY_VERSION_FILE}")
    INCLUDE(${REPOSITORY_VERSION_FILE})
    SET(REPOSITORY_MAJOR_VERSION ${${REPOSITORY_NAME}_MAJOR_VERSION})
    SET(REPOSITORY_MAJOR_MINOR_VERSION ${${REPOSITORY_NAME}_MAJOR_MINOR_VERSION})
    SET(REPOSITORY_VERSION_STRING ${${REPOSITORY_NAME}_VERSION_STRING})

    # Configure the file with everything set
    IF (TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE_DEBUG_DUMP)
      MESSAGE("-- Writing the file ${OUTPUT_VERSION_HEADER_FILE} ...")
    ENDIF()
    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_PACKAGE_ARCH_DIR}/Tribits_version.h.in
      ${OUTPUT_VERSION_HEADER_FILE})

    IF (ADD_INSTALL_TARGET)
      # Install version header file
      TRIBITS_INSTALL_HEADERS(HEADERS  ${OUTPUT_VERSION_HEADER_FILE})
    ENDIF()

  ENDIF()

ENDFUNCTION()


#
# Configure each of the Repositories version header files
#

FUNCTION(TRIBITS_REPOSITORY_CONFIGURE_ALL_VERSION_HEADER_FILES)
  #PRINT_VAR(ARGN)
  FOREACH(REPO ${ARGN})
    TRIBITS_GET_REPO_NAME_DIR(${REPO}  REPO_NAME  REPO_DIR)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Considering configuring version file for '${REPO_NAME}'")
    ENDIF()
    TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE( ${REPO_NAME}  ${REPO_DIR}  TRUE
      "${${PROJECT_NAME}_BINARY_DIR}/${REPO_DIR}/${REPO_NAME}_version.h")
  ENDFOREACH()
ENDFUNCTION()


#
# Macro that does the final set of package configurations
#

MACRO(TRIBITS_CONFIGURE_ENABLED_PACKAGES)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(CONFIGURE_PACKAGES_TIME_START_SECONDS)
  ENDIF()

  #
  # A) Global variable initialization
  #

  GLOBAL_NULL_SET(${PROJECT_NAME}_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PROJECT_NAME}_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PROJECT_NAME}_LIBRARIES)
  GLOBAL_NULL_SET(${PROJECT_NAME}_ETI_PACKAGES)

  #
  # B) Define the source and binary directories for all of the pacakges that
  # have been enbaled.  These are used to allow packages to refer to each
  # other even downstream packages (which is pretty messed up really).
  #

  SET(PACKAGE_IDX 0)
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})

    # Get all the package sources independent of whether they are enabled or not.
    # There are some messed up packages that grab parts out of unrelated
    # downstream packages that might not even be enabled.  To support this,
    # allow this.
    LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    #PRINT_VAR(${TRIBITS_PACKAGE}_SOURCE_DIR)

    TRIBITS_DETERMINE_IF_PROCESS_PACKAGE(${TRIBITS_PACKAGE}
       PROCESS_PACKAGE  PACKAGE_ENABLE_STR)

    IF (PROCESS_PACKAGE)

      IF (${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR)
        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR)
        ENDIF()
        IF(IS_ABSOLUTE ${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
          SET(${TRIBITS_PACKAGE}_BINARY_DIR ${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
        ELSE()
          SET(${TRIBITS_PACKAGE}_BINARY_DIR
            ${CMAKE_CURRENT_BINARY_DIR}/${${TRIBITS_PACKAGE}_SPECIFIED_BINARY_DIR})
        ENDIF()
      ELSE()
        SET(${TRIBITS_PACKAGE}_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_DIR})
      ENDIF()
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${TRIBITS_PACKAGE}_BINARY_DIR)
      ENDIF()

    ENDIF()

    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  ENDFOREACH()

  #
  # C) Loop over all of the packages and process their CMakeLists.txt files if
  # they are enabled or if any of their subpackages are enabled.
  #

  # Include these here so they don't need to be included in each package's
  # CMakeLists.txt files.
  INCLUDE(TribitsPackageMacros)
  INCLUDE(TribitsSubPackageMacros)
  INCLUDE(AddSubdirectories)

  SET(CONFIGURED_A_PACKAGE FALSE)
  SET(ENABLED_PACKAGE_LIBS_TARGETS)

  # Tell packages that are also repos they are being processed as a package.
  SET(TRIBITS_PROCESSING_PACKAGE TRUE)

  SET(PACKAGE_IDX 0)
  FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGES})

    TRIBITS_DETERMINE_IF_PROCESS_PACKAGE(${TRIBITS_PACKAGE}
      PROCESS_PACKAGE  PACKAGE_ENABLE_STR)

    IF (PROCESS_PACKAGE)

      MESSAGE("Processing enabled package: ${TRIBITS_PACKAGE} (${PACKAGE_ENABLE_STR})")

      IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

        IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
            AND
            ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
              OR ${TRIBITS_PACKAGE}_PACKAGE_CONFIGURE_TIMING )
          )
          TIMER_GET_RAW_SECONDS(PROCESS_THIS_PACKAGE_TIME_START_SECONDS)
        ENDIF()

        SET(PACKAGE_NAME ${TRIBITS_PACKAGE}) # Used in CMake code in downstream package
        SET(PARENT_PACKAGE_NAME ${TRIBITS_PACKAGE})
        STRING(TOUPPER "${PARENT_PACKAGE_NAME}" PARENT_PACKAGE_NAME_UC)

        IF (NOT EXISTS ${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt)
          MESSAGE(FATAL_ERROR
            "Error, the file ${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt does not exist!")
        ENDIF()

        IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          PRINT_VAR(${TRIBITS_PACKAGE}_SOURCE_DIR)
          PRINT_VAR(${TRIBITS_PACKAGE}_BINARY_DIR)
        ENDIF()

        SET(TRIBITS_PACKAGE_CMAKELIST_FILE
          "${${TRIBITS_PACKAGE}_SOURCE_DIR}/CMakeLists.txt")
        TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  ADD_SUBDIR
          "${TRIBITS_PACKAGE_CMAKELIST_FILE}")
        IF (NOT ${TRIBITS_PACKAGE}_SOURCE_DIR STREQUAL ${PROJECT_NAME}_SOURCE_DIR)
          ADD_SUBDIRECTORY(${${TRIBITS_PACKAGE}_SOURCE_DIR} ${${TRIBITS_PACKAGE}_BINARY_DIR})
	ELSE()
          INCLUDE("${TRIBITS_PACKAGE_CMAKELIST_FILE}")
        ENDIF()
        IF (NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS)
          MESSAGE(FATAL_ERROR
            "ERROR: Forgot to call TRIBITS_PACKAGE_POSTPROCESS() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}"
            )
        ENDIF()

        LIST(APPEND ENABLED_PACKAGE_LIBS_TARGETS ${TRIBITS_PACKAGE}_libs)
        LIST(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
        LIST(APPEND ${PROJECT_NAME}_LIBRARY_DIRS ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
        LIST(APPEND ${PROJECT_NAME}_LIBRARIES ${${TRIBITS_PACKAGE}_LIBRARIES})

        IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
            AND
            ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
              OR ${TRIBITS_PACKAGE}_PACKAGE_CONFIGURE_TIMING )
          )
          TIMER_GET_RAW_SECONDS(PROCESS_THIS_PACKAGE_TIME_STOP_SECONDS)
          TIMER_PRINT_REL_TIME(${PROCESS_THIS_PACKAGE_TIME_START_SECONDS}
            ${PROCESS_THIS_PACKAGE_TIME_STOP_SECONDS}
            "-- Total time to configure package ${TRIBITS_PACKAGE}")
        ENDIF()

      ENDIF()

      SET(CONFIGURED_A_PACKAGE TRUE)

    ENDIF()

    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  ENDFOREACH()


  #
  # D) Loop backwards over ETI packages if ETI is enabled
  #

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    # Do this regardless of whether project level ETI is enabled
    IF("${${PROJECT_NAME}_ETI_PACKAGES}" STREQUAL "")
      MESSAGE("\nNo ETI support requested by packages.\n")
    ELSE()
      #IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE("\nProcessing explicit instantiation support for enabled packages ...\n")
      #ENDIF()
      SET(REVERSE_ETI_LIST ${${PROJECT_NAME}_ETI_PACKAGES})
      LIST(REVERSE REVERSE_ETI_LIST)
      FOREACH(PACKAGE_NAME ${REVERSE_ETI_LIST})
        MESSAGE("Processing ETI support: ${PACKAGE_NAME}")
        IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
            AND
            ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
              OR ${PACKAGE_NAME}_PACKAGE_CONFIGURE_TIMING )
          )
          TIMER_GET_RAW_SECONDS(PROCESS_ETI_START_SECONDS)
        ENDIF()
        SET(ETIFILE ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/ExplicitInstantiationSupport.cmake)
        IF(NOT EXISTS "${ETIFILE}")
          MESSAGE(FATAL_ERROR "Could not find ${PACKAGE_NAME} ETI support file ${ETIFILE}")
        ENDIF()
        TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  INCLUDE  "${ETIFILE}")
        INCLUDE("${ETIFILE}")
        IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING
            AND
            ( ${PROJECT_NAME}_ENABLE_PACKAGE_CONFIGURE_TIMING
              OR ${PACKAGE_NAME}_PACKAGE_CONFIGURE_TIMING )
          )
          TIMER_GET_RAW_SECONDS(PROCESS_ETI_STOP_SECONDS)
          TIMER_PRINT_REL_TIME(${PROCESS_ETI_START_SECONDS}
            ${PROCESS_ETI_STOP_SECONDS}
            "-- Time to process ETI support for package ${PACKAGE_NAME}")
        ENDIF()
      ENDFOREACH()
    ENDIF()

  ENDIF()

  #
  # E) Check if no packages are enabled and if that is allowed
  #

  ADVANCED_SET( ${PROJECT_NAME}_ALLOW_NO_PACKAGES ON
    CACHE BOOL "Allow configuration to finish even if no packages are enabled")

  IF (NOT CONFIGURED_A_PACKAGE)
    IF (${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      SET(MSG_TYPE WARNING)
    ELSE()
      SET(MSG_TYPE ERROR)
    ENDIF()
    MESSAGE(
      "\n***"
      "\n*** ${MSG_TYPE}:  There were no packages configured so no libraries"
        " or tests/examples will be built!"
      "\n***\n"
      )
    IF (NOT ${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      MESSAGE(SEND_ERROR "Stopping configure!")
    ENDIF()
  ELSE()
    ASSERT_AND_TOUCH_DEFINED(${PROJECT_NAME}_ALLOW_NO_PACKAGES)
  ENDIF()

  #
  # F) Process the global variables and other cleanup
  #

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_INCLUDE_DIRS)
    REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_LIBRARY_DIRS)
    REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_LIBRARIES)

    # Add global 'libs' target
    IF(ENABLED_PACKAGE_LIBS_TARGETS)
      LIST(REVERSE ENABLED_PACKAGE_LIBS_TARGETS)
      # Make it so when no packages are enabled it is not a cmake error
      IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
        APPEND_SET(ENABLED_PACKAGE_LIBS_TARGETS last_lib)
      ENDIF()
      #PRINT_VAR(ENABLED_PACKAGE_LIBS_TARGETS)
      IF (NOT TARGET ${PROJECT_NAME}_libs)
        ADD_CUSTOM_TARGET(${PROJECT_NAME}_libs)
        ADD_DEPENDENCIES(${PROJECT_NAME}_libs ${ENABLED_PACKAGE_LIBS_TARGETS})
      ENDIF()
      ADD_CUSTOM_TARGET(libs)
      ADD_DEPENDENCIES(libs ${ENABLED_PACKAGE_LIBS_TARGETS})
    ENDIF()

  ENDIF()

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(CONFIGURE_PACKAGES_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${CONFIGURE_PACKAGES_TIME_START_SECONDS}
      ${CONFIGURE_PACKAGES_TIME_STOP_SECONDS}
      "\nTotal time to configure enabled packages")
  ENDIF()

ENDMACRO()


#
# Set up for packaging and distribution
#

MACRO(TRIBITS_SETUP_PACKAGING_AND_DISTRIBUTION)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    # Start the global timer
    TIMER_GET_RAW_SECONDS(CPACK_SETUP_TIME_START_SECONDS)
  ENDIF()

  # K.1) Run callback function for the base project.

  TRIBITS_PROJECT_DEFINE_PACKAGING_RUNNER()
  # The above must define the basic project settings for CPACK that are
  # specific to the project and should not be provided by the user.

  # K.2) Removing any packages or SE packages not enabled from the tarball

  IF (${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION)
    SET(_SE_OR_FULL_PACKAGES ${${PROJECT_NAME}_SE_PACKAGES})
    SET(_SE_OR_FULL_PACKAGE_DIRS ${${PROJECT_NAME}_SE_PACKAGE_DIRS})
  ELSE()
    SET(_SE_OR_FULL_PACKAGES ${${PROJECT_NAME}_PACKAGES})
    SET(_SE_OR_FULL_PACKAGE_DIRS ${${PROJECT_NAME}_PACKAGE_DIRS})
  ENDIF()

  TRIBITS_GET_NONENABLED_LIST(
    _SE_OR_FULL_PACKAGES  ${PROJECT_NAME}
    NON_ENABLED_SE_OR_FULL_PACKAGES  NUM_NON_ENABLED_SE_OR_FULL_PACKAGES)
  #PRINT_VAR(NON_ENABLED_SE_OR_FULL_PACKAGES)

  FOREACH(TRIBITS_PACKAGE ${NON_ENABLED_SE_OR_FULL_PACKAGES})

    # Determine if this is a package to not ignore
    FIND_LIST_ELEMENT(TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE
       ${TRIBITS_PACKAGE}  TRIBITS_PACKAGE_DONT_IGNORE)

    IF (NOT TRIBITS_PACKAGE_DONT_IGNORE)

      LIST(FIND _SE_OR_FULL_PACKAGES ${TRIBITS_PACKAGE} PACKAGE_IDX)
      LIST(GET _SE_OR_FULL_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
      # ToDo: Repalce the above O(N) LIST(FIND ...) with a O(1) lookup ...

      # Checking if we have a relative path to the package's files. Since the
      # exclude is a regular expression any "../" will be interpretted as <any
      # char><any char>/ which would never match the package's actual
      # directory. There isn't a direct way in cmake to convert a relative
      # path into an absolute path with string operations so as a way of
      # making sure that we get the correct path of the package we use a
      # find_path for the CMakeLists.txt file for the package. Since the
      # package has to have this file to work correctly it should be
      # guaranteed to be there.
      STRING(REGEX MATCH "[.][.]/" IS_RELATIVE_PATH ${PACKAGE_DIR})
      IF("${IS_RELATIVE_PATH}" STREQUAL "")
        SET(CPACK_SOURCE_IGNORE_FILES "${PROJECT_SOURCE_DIR}/${PACKAGE_DIR}/"
          ${CPACK_SOURCE_IGNORE_FILES})
      ELSE()
        FIND_PATH(ABSOLUTE_PATH  CMakeLists.txt  PATHS
          ${PROJECT_SOURCE_DIR}/${PACKAGE_DIR} NO_DEFAULT_PATH)
        IF("${ABSOLUTE_PATH}" STREQUAL "ABSOLUTE_PATH-NOTFOUND")
          MESSAGE(AUTHOR_WARNING "Relative path found for disabled package"
            " ${TRIBITS_PACKAGE} but package was missing a CMakeLists.txt file."
            " This disabled package will likely not be excluded from a source release")
        ENDIF()
        SET(CPACK_SOURCE_IGNORE_FILES ${ABSOLUTE_PATH} ${CPACK_SOURCE_IGNORE_FILES})
      ENDIF()
    ENDIF()

  ENDFOREACH()

  # Add excludes for VC files/dirs
  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /[.]git/
    [.]gitignore$
    )

  # Print the set of excluded files
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE OR
    ${PROJECT_NAME}_DUMP_CPACK_SOURCE_IGNORE_FILES
    )
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()

  # K.3) Set up install component dependencies

  TRIBITS_GET_ENABLED_LIST(
    ${PROJECT_NAME}_PACKAGES  ${PROJECT_NAME}
    ENABLED_PACKAGES  NUM_ENABLED)
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")

  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  # K.4) Resetting the name to avoid overwriting registery keys when installing

  IF(WIN32)
    SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-${${PROJECT_NAME}_VERSION}")
    IF (TPL_ENABLE_MPI)
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-mpi")
    ELSE ()
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-serial")
    ENDIF()
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()

  # K.5) Determine the source generator
  IF ("${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT}" STREQUAL "")
    SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ")
  ENDIF()
  SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR
    ${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT}
    CACHE STRING
    "The types of source generators to use for CPACK_SOURCE_GENERATOR.")
  SET(CPACK_SOURCE_GENERATOR ${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR})

  # K.6) Loop through the Repositories and run their callback functions.
  FOREACH(REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    TRIBITS_GET_REPO_NAME_DIR(${REPO}  REPO_NAME  REPO_DIR)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("Processing packaging call-backs for ${REPO_NAME}")
    ENDIF()
    TRIBITS_REPOSITORY_DEFINE_PACKAGING_RUNNER(${REPO_NAME})
  ENDFOREACH()

  # K.7) Include <Project>RepoVersion.txt if generated
  SET(PROJECT_REPO_VERSION_FILE
     "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}")
  IF (EXISTS "${PROJECT_REPO_VERSION_FILE}")
    FOREACH(SOURCE_GEN ${CPACK_SOURCE_GENERATOR})
      SET(CPACK_INSTALL_COMMANDS ${CPACK_INSTALL_COMMANDS}
        "${CMAKE_COMMAND} -E copy '${PROJECT_REPO_VERSION_FILE}' '${CMAKE_CURRENT_BINARY_DIR}/_CPack_Packages/Linux-Source/${SOURCE_GEN}/${CPACK_PACKAGE_NAME}-${${PROJECT_NAME}_VERSION}-Source/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}'")
    ENDFOREACH()
  ENDIF()

  # K.8) Finally process with CPack
  INCLUDE(CPack)

  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(CPACK_SETUP_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${CPACK_SETUP_TIME_START_SECONDS}  ${CPACK_SETUP_TIME_STOP_SECONDS}
      "Total time to set up for CPack packaging")
  ENDIF()

ENDMACRO()


#
# Setup for installation
#

MACRO(TRIBITS_SETUP_FOR_INSTALLATION)

  # Set up to install <Package>Config.cmake, <Project>Config.cmake, and export
  # makefiles.

  IF((${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
      OR ${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    AND NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    )

    INCLUDE(TribitsWriteClientExportFiles)

    TRIBITS_WRITE_PROJECT_CLIENT_EXPORT_FILES()

    IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
      # TEMPORARY: Install a compatibility copy of ${PROJECT_NAME}Config.cmake
      # where was previously installed to warn and load the new file.
      SET(COMPATIBILITY_CONFIG_INCLUDE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
      CONFIGURE_FILE(
        ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsConfigInclude.cmake.in
        ${COMPATIBILITY_CONFIG_INCLUDE}
        @ONLY
        )
      INSTALL(
        FILES ${COMPATIBILITY_CONFIG_INCLUDE}
        DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
        )
    ENDIF()

  ENDIF()

ENDMACRO()


#
# @MACRO: TRIBITS_EXCLUDE_FILES()
#
# Exclude package files/dirs from the source distribution by appending
# ``CPACK_SOURCE_IGNORE_FILES``.
#
# Usage::
#
#  TRIBITS_EXCLUDE_FILES(<file0> <file1> ...)
#
# This is called in the package's top-level `<packageDir>/CMakeLists.txt`_
# file and each file or directory name ``<filei>`` is actually interpreted by
# CMake/CPack as a regex that is prefixed by the project's and package's
# source directory names so as to not exclude files and directories of the
# same name and path from other packages.  If ``<filei>`` is an absolute path
# it is not prefixed but is appended to ``CPACK_SOURCE_IGNORE_FILES``
# unmodified.
#
# In general, do **NOT** put in excludes for files and directories that are
# not under this package's source tree.  If the given package is not enabled,
# then this command will never be called! For example, don't put in excludes
# for PackageB's files in PackageA's ``CMakeLists.txt`` file because if
# PackageB is enabled but PackageA is not, the excludes for PackageB will
# never get added to ``CPACK_SOURCE_IGNORE_FILES``.
#
# Also, be careful to note that the ``<filei>`` arguments are actually regexes
# and one must be very careful not understand how CPack will use these regexes
# to match files that get excluded from the tarball.  For more details, see
# `Creating Source Distributions`_.
#
MACRO(TRIBITS_EXCLUDE_FILES)

  SET(FILES_TO_EXCLUDE ${ARGN})

  # Need to add "/<project source dir>/<package dir>/" to each file to prevent
  # someone from trying to exclude a file like "readme" and having it
  # inadvertently exclude a file matching that name in another package.
  SET(MODIFIED_FILES_TO_EXCLUDE "")

  SET(${PROJECT_NAME}_SOURCE_PATH ${${PROJECT_NAME}_SOURCE_DIR})

  LIST(FIND ${PROJECT_NAME}_PACKAGES ${PACKAGE_NAME} PACKAGE_IDX)
  LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)

  FOREACH(FILE ${FILES_TO_EXCLUDE})
    #Ensure that if the full path was specified for the file that we don't add
    #"/<project source dir>/<package dir>/" again.
    SET(MATCH_STRING "${${PROJECT_NAME}_SOURCE_PATH}/${PACKAGE_DIR}")
    STRING(REGEX MATCH ${MATCH_STRING} MATCHED ${FILE} )
    IF(NOT MATCHED)
      LIST(APPEND MODIFIED_FILES_TO_EXCLUDE
        ${${PROJECT_NAME}_SOURCE_PATH}/${PACKAGE_DIR}/${FILE})
    ELSE()
      LIST(APPEND MODIFIED_FILES_TO_EXCLUDE ${FILE})
    ENDIF()
  ENDFOREACH()

#Leaving in for debugging purposes
#  MESSAGE("List of files being excluded for package ${PACKAGE_NAME}")
#  FOREACH(NEW_FILE ${MODIFIED_FILES_TO_EXCLUDE})
#    MESSAGE(${NEW_FILE})
#  ENDFOREACH()

  LIST(APPEND CPACK_SOURCE_IGNORE_FILES ${MODIFIED_FILES_TO_EXCLUDE})
  IF (NOT ${PROJECT_NAME}_BINARY_DIR STREQUAL ${PACKAGE_NAME}_BINARY_DIR)
    SET(CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES} PARENT_SCOPE)
  ENDIF()

ENDMACRO()


#
#  Macro for helping set up exclude files only for the packages that will not
#  be supporting autotools.
#
MACRO(TRIBITS_EXCLUDE_AUTOTOOLS_FILES) # PACKAGE_NAME LIST_RETURN)
  SET(AUTOTOOLS_FILES
    configure.ac$
    configure$
    Makefile.am$
    Makefile.in$
    bootstrap$
    .*[.]m4$
    config/
    )

  SET(FILES_TO_EXCLUDE)
  FOREACH(FILE ${AUTOTOOLS_FILES})
    LIST(APPEND FILES_TO_EXCLUDE ${FILE} \(.*/\)*${FILE})
  ENDFOREACH()

  TRIBITS_EXCLUDE_FILES(${FILES_TO_EXCLUDE})

ENDMACRO()
