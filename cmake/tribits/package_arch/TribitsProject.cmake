# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER
 

# Projects that change the location of the source need to consider this in
# their top-level CMakeLists.txt file
SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits"
  CACHE PATH
  "The base directory pointing to the TriBITS system."
  )
MARK_AS_ADVANCED(${PROJECT_NAME}_TRIBITS_DIR)

SET(CMAKE_MODULE_PATH
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/cmake
   ${${PROJECT_NAME}_TRIBITS_DIR}/utils
   ${${PROJECT_NAME}_TRIBITS_DIR}/package_arch
   ${${PROJECT_NAME}_TRIBITS_DIR}/config_tests
   ${${PROJECT_NAME}_TRIBITS_DIR}/modules
   ${${PROJECT_NAME}_TRIBITS_DIR}/installation
   )

IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("CMAKE_MODULE_PATH='${CMAKE_MODULE_PATH}'")
ENDIF() 
 
# Overrides that we have for CMake functions
INCLUDE(CMakeOverrides)

INCLUDE(TribitsConstants)
INCLUDE(TribitsGlobalMacros)
INCLUDE(AdvancedSet)
INCLUDE(AdvancedOption)
INCLUDE(TimingUtils)
INCLUDE(SetDefault)


#
# Defines a TriBITS project.
#
# Requires that PROJECT_NAME be defined before calling this macro.
#
# ToDo: Give documentation
#

MACRO(TRIBITS_PROJECT)

  #
  # A) Basic top-level TriBITS project stuff
  #
  
  MESSAGE("")
  MESSAGE("Configuring ${PROJECT_NAME} build directory")
  MESSAGE("")
  
  IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    MESSAGE(FATAL_ERROR "ERROR! "
      "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
      " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
      "\n${PROJECT_NAME} does not support in source builds!\n"
      "NOTE: You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under"
      " the source directory for ${PROJECT_NAME} or you will not be able to configure ${PROJECT_NAME} correctly!"
      "\nYou must now run something like:\n"
      "  $ rm -r CMakeCache.txt CMakeFiles/"
      "\n"
      "Please create a different directory and configure ${PROJECT_NAME} under that such as:\n"
      "  $ mkdir MY_BUILD\n"
      "  $ cd MY_BUILD\n"
      "  $ cmake [OPTIONS] .."
      )
  ENDIF()
  
  STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)
  SET(PROJECT_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
  SET(PROJECT_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
  PRINT_VAR(PROJECT_HOME_DIR)
  PRINT_VAR(PROJECT_BUILD_DIR)
  # Above, we put these in the cache so we can grep them out

  SET(${PROJECT_NAME}_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
  SET(${PROJECT_NAME}_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
  PRINT_VAR(${PROJECT_NAME}_HOME_DIR)
  PRINT_VAR(${PROJECT_NAME}_BUILD_DIR)
  
  MESSAGE("-- " "CMAKE_VERSION = ${CMAKE_VERSION}")
  
  IF (NOT ${PROJECT_NAME}_DEPS_HOME_DIR)
    SET(${PROJECT_NAME}_DEPS_HOME_DIR "${PROJECT_HOME_DIR}")
  ENDIF()

  SET(TRIBITS_PROJECT_CALLBACK_FILE
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${${PROJECT_NAME}_TRIBITS_PROJECT_CALLBACKS_FILE})
  IF(EXISTS ${TRIBITS_PROJECT_CALLBACK_FILE})
    INCLUDE(${TRIBITS_PROJECT_CALLBACK_FILE})
  ENDIF()

  
  TRIBITS_READ_IN_OPTIONS_FROM_FILE()

  #
  # A.1) Run misc unit tests that don't need anything else
  #
  # These below tests are a bit strange.  See
  # cmake/TestingUnitTests/CMakeLists.txt for details!
  #
  
  IF (${PROJECT_NAME}_INVOKE_TESTING_UNIT_TESTS)
    ADD_SUBDIRECTORY(cmake/TestingUnitTests/UnitTests)
    RETURN()
  ENDIF()
  
  #
  # A.2) Set up other stuff
  #
  
  INCLUDE(TribitsFindPythonInterp)
  TRIBITS_FIND_PYTHON()
  PRINT_VAR(PYTHON_EXECUTABLE)
  
  #
  # A.3) Set up version file that also sets other options as well
  #
  # NOTE: ${PROJECT_NAME}Version.cmake must be read *before* the global options are
  # read!
  #

  TRIBITS_CONFIGURE_VERSION_FILE(
    "${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}_version.h")
  
  # Since the version header file is now configured the root build
  # dir needs to be on the include path
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  
  #
  # B) Set up user options and global variables that will be used throughout
  #
  
  MESSAGE("")
  MESSAGE("Setting up major user options ...")
  MESSAGE("")
  
  TRIBITS_DEFINE_GLOBAL_OPTIONS()

  # Call-back defined by specific project
  TRIBITS_PROJECT_SETUP_EXTRA_OPTIONS()
  
  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    # Start the global timer
    TIMER_GET_RAW_SECONDS(GLOBAL_TIME_START_SECONDS)
  ENDIF()
  
  ADVANCED_OPTION(${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING
    "Shortcircut after dependency handling is complete"
    OFF )
  
  ADVANCED_OPTION(${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST
    "Skip the Fortran/C++ compatibility test"
    OFF )
  
  # Find an installed version of ${PROJECT_NAME} for installation testing
  # (the check that we are in installation mode is inside the macro)
  INCLUDE(TribitsInstallationTestingMacros)
  FIND_PROJECT_INSTALL()
  
  #
  # C) Read in ${PROJECT_NAME} packages and TPLs and process dependencies
  #
  
  TRIBITS_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()
  
  #
  # D) Apply logic to enable ${PROJECT_NAME} packages and tests
  #
  
  TRIBITS_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES()
  
  #
  # E) Stop if asked
  #
  
  IF (${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
    MESSAGE("")
    MESSAGE(SEND_ERROR "Shortcircuiting after dependency tracking ...")
    RETURN()
  ENDIF()
  
  # ToDo: rabartl: Remove the above once the unit tests have been refactored to
  # just run macros and not the entire system.
  
  #
  # F) Set up the environment on this computer
  #
  
  MESSAGE("")
  MESSAGE("Probing the environment ...")
  MESSAGE("")
  
  TRIBITS_SETUP_ENV()
  
  #
  # G) Go get the information for all enabled TPLS
  #
  
  MESSAGE("")
  MESSAGE("Getting information for all enabled TPLs ...")
  MESSAGE("")
  
  TRIBITS_PROCESS_ENABLED_TPLS()
  
  # OpenMP is similar to a TPL in some respects, but requires only compiler
  # flags to enable
  
  OPTION(${PROJECT_NAME}_ENABLE_OpenMP
    "Build with OpenMP support." OFF)
  
  #
  # H) Set up for testing with CTest and ${PROJECT_NAME} test harness
  #
  
  MESSAGE("")
  MESSAGE("Setting up testing support ...")
  MESSAGE("")
  
  INCLUDE(CTest)

  TRIBITS_PROJECT_SETUP_TESTING_SUPPORT()
  
  #
  # I) Add the 'dashboard' target
  #
  # NOTE: Must come after setting up for testing
  #
  
  TRIBITS_ADD_DASHBOARD_TARGET()
  
  
  #
  # J) Configure individual packages
  # 
  
  MESSAGE("")
  MESSAGE("Configuring individual enabled ${PROJECT_NAME} packages ...")
  MESSAGE("")
  
  TRIBITS_CONFIGURE_ENABLED_PACKAGES()
  
  
  #
  # K) Setup for packaging and distribution
  #
  
  TRIBITS_PROJECT_DEFINE_PACKAGING()
  
  
  #
  # L) Install-related commands
  #
  IF(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
    AND NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    )
  
    INCLUDE(TribitsPackageWritePackageConfig)
  
    TRIBITS_WRITE_CONFIG_FILE()
  
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
  
  
  #
  # M) Export the library dependencies. This will let client projects
  # refer to all TPLs used by ${PROJECT_NAME}. (KRL, 26 Nov 2009)
  #
  
  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    MESSAGE("")
    MESSAGE("Exporting library dependencies ...")
    MESSAGE("")
    EXPORT_LIBRARY_DEPENDENCIES( ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}LibraryDepends.cmake )
  ENDIF()
  
  #
  # P) Show final timing and end
  #
  
  MESSAGE("")
  MESSAGE("Finished configuring ${PROJECT_NAME}!")
  MESSAGE("")
  IF (${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING)
    TIMER_GET_RAW_SECONDS(GLOBAL_TIME_STOP_SECONDS)
    TIMER_PRINT_REL_TIME(${GLOBAL_TIME_START_SECONDS}  ${GLOBAL_TIME_STOP_SECONDS}
      "Total time to configure ${PROJECT_NAME}")
  ENDIF()
  
  TRIBITS_REMIND_ABOUT_UNCOMMITTED_DEPENDENCY_FILES()
  
ENDMACRO()
