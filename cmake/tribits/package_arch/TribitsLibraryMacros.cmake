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

INCLUDE(TribitsCreateClientTemplateHeaders)
INCLUDE(ParseVariableArguments)
INCLUDE(GlobalSet)
INCLUDE(AppendSet)
INCLUDE(AppendGlob)
INCLUDE(AppendGlobalSet)
INCLUDE(AppendStringVar)
INCLUDE(PrependGlobalSet)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(TribitsGeneralMacros)
INCLUDE(SetAndIncDirs)


#
# Macro that configures the package's main config.h file
#

FUNCTION(TRIBITS_ADD_CONFIG_DEFINE DEFINE)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPackage ${PARENT_PACKAGE_NAME}: adding compiler"
      " define to config file: ${DEFINE}")
  ENDIF()
  GLOBAL_SET(${PARENT_PACKAGE_NAME}_CONFIG_DEFINES
    "${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}\n#define ${DEFINE}")
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("--${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}")
  ENDIF()
ENDFUNCTION()


#
# Macro that configures the package's main config.h file
#

FUNCTION(TRIBITS_CONFIGURE_FILE PACKAGE_NAME_CONFIG_FILE)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  ENDIF()

  # Get an upper case version of the parent package name
  STRING(TOUPPER "${PARENT_PACKAGE_NAME}" PARENT_PACKAGE_NAME_UC)

  # Set up the deprecated attribute if showing deprecated warnings
  IF (${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS)
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#ifndef ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  endif\n"
      "#endif\n"
      )
  ELSE()
    SET(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED")
  ENDIF()

  IF (${PARENT_PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    APPEND_STRING_VAR(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "\n#define ${PARENT_PACKAGE_NAME_UC}_HIDE_DEPRECATED_CODE")
  ENDIF()

  # Set up the macro to create the define for time monitor
  SET(TIME_MONITOR_DEFINE_NAME ${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR)
  SET(FUNC_TIME_MONITOR_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR)
  SET(FUNC_TIME_MONITOR_DIFF_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR_DIFF)
  IF (${PARENT_PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR)
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#ifndef ${FUNC_TIME_MONITOR_MACRO_NAME}\n"
      "#  define ${TIME_MONITOR_DEFINE_NAME}\n"
      "#  define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME) \\\\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, ${PARENT_PACKAGE_NAME_UC})\n"
      "#  define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF) \\\\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)\n"
      "#endif\n"
      )
  ELSE()
    MULTILINE_SET(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME)\n"
      "#define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF)\n"
      )
  ENDIF()

  # Configure the file
  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

ENDFUNCTION()


#
# Macro used to add a package library
#
# ToDo: Add documentation!
#

FUNCTION(TRIBITS_ADD_LIBRARY LIBRARY_NAME)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_LIBRARY: ${LIBRARY_NAME}")
    IF(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
      MESSAGE("\n${PACKAGE_NAME}_LIBRARIES In installation testing mode,"
        " libraries will be found instead of created.")
    ENDIF()
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  PARSE_ARGUMENTS(
    PARSE #prefix
    "HEADERS;NOINSTALLHEADERS;SOURCES;DEPLIBS;IMPORTEDLIBS;DEFINES" # Lists
    "TESTONLY;NO_INSTALL_LIB_OR_HEADERS;CUDALIBRARY" #Options
    ${ARGN} # Remaining arguments passed in
    )

  #if we are doing installation testing we want to skip adding libraries, unless
  #they are test only libraries which are not installed.
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_TESTONLY)

    # Add the link directory for this library.

    SET_PROPERTY(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${CMAKE_CURRENT_BINARY_DIR})

    # NOTE: Above, this link path not really used here for anything.
    # Instead it is just added to the other set link library directories
    # that are already set.  These link directories are then extracted
    # and stored into stored in ${PACKAGE_NAME}_LIBRARY_DIRS.

    # Add whatever include directories have been defined so far

    INCLUDE_DIRECTORIES(AFTER ${${PACKAGE_NAME}_INCLUDE_DIRS})

    # Add whatever link directories have been added so far

    SET_PROPERTY(DIRECTORY  APPEND  PROPERTY  PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_LIBRARY_DIRS})

    # Local varaible to hold all of the libraries that will be directly linked
    # to this library.
    SET(LINK_LIBS)

    # Add dependent libraries passed directly in

    IF (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "DEPLIBS = ${PARSE_DEPLIBS}")
    ENDIF()
    IF (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    ENDIF()

    IF (PARSE_DEPLIBS)
      APPEND_SET(LINK_LIBS ${PARSE_DEPLIBS})
    ENDIF()
    IF (PARSE_IMPORTEDLIBS)
      APPEND_SET(LINK_LIBS ${PARSE_IMPORTEDLIBS})
    ENDIF()

    #
    # We only want to link to the dependent package and TPL libraries when we need
    # to.  We only need to link to these dependent libraries when this is the first
    # library being created for this package or if this library does not depend
    # on other libraries created for this package.  Otherwise, we don't need to
    # add the include directories or link libraries because a dependent lib
    # specified in PARSE_DEP_LIBS already has everything that we need.
    #
    # We also need to make special considerations for test libraries since
    # things need to be handled a little bit differently (but not much).  In the
    # case of test libaries, we need to also pull the test-only dependencies.
    # In this case, we will always assume that we will add in the test
    # libraries.
    #

    SET(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

    IF (PARSE_DEPLIBS AND NOT PARSE_TESTONLY)
      FOREACH(DEPLIB ${PARSE_DEPLIBS})
        LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${DEPLIB} DEPLIB_IDX)
        IF (NOT DEPLIB_IDX EQUAL -1)
          # The library being created here is dependent on another of this
          # package's libraries so there is no need to add in this package's
          # dependent package and TPL libraries.
          SET(ADD_DEP_PACKAGE_AND_TPL_LIBS FALSE)
        ENDIF()
      ENDFOREACH()
    ELSE()
      # If there are no dependent libs passed in, then this library can not
      # possiblly depend on the package's other libraries so we must link to
      # the dependent libraries in dependent libraries and TPLs.
    ENDIF()

    IF (ADD_DEP_PACKAGE_AND_TPL_LIBS)

      IF (NOT PARSE_TESTONLY)
        SET(LIB_OR_TEST_ARG LIB)
      ELSE()
        SET(LIB_OR_TEST_ARG TEST)
      ENDIF()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "\nPulling in header and libraries dependencies"
          " for ${LIB_OR_TEST_ARG} ...")
      ENDIF()

      #
      # Call INCLUDE_DIRECTORIES() and LINK_DIRECTORIES(...) for all upstream
      # dependent Packages and TPLs and accumulate the libraries to link against.
      #
      # NOTE: Adding these directories serves two purposes.  First, so that the includes
      # get added the the sources that get built for this library.  Second, so
      # that list full list of include directories can be extracted as a
      # propery and set on ${PACKAGE_NAME}_INCLUDE_DIRS
      #

      TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}  LINK_LIBS)

      TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
        ${PACKAGE_NAME}  ${LIB_OR_TEST_ARG}  LINK_LIBS)

    ENDIF()

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(LINK_LIBS)
    ENDIF()

    # Add the library and all the dependencies

    IF (PARSE_DEFINES)
      ADD_DEFINITIONS(${PARSE_DEFINES})
    ENDIF()

    IF (NOT PARSE_CUDALIBRARY)
      ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ELSE()
      CUDA_ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ENDIF()

    SET_PROPERTY(TARGET ${LIBRARY_NAME} APPEND PROPERTY
      LABELS ${PACKAGE_NAME})

    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
    PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

    TARGET_LINK_LIBRARIES(${LIBRARY_NAME}  ${LINK_LIBS})

    # Add to the install target

    SET(INSTALL_LIB ON)
    SET(INSTALL_HEADERS ON)
    SET(APPEND_LIB_AND_HEADERS_TO_PACKAGE ON)

    IF (PARSE_TESTONLY)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping installation hooks for this library"
          " because 'TESTONLY' was passed in ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
      SET(APPEND_LIB_AND_HEADERS_TO_PACKAGE OFF)
    ELSEIF (PARSE_NO_INSTALL_LIB_OR_HEADERS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping installation hooks for this library"
          " because 'NO_INSTALL_LIB_OR_HEADERS' was passed in ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
    ELSEIF (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND NOT BUILD_SHARED_LIBS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping installation of headers and libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=FALSE  ...")
      ENDIF()
      SET(INSTALL_LIB OFF)
      SET(INSTALL_HEADERS OFF)
    ELSEIF (NOT ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS AND BUILD_SHARED_LIBS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping installation of headers but installing libraries"
          " because ${PROJECT_NAME}_INSTALL_LIBRARIES_AND_HEADERS=FALSE and BUILD_SHARED_LIBS=TRUE  ...")
      ENDIF()
      SET(INSTALL_HEADERS OFF)
    ENDIF()

    IF (INSTALL_LIB OR INSTALL_HEADERS)
      SET_PROPERTY(GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS ON)
      SET_PROPERTY(GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS ON)
    ENDIF()

    IF (INSTALL_LIB)
      SET_PROPERTY(TARGET ${LIBRARY_NAME} PROPERTY INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")
      INSTALL(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PROJECT_NAME}
          RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
          LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
          ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
          COMPONENT ${PACKAGE_NAME}
        )
    ENDIF()

    IF (INSTALL_HEADERS)
      INSTALL(
        FILES ${PARSE_HEADERS}
        DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )
    ENDIF()

    # Append the new include dirs, library dirs, and libraries to this package's lists

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    IF (APPEND_LIB_AND_HEADERS_TO_PACKAGE)

      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
      IF (PARSE_IMPORTEDLIBS)
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${PARSE_IMPORTEDLIBS})
      ENDIF()
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${LIBRARY_NAME})

      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_INCLUDE_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARY_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARIES)

    ELSE()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping augmentation of package's lists of include"
          " directories and libraries! ...")
      ENDIF()

      GLOBAL_SET(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${LIBRARY_NAME}_INCLUDE_DIRS)
      ENDIF()

    ENDIF()
  ENDIF() #if not in installation testing mode

  IF (${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)

    LIST(FIND ${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST ${PACKAGE_NAME}
      ${PACKAGE_NAME}_WAS_INSTALLED)
    IF(${${PACKAGE_NAME}_WAS_INSTALLED} EQUAL -1)
      MESSAGE(FATAL_ERROR
        "The package ${PACKAGE_NAME} was not installed with ${PROJECT_NAME}!"
        "  Please disable package ${PACKAGE_NAME} or install it.")
    ENDIF()

    INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING  BEFORE
       ${${PACKAGE_NAME}_INSTALLATION_INCLUDE_DIRS}
       ${${TRIBITS_PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS})
    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
      ${${PACKAGE_NAME}_INSTALLATION_LIBRARY_DIRS})

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT PACKAGE_LIBRARY_DIRS)

    GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES    ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})

  ENDIF() #instalation testing mode

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

ENDFUNCTION()
