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

INCLUDE(TribitsPackageSetupCompilerFlags)
INCLUDE(TribitsWriteClientExportFiles)
INCLUDE(TribitsGeneralMacros)

INCLUDE(ParseVariableArguments)
INCLUDE(GlobalNullSet)
INCLUDE(AppendGlobalSet)
INCLUDE(PrintVar)
INCLUDE(PrependSet)
INCLUDE(PrependGlobalSet)
INCLUDE(RemoveGlobalDuplicates)

INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddAdvancedTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of file!
###


#
# Utility macros
#


#
# Macro that defines the package architecture system varaibles used to link
# different SE packages together
#
# See README.DEPENDENCIES for information on what these varaibles mean and how
# they are used.
#
MACRO(TRIBITS_DEFINE_LINKAGE_VARS PACKAGE_NAME_IN)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME_IN}_LIBRARIES)
  GLOBAL_SET(${PACKAGE_NAME_IN}_HAS_NATIVE_LIBRARIES FALSE)
ENDMACRO()


#
# Macro that defines varaibles that create global targets
#
MACRO(TRIBITS_DEFINE_TARGET_VARS PARENT_PACKAGE_NAME_IN)
  GLOBAL_NULL_SET(${PARENT_PACKAGE_NAME_IN}_LIB_TARGETS)
  GLOBAL_NULL_SET(${PARENT_PACKAGE_NAME_IN}_ALL_TARGETS)
ENDMACRO()

#
# Set up some common varaibles used in the creation of an SE package
#

MACRO(TRIBITS_SET_COMMON_VARS PACKAGE_NAME_IN)

  STRING(TOUPPER ${PACKAGE_NAME_IN} PACKAGE_NAME_UC)

  # Write TRIBITS_PACKAGE versions of common variables
  SET(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this ${PROJECT_NAME} package is in
  FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  STRING(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

ENDMACRO()


#
# @MACRO: TRIBITS_PACKAGE_DECL()
#
# Macro called at the very beginning of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file when a package has subpackages.
#
# Usage::
#
#   TRIBITS_PACKAGE_DECL(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# The arguments are:
#
#   ``<packageName>``
#
#     Gives the name of the Package, mostly just for checking and
#     documentation purposes.  This must match the name of the package
#     provided in the `<repoDir>/PackagesList.cmake`_ or an error is issued.
#
#   ``ENABLE_SHADOWING_WARNINGS``
#
#     If specified, then shadowing warnings for the package's sources will be
#     turned on for supported platforms/compilers.  The default is for
#     shadowing warnings to be turned off.  Note that this can be overridden
#     globally by setting the cache variable
#     ``${PROJECT_NAME}_ENABLE_SHADOWING_WARNINGS``.
#
#   ``DISABLE_STRONG_WARNINGS``
#
#     If specified, then all strong warnings for the package's sources will be
#     turned off, if they are not already turned off by global cache
#     variables.  Strong warnings are turned on by default in development
#     mode.
#
#   ``CLEANED``
#
#     If specified, then warnings will be promoted to errors for compiling the
#     package's sources for all defined warnings.
#
#   ``DISABLE_CIRCULAR_REF_DETECTION_FAILURE``
#
#     If specified, then the standard grep looking for RCPNode circular
#     references in `TRIBITS_ADD_TEST()`_ and `TRIBITS_ADD_ADVANCED_TEST()`_
#     that causes tests to fail will be disabled.  Note that if these warnings
#     are being produced then it means that the test is leaking memory and
#     user like may also be leaking memory.
#
# There are several side-effects of calling this macro:
#
# * The variables ``${PACKAGE_NAME}_LIB_TARGETS`` (lists all of the package's
#   targets) and ``${PACKAGE_NAME}_ALL_TARGETS`` (lists all of the package's
#   libraries) and are initialized to empty.
#
# * The local variables ``PACKAGE_SOURCE_DIR`` and ``PACKAGE_BINARY_DIR`` are
#   set for this package's use in its CMakeLists.txt files.
#
# * Package-specific compiler options are set up in package-scope (i.e., the
#   package's subdirs) in ``CMAKE_<LANG>_FLAG``.
#
# * This packages's cmake subdir ``${PACKAGE_SOURCE_DIR}/cmake`` is added to
#   ``CMAKE_MODULE_PATH`` locally so that the package's try-compile modules
#   can be read in with just a raw ``INCLUDE()`` leaving off the full path and
#   the ``*.cmake`` extension.
#
# If the package does not have subpackages, just call `TRIBITS_PACKAGE()`_
# which calls this macro.
#
MACRO(TRIBITS_PACKAGE_DECL PACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_PACKAGE_DECL: ${PACKAGE_NAME_IN}")
  ENDIF()

  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    ""
    #options
    "CLEANED;ENABLE_SHADOWING_WARNINGS;DISABLE_STRONG_WARNINGS;DISABLE_CIRCULAR_REF_DETECTION_FAILURE"
    ${ARGN}
    )

  #
  # B) Assert that the global and local package names are the same!
  #

  IF (DEFINED PACKAGE_NAME)
    IF (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME})
      MESSAGE(FATAL_ERROR "Error, the package-defined package name"
        " '${PACKAGE_NAME_IN}' is not the same as the package name"
        " defined at the global level '${PACKAGE_NAME}'")
    ENDIF()
  ENDIF()

  #
  # C) Set up the CMake support for this ${PROJECT_NAME} package and define some
  # top-level varaibles.
  #

  TRIBITS_SET_COMMON_VARS(${PACKAGE_NAME_IN})

  SET(${PACKAGE_NAME_IN}_DISABLE_STRONG_WARNINGS OFF
     CACHE BOOL
     "If set to true, then strong warnings for package ${PACKAGE_NAME_IN} will be disabled."
     )

  # Set up the compile flags for the package
  TRIBITS_SETUP_COMPILER_FLAGS(${PACKAGE_NAME_IN})

  # Set up circular reference detection test failure
  IF (PARSE_DISABLE_CIRCULAR_REF_DETECTION_FAILURE)
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)
  ELSE()
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE ON)
  ENDIF()

  # Set up parent package linkage varaibles
  TRIBITS_DEFINE_TARGET_VARS(${PACKAGE_NAME})

  #
  # Append the local package's cmake directory in order to help pull in
  # configure-time testing macros
  #

  PREPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

ENDMACRO()


#
# @MACRO: TRIBITS_PACKAGE_DEF()
#
# Macro called in `<packageDir>/CMakeLists.txt`_ after subpackages are
# processed in order to handle the libraries, tests, and examples of the
# parent package.
#
# Usage::
#
#   TRIBITS_PACKAGE_DEF()
#
# If the package does not have subpackages, just call `TRIBITS_PACKAGE()`_
# which calls this macro.
#
# This macro has several side effects:
#
# * The variable ``PACKAGE_NAME`` is set in the local scope for usage by the
#   package's ``CMakeLists.txt`` files.
#
# * The intra-package dependency variables (i.e. list of include directories,
#   list of libraries, etc.) are initialized to empty.
#
MACRO(TRIBITS_PACKAGE_DEF)

  # Reset since it was changed by the subpackages
  SET(PACKAGE_NAME ${PARENT_PACKAGE_NAME})

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_PACKAGE_DEF: ${PACKAGE_NAME}")
  ENDIF()

  IF (NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("\n${PACKAGE_NAME} not enabled so exiting package processing")
    ENDIF()
    RETURN()
  ENDIF()

  # Reset in case were changed by subpackages
  TRIBITS_SET_COMMON_VARS(${PACKAGE_NAME})

  # Define package linkage varaibles
  TRIBITS_DEFINE_LINKAGE_VARS(${PACKAGE_NAME})

ENDMACRO()


#
# @MACRO: TRIBITS_PACKAGE()
#
# Macro called at the very beginning of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   TRIBITS_PACKAGE(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# See `TRIBITS_PACKAGE_DECL()`_ for the documentation for the arguments and
# `TRIBITS_PACKAGE_DECL()`_ and `TRIBITS_PACKAGE()`_ for a description the
# side-effects (and variables set) after calling this macro.
#
MACRO(TRIBITS_PACKAGE PACKAGE_NAME_IN)
  TRIBITS_PACKAGE_DECL(${PACKAGE_NAME_IN} ${ARGN})
  TRIBITS_PACKAGE_DEF()
ENDMACRO()


#
# @MACRO: TRIBITS_ADD_TEST_DIRECTORIES()
#
# Macro called to add a set of test directories for an SE package.
#
# Usage::
#
#    TRIBITS_ADD_TEST_DIRECTORIES(<dir1> <dir2> ...)
#
# This macro only needs to be called from the top most ``CMakeLists.txt`` file
# for which all subdirectories are all "tests".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# Currently, all this macro does macro is to call ``ADD_SUBDIRECTORY(<diri>)``
# if ``${PACKAGE_NAME}_ENABLE_TESTS`` or
# ``${PARENT_PACKAGE_NAME}_ENABLE_TESTS`` are ``TRUE``. However, this macro
# may be extended in the future in order to modify behavior related to adding
# tests and examples in a uniform way.
#
MACRO(TRIBITS_ADD_TEST_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_TESTS OR ${PARENT_PACKAGE_NAME}_ENABLE_TESTS)
    FOREACH(TEST_DIR ${ARGN})
      TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  ADD_SUBDIR
        "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/CMakeLists.txt")
      ADD_SUBDIRECTORY(${TEST_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Macros to add common options to add to an SE package
#


#
# @MACRO: TRIBITS_ADD_DEBUG_OPTION()
#
# Add the standard cache variable option ``${PACKAGE_NAME}_ENABLE_DEBUG`` for
# the package.
#
# Usage::
#
#   TRIBITS_ADD_DEBUG_OPTION()
#
# This option is given the default ``${${PROJECT_NAME}_ENABLE_DEBUG}`` and if
# true, will set the variable ``HAVE_${PACKAGE_NAME_UC}_DEBUG`` (to be used in
# the package's configured header file).  This macro is typically called in
# the package's `<packageDir>/CMakeLists.txt`_ file.
#
MACRO(TRIBITS_ADD_DEBUG_OPTION)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_DEBUG
    HAVE_${PACKAGE_NAME_UC}_DEBUG
    "Enable a host of runtime debug checking."
    ${${PROJECT_NAME}_ENABLE_DEBUG}
    )
ENDMACRO()


MACRO(TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION)
  OPTION(
    ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
     "Enable Teuchos time monitors for package ${PACKAGE_NAME}"
    ${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}
    )
ENDMACRO()


#
# @MACRO: TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()
#
# Add the standard option ``${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS`` for the
# package.
#
# Usage::
#
#   TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()
#
# This macro should be called in the package's <packageDir>/CMakeLists.txt`_
# file.  This option is given the default value
# ``${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}``.  This option is then looked
# for in `TRIBITS_CONFIGURE_FILE()`_ to add macros to add deprecated warnings
# to deprecated parts of a package.
#
MACRO(TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION)
  ADVANCED_SET(
    ${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS  ${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}
    CACHE BOOL
    "Show warnings about deprecated code in ${PACKAGE_NAME}"
    )
  ADVANCED_SET(
    ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE  ${${PROJECT_NAME}_HIDE_DEPRECATED_CODE}
    CACHE BOOL
    "Fully exclude deprecated code in ${PACKAGE_NAME}"
    )
ENDMACRO()


MACRO(TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
    HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
    "Enable the use of explicit template instantiation."
    ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
    )
ENDMACRO()


MACRO(TRIBITS_ADD_ETI_SUPPORT)
  APPEND_GLOBAL_SET(${PROJECT_NAME}_ETI_PACKAGES ${PACKAGE_NAME})
  GLOBAL_NULL_SET(${PACKAGE_NAME}_ETI_LIBRARYSET)
ENDMACRO()


#
# @MACRO: TRIBITS_ADD_EXAMPLE_DIRECTORIES()
#
# Macro called to conditionally add a set of example directories for an SE
# package.
#
# Usage::
#
#    TRIBITS_ADD_EXAMPLE_DIRECTORIES(<dir1> <dir2> ...)
#
# This macro typically is called from the top-level
# `<packageDir>/CMakeLists.txt`_ file for which all subdirectories are all
# "examples" according to standard package layout.
#
# This macro can be called several times within a package as desired to break
# up example directories any way one would like.
#
# Currently, all it does macro does is to call ``ADD_SUBDIRECTORY(<diri>)`` if
# ``${PACKAGE_NAME}_ENABLE_EXAMPLES`` or
# ``${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES`` are true. However, this macro may
# be extended in the future in order to modify behavior related to adding
# tests and examples in a uniform way.
#
MACRO(TRIBITS_ADD_EXAMPLE_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES)
    FOREACH(EXAMPLE_DIR ${ARGN})
      TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  ADD_SUBDIR
        "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_DIR}/CMakeLists.txt")
      ADD_SUBDIRECTORY(${EXAMPLE_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Utility function that sets up package linkage linkage variables in case the
# package has no libraries.
#

FUNCTION(TRIBITS_PACKAGE_FINALIZE_DEPENDENCY_VARS)

  IF(${PACKAGE_NAME}_SUBPACKAGES)

    # A package with subpackages should get all of its dependency vars from
    # its enabled subpackages.

    SET(PARENT_PACKAGE_INCLUDE_DIRS)
    SET(PARENT_PACKAGE_LIBRARY_DIRS)
    SET(PARENT_PACKAGE_LIBRARIES)

    SET(SUBPACKAGE_IDX 0)
    FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        PREPEND_SET(PARENT_PACKAGE_INCLUDE_DIRS
          ${${SUBPACKAGE_FULLNAME}_INCLUDE_DIRS})
        PREPEND_SET(PARENT_PACKAGE_LIBRARY_DIRS
          ${${SUBPACKAGE_FULLNAME}_LIBRARY_DIRS})
        PREPEND_SET(PARENT_PACKAGE_LIBRARIES
          ${${SUBPACKAGE_FULLNAME}_LIBRARIES})
      ENDIF()

      MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    ENDFOREACH()

    IF (PARENT_PACKAGE_INCLUDE_DIRS)
      LIST(REMOVE_DUPLICATES PARENT_PACKAGE_INCLUDE_DIRS)
    ENDIF()
    IF (PARENT_PACKAGE_LIBRARY_DIRS)
      LIST(REMOVE_DUPLICATES PARENT_PACKAGE_LIBRARY_DIRS)
    ENDIF()
    # NOTE: Above, in the rare case that none of the subpackages contain any
    # libraries or any include directories, we need to not call
    # LIST(REMOVE_DUPLICATES ...).

    # NOTE: There can't be any dupicate libraries in PARENT_PACKAGE_LIBRARIES
    # so no need to remove them.

    GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS "${PARENT_PACKAGE_INCLUDE_DIRS}")
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS "${PARENT_PACKAGE_LIBRARY_DIRS}")
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES "${PARENT_PACKAGE_LIBRARIES}")

  ELSEIF(NOT ${PACKAGE_NAME}_INCLUDE_DIRS)

    # No libraries have been defined for this package so we are going to set
    # them based on this package's dependencies.

    TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)

    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${LINK_LIBS})

  ENDIF()

ENDFUNCTION()


#
# Helper macro for [SUB]TRIBITS_PACKAGE_POSTPROCESS()
#
MACRO(TRIBITS_PACKAGE_POSTPROCESS_COMMON)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_PACKAGE_POSTPROCESS_COMMON: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES OR
    ${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES
    )
    # Create the configure file so external projects can find packages with a
    # call to find_package(<package_name>).  This also creates the
    # Makefile.export.* files.
    TRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES(${PACKAGE_NAME})
  ENDIF()

  SET(${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE TRUE
    CACHE INTERNAL "")

ENDMACRO()


#
# @MACRO: TRIBITS_PACKAGE_POSTPROCESS()
#
# Macro called at the very end of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file that performs some critical
# post-processing activities.
#
# Usage::
#
#   TRIBITS_PACKAGE_POSTPROCESS()
#
# NOTE: It is unfortunate that this macro must be called in a packages's
# top-level ``CMakeLists.txt`` file but limitations of the CMake language make
# it necessary to do so.
#
MACRO(TRIBITS_PACKAGE_POSTPROCESS)

  # Only parent packages have the targets (${PACKAGE_NAME}_libs and
  # (${PACKAGE_NAME}_all
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_PACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_LIB_TARGETS)
    PRINT_VAR(${PACKAGE_NAME}_ALL_TARGETS)
  ENDIF()
  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_libs DEPENDS ${${PACKAGE_NAME}_LIB_TARGETS})
  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_all DEPENDS ${${PACKAGE_NAME}_ALL_TARGETS})

  TRIBITS_PACKAGE_FINALIZE_DEPENDENCY_VARS()
  TRIBITS_PACKAGE_POSTPROCESS_COMMON()

  # NOTE: This package is only able

ENDMACRO()


#
# @MACRO: TRIBITS_PROCESS_SUBPACKAGES()
#
# Macro that processes the `TriBITS Subpackages`_ for a parent `TriBITS
# package`_ for packages that are broken down into subpackages.  This is
# called in the parent packages top-level `<packageDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   TRIBITS_PROCESS_SUBPACKAGES()
#
# This macro must be called after `TRIBITS_PACKAGE_DECL()`_ but before
# `TRIBITS_PACKAGE_DEF()`_.
#
MACRO(TRIBITS_PROCESS_SUBPACKAGES)

  #MESSAGE("TRIBITS_PROCESS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  SET(SUBPACKAGE_IDX 0)
  FOREACH(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

    #MESSAGE("")
    #PRINT_VAR(SUBPACKAGE_IDX)
    #PRINT_VAR(TRIBITS_SUBPACKAGE)

    SET(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
    SET(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #PRINT_VAR(SUBPACKAGE_FULLNAME)

    IF (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})

      LIST(GET ${PARENT_PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
      #PRINT_VAR(SUBPACKAGE_DIR)

      IF (NOT ${PROJECT_NAME}_BINARY_DIR STREQUAL ${PARENT_PACKAGE_NAME}_BINARY_DIR)
        DUAL_SCOPE_SET(${SUBPACKAGE_FULLNAME}_BINARY_DIR 
          ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/${SUBPACKAGE_DIR})
      ELSE()
        SET(${SUBPACKAGE_FULLNAME}_BINARY_DIR 
          ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/${SUBPACKAGE_DIR})
      ENDIF()

      TRIBITS_TRACE_FILE_PROCESSING(PACKAGE  ADD_SUBDIR
        "${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}/CMakeLists.txt")
      ADD_SUBDIRECTORY(${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}
        ${${SUBPACKAGE_FULLNAME}_BINARY_DIR})

    ENDIF()

    MATH(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

  ENDFOREACH()

ENDMACRO()


##################################################################
#
#                    NOTES TO DEVELOPERS
#
# Don't even attempt to touch the logic that goes into setting up and
# modifying the variables:
#
#   ${PACKAGE_NAME}_INCLUDE_DIRS
#   ${PACKAGE_NAME}_LIBRARY_DIRS
#   ${PACKAGE_NAME}_LIBRARIES
#   ${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES
#   ${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES
#   ${PARENT_PACKAGE_NAME}_LIB_TARGETS
#   ${PARENT_PACKAGE_NAME}_ALL_TARGETS
#
# without carefully studying the documentation in README.DEPENENCIES and then
# carefully studying all of the code and issues that modify these variables!
#
# ToDo: Write some good unit tests that pin down the behavior of all of this!
#
