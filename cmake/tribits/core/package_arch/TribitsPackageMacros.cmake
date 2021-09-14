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

include(TribitsPackageSetupCompilerFlags)
include(TribitsWriteClientExportFiles)
include(TribitsGeneralMacros)

include(CMakeParseArguments)
include(GlobalNullSet)
include(AppendGlobalSet)
include(PrintVar)
include(PrependSet)
include(PrependGlobalSet)
include(RemoveGlobalDuplicates)

include(TribitsAddOptionAndDefine)
include(TribitsLibraryMacros)
include(TribitsAddExecutable)
include(TribitsAddExecutableAndTest)
include(TribitsAddTest)
include(TribitsAddAdvancedTest)
include(TribitsCopyFilesToBinaryDir)
include(TribitsReportInvalidTribitsUsage)


###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of file!
###


#
# Utility macros
#


#
# Macro that defines the package architecture system variables used to link
# different SE packages together
#
# See README.DEPENDENCIES for information on what these variables mean and how
# they are used.
#
macro(tribits_define_linkage_vars PACKAGE_NAME_IN)
  global_null_set(${PACKAGE_NAME_IN}_INCLUDE_DIRS)
  global_null_set(${PACKAGE_NAME_IN}_LIBRARY_DIRS)
  global_null_set(${PACKAGE_NAME_IN}_LIBRARIES)
  global_set(${PACKAGE_NAME_IN}_HAS_NATIVE_LIBRARIES_TO_INSTALL FALSE)
endmacro()


#
# Macro that defines variables that create global targets
#
macro(tribits_define_target_vars PARENT_PACKAGE_NAME_IN)
  global_null_set(${PARENT_PACKAGE_NAME_IN}_LIB_TARGETS)
  global_null_set(${PARENT_PACKAGE_NAME_IN}_ALL_TARGETS)
endmacro()

#
# Set up some common variables used in the creation of an SE package
#

macro(tribits_set_common_vars PACKAGE_NAME_IN)

  string(TOUPPER ${PACKAGE_NAME_IN} PACKAGE_NAME_UC)

  # Write TRIBITS_PACKAGE versions of common variables
  set(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  set(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this ${PROJECT_NAME} package is in
  file(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  string(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

endmacro()


#
# @MACRO: tribits_package_decl()
#
# Macro called at the very beginning of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file when a package has subpackages.
#
# Usage::
#
#   tribits_package_decl(
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
#     references in `tribits_add_test()`_ and `tribits_add_advanced_test()`_
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
#   can be read in with just a raw ``include()`` leaving off the full path and
#   the ``*.cmake`` extension.
#
# If the package does not have subpackages, just call `tribits_package()`_
# which calls this macro.
#
macro(tribits_package_decl PACKAGE_NAME_IN)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_PACKAGE_DECL: ${PACKAGE_NAME_IN}")
  endif()

  if (CURRENTLY_PROCESSING_SUBPACKAGE)
    tribits_report_invalid_tribits_usage(
      "Cannot call tribits_package_decl() in a subpackage."
      " Use tribits_subpackage() instead"
      " error in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
  endif()

  if(${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
    tribits_report_invalid_tribits_usage(
      "tribits_package_decl() called more than once in Package ${PACKAGE_NAME}"
      " This may be because tribits_package_decl() was explicitly called more than once or"
      " TRIBITS_PACKAGE_DECL was called after TRIBITS_PACKAGE. You do not need both."
      " If your package has subpackages then do not call tribits_package() instead call:"
      " tribits_pacakge_decl() then tribits_process_subpackages() then tribits package_def()"
    )
  endif()

  # Set flag to check that macros are called in the correct order
  set(${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED TRUE)

  #
  # A) Parse the input arguments
  #

  cmake_parse_arguments(
    #prefix
    PARSE
    #options
    "CLEANED;ENABLE_SHADOWING_WARNINGS;DISABLE_STRONG_WARNINGS;DISABLE_CIRCULAR_REF_DETECTION_FAILURE"
    #one_value_keywords
    ""
    #multi_value_keywords
    ""
    ${ARGN}
    )

  tribits_check_for_unparsed_arguments()

  #
  # B) Assert that the global and local package names are the same!
  #

  if (DEFINED PACKAGE_NAME)
    if (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME})
      message(FATAL_ERROR "Error, the package-defined package name"
        " '${PACKAGE_NAME_IN}' is not the same as the package name"
        " defined at the global level '${PACKAGE_NAME}'")
    endif()
  endif()

  #
  # C) Set up the CMake support for this ${PROJECT_NAME} package and define some
  # top-level variables.
  #

  tribits_set_common_vars(${PACKAGE_NAME_IN})

  set(${PACKAGE_NAME_IN}_DISABLE_STRONG_WARNINGS OFF
     CACHE BOOL
     "If set to true, then strong warnings for package ${PACKAGE_NAME_IN} will be disabled."
     )

  # Set up the compile flags for the package
  tribits_setup_compiler_flags(${PACKAGE_NAME_IN})

  # Set up circular reference detection test failure
  if (PARSE_DISABLE_CIRCULAR_REF_DETECTION_FAILURE)
    set(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)
  else()
    set(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE ON)
  endif()

  # Set up parent package linkage variables
  tribits_define_target_vars(${PACKAGE_NAME})

  # Define this as a CMake/CTest "Subproject"
  set_directory_properties(PROPERTIES LABELS "${PACKAGE_NAME}")

  #
  # Append the local package's cmake directory in order to help pull in
  # configure-time testing macros
  #

  prepend_set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

endmacro()


#
# @MACRO: tribits_package_def()
#
# Macro called in `<packageDir>/CMakeLists.txt`_ after subpackages are
# processed in order to handle the libraries, tests, and examples of the
# parent package.
#
# Usage::
#
#   tribits_package_def()
#
# If the package does not have subpackages, just call `tribits_package()`_
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
macro(tribits_package_def)

  # check that this is not being called from a subpackage
  if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
    if (CURRENTLY_PROCESSING_SUBPACKAGE)
      tribits_report_invalid_tribits_usage(
        "Cannot call tribits_package_def() in a subpackage."
        " Use tribits_subpackage() instead"
        " error in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()
  endif()

  # Reset since it was changed by the subpackages
  set(PACKAGE_NAME ${PARENT_PACKAGE_NAME})

  # check that this is not called morethan once in a package
  if (${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED)
    tribits_report_invalid_tribits_usage(
      "tribits_package_def() was called more than once in"
      "${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_PACKAGE_DEF: ${PACKAGE_NAME}")
  endif()

  if (NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("\n${PACKAGE_NAME} not enabled so exiting package processing")
    endif()
    return()
  endif()

  # Reset in case were changed by subpackages
  tribits_set_common_vars(${PACKAGE_NAME})

  # Define package linkage variables
  tribits_define_linkage_vars(${PACKAGE_NAME})

  set(${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED TRUE)

endmacro()


#
# @MACRO: tribits_package()
#
# Macro called at the very beginning of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   tribits_package(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# See `tribits_package_decl()`_ for the documentation for the arguments and
# `tribits_package_decl()`_ and `tribits_package()`_ for a description the
# side-effects (and variables set) after calling this macro.
#
macro(tribits_package PACKAGE_NAME_IN)

  if (CURRENTLY_PROCESSING_SUBPACKAGE)
    if (NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Cannot call tribits_package() in a subpackage."
        " Use tribits_subpackage() instead"
        " error in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()
  endif()

  if(${PACKAGE_NAME}_SUBPACKAGES)
    tribits_report_invalid_tribits_usage(
      "This package has subpackages so you cannot use tribits_package()"
      "\n Instead use the following call order:"
      "\n tribits_project_decl(${PACKAGE_NAME})"
      "\n tribits_process_subpackages()"
      "\n [do other things you want to do]"
      "\n tribits_package_def()"
      "\n tribits_package_postprocess()" )
  endif()

  if(${PACKAGE_NAME}_TRIBITS_PACKAGE_CALLED)
    tribits_report_invalid_tribits_usage(
      "Package ${PACKAGE_NAME} declared more than once!")
  endif()

  set(${PACKAGE_NAME}_TRIBITS_PACKAGE_CALLED TRUE)

  tribits_package_decl(${PACKAGE_NAME_IN} ${ARGN})
  tribits_package_def()
endmacro()


#
# @MACRO: tribits_add_test_directories()
#
# Macro called to add a set of test directories for an SE package.
#
# Usage::
#
#    tribits_add_test_directories(<dir1> <dir2> ...)
#
# This macro only needs to be called from the top most ``CMakeLists.txt`` file
# for which all subdirectories are all "tests".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# Currently, all this macro does macro is to call ``add_subdirectory(<diri>)``
# if ``${PACKAGE_NAME}_ENABLE_TESTS`` or
# ``${PARENT_PACKAGE_NAME}_ENABLE_TESTS`` are ``TRUE``. However, this macro
# may be extended in the future in order to modify behavior related to adding
# tests and examples in a uniform way.
#
macro(tribits_add_test_directories)

  if (CURRENTLY_PROCESSING_SUBPACKAGE)

    # This is a subpackage being processed

    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_subpackage() before tribits_add_test_directories()"
        " in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_test_directories() before"
        " tribits_subpackage_postprocess() in"
        " ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

  else()

    # This is a package being processed

    if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_package() or tribits_package_decl() before"
        " tribits_add_test_directories() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

    if(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_test_directories() before "
        " tribits_package_postprocess() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

  endif()

  if(${PACKAGE_NAME}_ENABLE_TESTS OR ${PARENT_PACKAGE_NAME}_ENABLE_TESTS)
    foreach(TEST_DIR ${ARGN})
      tribits_trace_file_processing(PACKAGE  ADD_SUBDIR
        "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/CMakeLists.txt")
      add_subdirectory(${TEST_DIR})
    endforeach()
  endif()

endmacro()


#
# Macros to add common options to add to an SE package
#


#
# @MACRO: tribits_add_debug_option()
#
# Add the standard cache variable option ``${PACKAGE_NAME}_ENABLE_DEBUG`` for
# the package.
#
# Usage::
#
#   tribits_add_debug_option()
#
# This option is given the default ``${${PROJECT_NAME}_ENABLE_DEBUG}`` and if
# true, will set the variable ``HAVE_${PACKAGE_NAME_UC}_DEBUG`` (to be used in
# the package's configured header file).  This macro is typically called in
# the package's `<packageDir>/CMakeLists.txt`_ file.
#
macro(tribits_add_debug_option)
  tribits_add_option_and_define(
    ${PACKAGE_NAME}_ENABLE_DEBUG
    HAVE_${PACKAGE_NAME_UC}_DEBUG
    "Enable a host of runtime debug checking."
    ${${PROJECT_NAME}_ENABLE_DEBUG}
    )
endmacro()


macro(tribits_add_enable_teuchos_time_monitor_option)
  option(
    ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
     "Enable Teuchos time monitors for package ${PACKAGE_NAME}"
    ${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}
    )
endmacro()


#
# @MACRO: tribits_add_show_deprecated_warnings_option()
#
# Add the standard option ``${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS`` for the
# package.
#
# Usage::
#
#   tribits_add_show_deprecated_warnings_option()
#
# This macro should be called in the package's <packageDir>/CMakeLists.txt`_
# file.  This option is given the default value
# ``${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}``.  This option is then looked
# for in `tribits_configure_file()`_ to add macros to add deprecated warnings
# to deprecated parts of a package.
#
macro(tribits_add_show_deprecated_warnings_option)
  advanced_set(
    ${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS  ${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}
    CACHE BOOL
    "Show warnings about deprecated code in ${PACKAGE_NAME}"
    )
  advanced_set(
    ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE  ${${PROJECT_NAME}_HIDE_DEPRECATED_CODE}
    CACHE BOOL
    "Fully exclude deprecated code in ${PACKAGE_NAME}"
    )
endmacro()


macro(tribits_add_explicit_instantiation_option)
  tribits_add_option_and_define(
    ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
    HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
    "Enable the use of explicit template instantiation."
    ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
    )
endmacro()


macro(tribits_add_eti_support)
  append_global_set(${PROJECT_NAME}_ETI_PACKAGES ${PACKAGE_NAME})
  global_null_set(${PACKAGE_NAME}_ETI_LIBRARYSET)
endmacro()


#
# @MACRO: tribits_add_example_directories()
#
# Macro called to conditionally add a set of example directories for an SE
# package.
#
# Usage::
#
#    tribits_add_example_directories(<dir1> <dir2> ...)
#
# This macro typically is called from the top-level
# `<packageDir>/CMakeLists.txt`_ file for which all subdirectories are all
# "examples" according to standard package layout.
#
# This macro can be called several times within a package as desired to break
# up example directories any way one would like.
#
# Currently, all it does macro does is to call ``add_subdirectory(<diri>)`` if
# ``${PACKAGE_NAME}_ENABLE_EXAMPLES`` or
# ``${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES`` are true. However, this macro may
# be extended in the future in order to modify behavior related to adding
# tests and examples in a uniform way.
#
macro(tribits_add_example_directories)

  if (CURRENTLY_PROCESSING_SUBPACKAGE)

    # This is a subpackage being processed
    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_subpackage() before tribits_add_example_directories()"
        " in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_example_directories() before "
        " tribits_subpackage_postprocess() in ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()

  else()

    # This is a package being processed
    if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_package() or tribits_package_decl() before"
        " tribits_add_example_directories() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

    if(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Must call tribits_add_example_directories() before "
        " tribits_package_postprocess() in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
    endif()

  endif()

  if(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES)
    foreach(EXAMPLE_DIR ${ARGN})
      tribits_trace_file_processing(PACKAGE  ADD_SUBDIR
        "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_DIR}/CMakeLists.txt")
      add_subdirectory(${EXAMPLE_DIR})
    endforeach()
  endif()

endmacro()


#
# Utility function that sets up package linkage linkage variables in case the
# package has no libraries.
#

function(tribits_package_finalize_dependency_vars)

  if(${PACKAGE_NAME}_SUBPACKAGES)

    # A package with subpackages should get all of its dependency vars from
    # its enabled subpackages.

    set(PARENT_PACKAGE_INCLUDE_DIRS)
    set(PARENT_PACKAGE_LIBRARY_DIRS)
    set(PARENT_PACKAGE_LIBRARIES)

    set(SUBPACKAGE_IDX 0)
    foreach(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      set(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      set(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      if (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        prepend_set(PARENT_PACKAGE_INCLUDE_DIRS
          ${${SUBPACKAGE_FULLNAME}_INCLUDE_DIRS})
        prepend_set(PARENT_PACKAGE_LIBRARY_DIRS
          ${${SUBPACKAGE_FULLNAME}_LIBRARY_DIRS})
        prepend_set(PARENT_PACKAGE_LIBRARIES
          ${${SUBPACKAGE_FULLNAME}_LIBRARIES})
      endif()

      math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    endforeach()

    if (PARENT_PACKAGE_INCLUDE_DIRS)
      list(REMOVE_DUPLICATES PARENT_PACKAGE_INCLUDE_DIRS)
    endif()
    if (PARENT_PACKAGE_LIBRARY_DIRS)
      list(REMOVE_DUPLICATES PARENT_PACKAGE_LIBRARY_DIRS)
    endif()
    # NOTE: Above, in the rare case that none of the subpackages contain any
    # libraries or any include directories, we need to not call
    # list(REMOVE_DUPLICATES ...).

    # NOTE: There can't be any duplicate libraries in PARENT_PACKAGE_LIBRARIES
    # so no need to remove them.

    global_set(${PACKAGE_NAME}_INCLUDE_DIRS "${PARENT_PACKAGE_INCLUDE_DIRS}")
    global_set(${PACKAGE_NAME}_LIBRARY_DIRS "${PARENT_PACKAGE_LIBRARY_DIRS}")
    global_set(${PACKAGE_NAME}_LIBRARIES "${PARENT_PACKAGE_LIBRARIES}")

  elseif(NOT ${PACKAGE_NAME}_INCLUDE_DIRS)

    # No libraries have been defined for this package so we are going to set
    # them based on this package's dependencies.

    tribits_sort_and_append_package_include_and_link_dirs_and_libs(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)

    tribits_sort_and_append_tpl_include_and_link_dirs_and_libs(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)

    get_directory_property(INCLUDE_DIRS_CURRENT  INCLUDE_DIRECTORIES)
    get_directory_property(LIBRARY_DIRS_CURRENT  PACKAGE_LIBRARY_DIRS)

    prepend_global_set(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCLUDE_DIRS_CURRENT})
    prepend_global_set(${PACKAGE_NAME}_LIBRARY_DIRS  ${LIBRARY_DIRS_CURRENT})
    prepend_global_set(${PACKAGE_NAME}_LIBRARIES  ${LINK_LIBS})

  endif()

endfunction()


#
# Helper macro for [SUB]tribits_package_postprocess()
#
macro(tribits_package_postprocess_common)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_PACKAGE_POSTPROCESS_COMMON: ${PACKAGE_NAME}")
    print_var(${PACKAGE_NAME}_INCLUDE_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARY_DIRS)
    print_var(${PACKAGE_NAME}_LIBRARIES)
  endif()

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES OR
    ${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES
    )
    # Create the configure file so external projects can find packages with a
    # call to find_package(<package_name>).  This also creates the
    # Makefile.export.* files.
    tribits_write_package_client_export_files(${PACKAGE_NAME})
  endif()

  set(${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE TRUE
    CACHE INTERNAL "")

endmacro()


#
# @MACRO: tribits_package_postprocess()
#
# Macro called at the very end of a package's top-level
# `<packageDir>/CMakeLists.txt`_ file that performs some critical
# post-processing activities.
#
# Usage::
#
#   tribits_package_postprocess()
#
# NOTE: It is unfortunate that this macro must be called in a packages's
# top-level ``CMakeLists.txt`` file but limitations of the CMake language make
# it necessary to do so.
#
macro(tribits_package_postprocess)

  # check that this is not being called from inside a subpackage
  if (CURRENTLY_PROCESSING_SUBPACKAGE)
    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Cannot call tribits_package_postprocess() in a subpackage."
        " Use tribits_subpackage_postprocess() instead"
        " ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
    endif()
  endif()

  if(${PACKAGE_NAME}_SUBPACKAGES)
     
    # This is a package that has subpackages
    if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED OR
       NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED OR
       NOT ${PACKAGE_NAME}_TRIBITS_PROCESS_SUBPACKAGES_CALLED )

      tribits_report_invalid_tribits_usage(
	"Must call tribits_package_decl(), tribits_process_subpackages()"
        " and tribits_package_def() before tribits_package_postprocess()."
        "  Because this package has subpackages you cannot use tribits_package()"
        " you must call these in the following order:"
        " tribits_package_decl()"
        " tribits_process_subpackages()"
        " tribits_package_def()"
        " tribits_package_postprocess()"
        " in: "
        "${TRIBITS_PACKAGE_CMAKELIST_FILE}"
        )
    endif()

  else()

    # This is a package without subpackages

    if (
	(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_CALLED)
	AND
	(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED)
      )
      tribits_report_invalid_tribits_usage(
        "Must call tribits_package() or tribits_package_def() before"
	" tribits_package_postprocess()"
	" at the top of the file:\n"
	"  ${TRIBITS_PACKAGE_CMAKELIST_FILE}"
	)
    endif()

  endif()
  
  if(NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED)
    tribits_report_invalid_tribits_usage(
      "Must call tribits_package() before tribits_package_postprocess()." 
      "  Or if your package has subpackages you must first call tribits_package_decl(),"
      " then tribits_process_subpackages(), then tribits_package_def(), then"
      " tribits_package_postprocess() in"
      " ${TRIBITS_PACKAGE_CMAKELIST_FILE}"
      )
  endif()

  if(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
    tribits_report_invalid_tribits_usage(
      "tribits_package_postprocess() has been called more than once in"
      " ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
  endif()

  # Only parent packages have the targets (${PACKAGE_NAME}_libs and
  # (${PACKAGE_NAME}_all
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_PACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    print_var(${PACKAGE_NAME}_LIB_TARGETS)
    print_var(${PACKAGE_NAME}_ALL_TARGETS)
  endif()
  add_custom_target(${PACKAGE_NAME}_libs DEPENDS ${${PACKAGE_NAME}_LIB_TARGETS})
  add_custom_target(${PACKAGE_NAME}_all DEPENDS ${${PACKAGE_NAME}_ALL_TARGETS})

  tribits_package_finalize_dependency_vars()
  tribits_package_postprocess_common()

  if (${PACKAGE_NAME}_SOURCE_DIR STREQUAL ${PROJECT_NAME}_SOURCE_DIR)
    set(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS TRUE)
  else()
    set(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS TRUE PARENT_SCOPE)
  endif()

  set(${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED TRUE)

endmacro()


#
# @MACRO: tribits_process_subpackages()
#
# Macro that processes the `TriBITS Subpackages`_ for a parent `TriBITS
# package`_ for packages that are broken down into subpackages.  This is
# called in the parent packages top-level `<packageDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   tribits_process_subpackages()
#
# This macro must be called after `tribits_package_decl()`_ but before
# `tribits_package_def()`_.
#
macro(tribits_process_subpackages)

  if (CURRENTLY_PROCESSING_SUBPACKAGE)
    tribits_report_invalid_tribits_usage(
      "Cannot call tribits_process_subpackages() in a subpackage."
      " subpackages cannot contain other subpackages"
      " ${CURRENT_SUBPACKAGE_CMAKELIST_FILE}")
  endif()

  if (${PACKAGE_NAME}_TRIBITS_PACKAGE_POSTPROCESS_CALLED)
    tribits_report_invalid_tribits_usage(
      "Must call tribits_process_subpackages() before tribits_package_postprocess()"
      " in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")    
  endif()

  if (NOT ${PACKAGE_NAME}_TRIBITS_PACKAGE_DECL_CALLED)
    tribits_report_invalid_tribits_usage(
      "Must call tribits_package_decl() before tribits_process_subpackages()"
      " in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
  endif()

  if (${PACKAGE_NAME}_TRIBITS_PACKAGE_DEF_CALLED)
    tribits_report_invalid_tribits_usage(
      "Must call tribits_package_def() after tribits_process_subpackages()"
      " in ${TRIBITS_PACKAGE_CMAKELIST_FILE}")
  endif()

  if (NOT ${PARENT_PACKAGE_NAME}_SUBPACKAGES)
    tribits_report_invalid_tribits_usage(
      "The TriBITS Package '${PACKAGE_NAME}' does not have any subpackages."
      "  Therefore, you are not allowed to call tribits_process_subpackages()!")
  endif()

  set(SUBPACKAGE_IDX 0)
  foreach(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

    #message("")
    #print_var(SUBPACKAGE_IDX)
    #print_var(TRIBITS_SUBPACKAGE)

    set(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
    set(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #print_var(SUBPACKAGE_FULLNAME)

    if (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})

      list(GET ${PARENT_PACKAGE_NAME}_SUBPACKAGE_DIRS ${SUBPACKAGE_IDX} SUBPACKAGE_DIR)
      #print_var(SUBPACKAGE_DIR)

      if (NOT ${PROJECT_NAME}_BINARY_DIR STREQUAL ${PARENT_PACKAGE_NAME}_BINARY_DIR)
        dual_scope_set(${SUBPACKAGE_FULLNAME}_BINARY_DIR 
          ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/${SUBPACKAGE_DIR})
      else()
        set(${SUBPACKAGE_FULLNAME}_BINARY_DIR 
          ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/${SUBPACKAGE_DIR})
      endif()

      set(CURRENT_SUBPACKAGE_CMAKELIST_FILE
        "${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}/CMakeLists.txt")
      tribits_trace_file_processing(PACKAGE  ADD_SUBDIR
        ${CURRENT_SUBPACKAGE_CMAKELIST_FILE} )
      set(CURRENTLY_PROCESSING_SUBPACKAGE ${SUBPACKAGE_FULLNAME}) 
      add_subdirectory(${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}
        ${${SUBPACKAGE_FULLNAME}_BINARY_DIR})

    endif()

    math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

  endforeach()
  
        set(CURRENTLY_PROCESSING_SUBPACKAGE FALSE) 
  set(${PACKAGE_NAME}_TRIBITS_PROCESS_SUBPACKAGES_CALLED TRUE)

endmacro()


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
#   ${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES_TO_INSTALL
#   ${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES
#   ${PARENT_PACKAGE_NAME}_LIB_TARGETS
#   ${PARENT_PACKAGE_NAME}_ALL_TARGETS
#
# without carefully studying the documentation in README.DEPENDENCIES and then
# carefully studying all of the code and issues that modify these variables!
#
# ToDo: Write some good unit tests that pin down the behavior of all of this!
#
