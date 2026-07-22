# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(TribitsSetupStrongCompileWarnings)
include(PrependCmndlineArgs)
include(DualScopeAppendCmndlineArgs)


#
# Helper macros and functions
#


macro(tribits_apply_warnings_as_error_flags_lang LANG)
  prepend_cmndline_args(CMAKE_${LANG}_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}")
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message(STATUS "Setting up for ${LANG} warnings as errors just in this package ...")
    print_var(CMAKE_${LANG}_FLAGS)
  endif()
endmacro()


macro(tribits_set_package_compiler_lang_flags LANG)

  #message("Entering tribits_set_package_compiler_lang_flags(${LANG})")
  #print_var(${PROJECT_NAME}_ENABLE_STRONG_${LANG}_COMPILE_WARNINGS)

  if (${PACKAGE_NAME}_${LANG}_FLAGS)
    dual_scope_append_cmndline_args(CMAKE_${LANG}_FLAGS
      "${${PACKAGE_NAME}_${LANG}_FLAGS}")
  endif()

endmacro()


function(tribits_print_package_compiler_lang_flags LANG SUFFIX)
    message("-- " "${PACKAGE_NAME}: CMAKE_${LANG}_FLAGS${SUFFIX}=\"${CMAKE_${LANG}_FLAGS${SUFFIX}}\"")
endfunction()


# Function that appends package-specific compiler flags for each language
#
function(tribits_append_package_specific_compiler_flags)

  #message("Entering tribits_append_package_specific_compiler_flags() for ${PACKAGE_NAME}")

  # C compiler options
  assert_defined(${PROJECT_NAME}_ENABLE_C CMAKE_C_COMPILER_ID)
  if (${PROJECT_NAME}_ENABLE_C)
    tribits_set_package_compiler_lang_flags(C)
  endif()

  # C++ compiler options
  assert_defined(${PROJECT_NAME}_ENABLE_CXX CMAKE_CXX_COMPILER_ID)
  if (${PROJECT_NAME}_ENABLE_CXX)
    tribits_set_package_compiler_lang_flags(CXX)
  endif()

  # Fortran compiler options
  assert_defined(${PROJECT_NAME}_ENABLE_Fortran)
  if (${PROJECT_NAME}_ENABLE_Fortran)
    tribits_set_package_compiler_lang_flags(Fortran)
  endif()

endfunction()


# Function that prints out all of the compiler flags for a package
#
function(tribits_print_package_compiler_flags)

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE OR ${PROJECT_NAME}_PRINT_PACKAGE_COMPILER_FLAGS)

    string(TOUPPER "${CMAKE_BUILD_TYPE}" upperBuildType)
    set(buildNameSuffix "_${upperBuildType}")

    # C compiler options
    if (${PROJECT_NAME}_ENABLE_C)
      tribits_print_package_compiler_lang_flags(C "")
      tribits_print_package_compiler_lang_flags(C ${buildNameSuffix})
    endif()

    # C++ compiler options
    if (${PROJECT_NAME}_ENABLE_CXX)
      tribits_print_package_compiler_lang_flags(CXX "")
      tribits_print_package_compiler_lang_flags(CXX ${buildNameSuffix})
    endif()

    # Fortran compiler options
    if (${PROJECT_NAME}_ENABLE_Fortran)
      tribits_print_package_compiler_lang_flags(Fortran "")
      tribits_print_package_compiler_lang_flags(Fortran ${buildNameSuffix})
    endif()

  endif()

endfunction()


# Macro that sets up compiler flags for a top-level package (not subpackage)
#
# This CMake code is broken out in order to allow it to be unit tested.
#
macro(tribits_setup_compiler_flags  PACKAGE_NAME_IN)

  # Set up strong warning flags

  if (NOT PARSE_DISABLE_STRONG_WARNINGS AND NOT ${PACKAGE_NAME_IN}_DISABLE_STRONG_WARNINGS)
    tribits_setup_strong_compile_warnings(${PARSE_ENABLE_SHADOWING_WARNINGS})
  endif()

  # Set up for warnings as errors if requested

  assert_defined(PARSE_CLEANED)

  assert_defined(${PROJECT_NAME}_ENABLE_C ${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS)
  if (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS)
    tribits_apply_warnings_as_error_flags_lang(C)
  endif()

  assert_defined(${PROJECT_NAME}_ENABLE_CXX ${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS)
  if (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS)
    tribits_apply_warnings_as_error_flags_lang(CXX)
  endif()

  tribits_append_package_specific_compiler_flags()

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("Final package compiler flags:")
  endif()
  tribits_print_package_compiler_flags()

endmacro()
