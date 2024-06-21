# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsDefineStandardCompileVars)
include(DualScopePrependCmndlineArgs)


#
# Helper macros and functions
#


macro(tribits_set_language_strong_warning_flags LANG)

  #message("Entering tribits_set_language_strong_warning_flags(${LANG})")
  #print_var(${PROJECT_NAME}_ENABLE_STRONG_${LANG}_COMPILE_WARNINGS)

  dual_scope_prepend_cmndline_args(CMAKE_${LANG}_FLAGS
    "${${LANG}_STRONG_COMPILE_WARNING_FLAGS}")

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message(STATUS "Adding strong ${LANG} warning flags \"${${LANG}_STRONG_COMPILE_WARNING_FLAGS}\"")
    print_var(CMAKE_${LANG}_FLAGS)
  endif()

endmacro()


function(tribits_setup_strong_compile_warnings  ENABLE_SHADOWING_WARNINGS)

  #message("Entering tribits_setup_strong_compile_warnings(${ENABLE_SHADOWING_WARNINGS})")

  #
  # Setup and general flags
  #

  tribits_define_standard_compile_flags_vars(${ENABLE_SHADOWING_WARNINGS})

  #
  # C compiler options
  #

  assert_defined(${PROJECT_NAME}_ENABLE_C CMAKE_C_COMPILER_ID)
  if (${PROJECT_NAME}_ENABLE_C
    AND CMAKE_C_COMPILER_ID STREQUAL "GNU"
    AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS
    )
    tribits_set_language_strong_warning_flags(C)
  endif()

  #
  # C++ compiler options
  #

  assert_defined(${PROJECT_NAME}_ENABLE_CXX CMAKE_CXX_COMPILER_ID)
  if (${PROJECT_NAME}_ENABLE_CXX
    AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS
    )
    tribits_set_language_strong_warning_flags(CXX)
  endif()

  #
  # Fortran compiler options
  #

  assert_defined(${PROJECT_NAME}_ENABLE_Fortran)
  if (${PROJECT_NAME}_ENABLE_Fortran AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # ToDo: Add Fortran warnings?
  endif()

endfunction()
