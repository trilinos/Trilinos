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
