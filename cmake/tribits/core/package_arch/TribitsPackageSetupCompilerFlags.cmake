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


INCLUDE(TribitsSetupStrongCompileWarnings)
INCLUDE(PrependCmndlineArgs)
INCLUDE(DualScopeAppendCmndlineArgs)


#
# Helper macros and functions
#


MACRO(TRIBITS_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG LANG)
  PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}")
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Setting up for ${LANG} warnings as errors just in this package ...")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()
ENDMACRO()


MACRO(TRIBITS_SET_PACKAGE_LANGUAGE_FLAGS LANG)

  #MESSAGE("Entering TRIBITS_SET_PACKAGE_LANGUAGE_FLAGS(${LANG})")
  #PRINT_VAR(${PROJECT_NAME}_ENABLE_STRONG_${LANG}_COMPILE_WARNINGS)

  IF (${PACKAGE_NAME}_${LANG}_FLAGS)
    DUAL_SCOPE_APPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS
      "${${PACKAGE_NAME}_${LANG}_FLAGS}")
  ENDIF()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding strong ${LANG} warning flags \"${${LANG}_STRONG_COMPILE_WARNING_FLAGS}\"")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()

ENDMACRO()


FUNCTION(TRIBITS_SETUP_ADD_PACKAGE_COMPILE_FLAGS)

  #MESSAGE("Entering TRIBITS_SETUP_ADD_PACKAGE_COMPILE_FLAGS()")

  #
  # C compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C CMAKE_C_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_C)
    TRIBITS_SET_PACKAGE_LANGUAGE_FLAGS(C)
  ENDIF()

  #
  # C++ compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX CMAKE_CXX_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_CXX)
    TRIBITS_SET_PACKAGE_LANGUAGE_FLAGS(CXX)
  ENDIF()

  #
  # Fortran compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_Fortran)
    TRIBITS_SET_PACKAGE_LANGUAGE_FLAGS(Fortran)
  ENDIF()

ENDFUNCTION()













#
# Macro that sets up compiler flags for a package
#
# This CMake code is broken out in order to allow it to be unit tested.
#

MACRO(TRIBITS_SETUP_COMPILER_FLAGS  PACKAGE_NAME_IN)

  # Set up strong warning flags

  IF (NOT PARSE_DISABLE_STRONG_WARNINGS AND NOT ${PACKAGE_NAME_IN}_DISABLE_STRONG_WARNINGS)
    TRIBITS_SETUP_STRONG_COMPILE_WARNINGS(${PARSE_ENABLE_SHADOWING_WARNINGS})
  ENDIF()

  # Set up for warnings as errors if requested

  ASSERT_DEFINED(PARSE_CLEANED)

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C ${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS)
    TRIBITS_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG(C)
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX ${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS)
    TRIBITS_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG(CXX)
  ENDIF()

  # Append package specific options
  TRIBITS_SETUP_ADD_PACKAGE_COMPILE_FLAGS()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Final compiler flags:")
    PRINT_VAR(CMAKE_CXX_FLAGS)
    PRINT_VAR(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
    PRINT_VAR(CMAKE_C_FLAGS)
    PRINT_VAR(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE})
    PRINT_VAR(CMAKE_Fortran_FLAGS)
    PRINT_VAR(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE})
  ENDIF()

ENDMACRO()
