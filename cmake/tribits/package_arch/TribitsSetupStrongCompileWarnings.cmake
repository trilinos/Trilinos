# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

INCLUDE(TribitsDefineStandardCompileVars)
INCLUDE(DualScopePrependCmndlineArgs)


#
# Helper macros and functions
#


MACRO(TRIBITS_SET_LANGUAGE_STRONG_WARNING_FLAGS LANG)

  #MESSAGE("Entering TRIBITS_SET_LANGUAGE_STRONG_WARNING_FLAGS(${LANG})")
  #PRINT_VAR(${PROJECT_NAME}_ENABLE_STRONG_${LANG}_COMPILE_WARNINGS)
  
  DUAL_SCOPE_PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS
    "${${LANG}_STRONG_COMPILE_WARNING_FLAGS}")
  
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding strong ${LANG} warning flags \"${${LANG}_STRONG_COMPILE_WARNING_FLAGS}\"")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()

ENDMACRO()


FUNCTION(TRIBITS_SETUP_STRONG_COMPILE_WARNINGS  ENABLE_SHADOWING_WARNINGS)

  #MESSAGE("Entering TRIBITS_SETUP_STRONG_COMPILE_WARNINGS(${ENABLE_SHADOWING_WARNINGS})")

  #
  # Setup and general flags
  #

  TRIBITS_DEFINE_STANDARD_COMPILE_FLAGS_VARS(${ENABLE_SHADOWING_WARNINGS})

  #
  # C compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C CMAKE_C_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_C
    AND CMAKE_C_COMPILER_ID STREQUAL "GNU"
    AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS
    )
    TRIBITS_SET_LANGUAGE_STRONG_WARNING_FLAGS(C)
  ENDIF()

  #
  # C++ compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX CMAKE_CXX_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_CXX
    AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS
    )
    TRIBITS_SET_LANGUAGE_STRONG_WARNING_FLAGS(CXX)
  ENDIF()

  #
  # Fortran compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_Fortran AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # ToDo: Add Fortran warnings?
  ENDIF()
  
ENDFUNCTION()
