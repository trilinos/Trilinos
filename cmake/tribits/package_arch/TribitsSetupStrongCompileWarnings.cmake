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
