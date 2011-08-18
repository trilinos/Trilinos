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


INCLUDE(PackageArchSetupStrongCompileWarnings)
INCLUDE(PrependCmndlineArgs)



MACRO(PACKAGE_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG LANG)
  PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS}") 
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Setting up for ${LANG} warnings as errors just in this package ...")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()
ENDMACRO()


#
# Macro that sets up compiler flags for a package
#
# This CMake code is broken out in order to allow it to be unit tested.
#

MACRO(PACKAGE_SETUP_COMPILER_FLASGS)

  # Set up strong warning flags

  IF (NOT PARSE_DISABLE_STRONG_WARNINGS)
    PACKAGE_ARCH_SETUP_STRONG_COMPILE_WARNINGS(${PARSE_ENABLE_SHADOWING_WARNINGS})
  ENDIF()

  # Set up for warnings as errors if requested

  ASSERT_DEFINED(PARSE_CLEANED)

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C ${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS)
    PACKAGE_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG(C)
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX ${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS)
    PACKAGE_APPLY_WARNINGS_AS_ERROR_FLAGS_LANG(CXX)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Final compiler flags")
    PRINT_VAR(CMAKE_CXX_FLAGS)
    PRINT_VAR(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
    PRINT_VAR(CMAKE_C_FLAGS)
    PRINT_VAR(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE})
    PRINT_VAR(CMAKE_Fortran_FLAGS)
    PRINT_VAR(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE})
  ENDIF()

ENDMACRO()
