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
# Helper macros
#

MACRO(TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS LANG BUILDTYPE)

  #PRINT_VAR(CMAKE_${LANG}_FLAGS_${BUILDTYPE})
  #PRINT_VAR(GENERAL_${BUILDTYPE}_FLAGS)

  # Set the default CMAKE_${LANG}_FLAGS_${BUILDTYPE} to empty to override the
  # default that CMake gives!  We want to define these for ourselves.
  SET(CMAKE_${LANG}_FLAGS_${BUILDTYPE} "")

  IF (GENERAL_${BUILDTYPE}_FLAGS)
    DUAL_SCOPE_PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS_${BUILDTYPE}
      "${GENERAL_${BUILDTYPE}_FLAGS}")
  ENDIF()
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding ${LANG} ${BUILDTYPE} flags"
      " \"${GENERAL_${BUILDTYPE}_FLAGS}\"")
    PRINT_VAR(CMAKE_${LANG}_FLAGS_${BUILDTYPE})
  ENDIF()

ENDMACRO()


MACRO(TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS_OVERRIDE LANG BUILDTYPE)

  SET(CMAKE_${LANG}_FLAGS_${BUILDTYPE}_OVERRIDE  ""
    CACHE  STRING
    "If set to non-empty, will override CMAKE_${LANG}_FLAGS_${BUILDTYPE}")

  IF (CMAKE_${LANG}_FLAGS_${BUILDTYPE}_OVERRIDE)
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(
        "-- Overriding CMAKE_${LANG}_FLAGS_${BUILDTYPE}:\n"
        "--  from\n"
        "--    ${CMAKE_${LANG}_FLAGS_${BUILDTYPE}}")
    ENDIF()
    DUAL_SCOPE_SET(CMAKE_${LANG}_FLAGS_${BUILDTYPE}
      ${CMAKE_${LANG}_FLAGS_${BUILDTYPE}_OVERRIDE})
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(
        "--  to\n"
        "--    ${CMAKE_${LANG}_FLAGS_${BUILDTYPE}}")
    ENDIF()
  ENDIF()

  # NOTE: Above, we can't just let users set CMAKE_${LANG}_FLAGS_${BUILDTYPE}
  # in the cache because CMake overrides it.  We have ot add this override
  # option to allow them to override it.

ENDMACRO()


MACRO(TRIBITS_SET_LANGUAGE_ALL_BUILDTYPES_FLAGS_OVERRIDE LANG)

  FOREACH(BUILDTYPE ${CMAKE_BUILD_TYPES_LIST})
    TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS_OVERRIDE(${LANG} ${BUILDTYPE})
  ENDFOREACH()

ENDMACRO()


MACRO(TRIBITS_SET_LANGUAGE_GENERAL_FLAGS LANG)
 
  DUAL_SCOPE_PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS "${GENERAL_BUILD_FLAGS}")
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding general ${LANG} flags \"${${GENERAL_BUILD_FLAGS}}\"")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()

ENDMACRO()


MACRO(TRIBITS_SET_LANGUAGE_COVERAGE_FLAGS LANG)

  DUAL_SCOPE_PREPEND_CMNDLINE_ARGS(CMAKE_${LANG}_FLAGS
   "${COVERAGE_OPTIONS}")
  IF(COVERAGE_OPTIONS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding coverage ${LANG} flags \"${COVERAGE_OPTIONS}\"")
    PRINT_VAR(CMAKE_${LANG}_FLAGS)
  ENDIF()

ENDMACRO()


#
# Function that sets up strong compile options for the primary
# development platform (i.e. gcc)
#
# NOTE: The compiler flags in the cache, which may have been set by
# the user, are not disturbed in this function.  Instead variables in
# the parent base scope are set.
#
# NOTE: This function should be called only once before adding any libraries,
# executables, etc.
#

FUNCTION(TRIBITS_SETUP_BASIC_COMPILE_LINK_FLAGS)

  #
  # Setup and general flags
  #

  TRIBITS_DEFINE_STANDARD_COMPILE_FLAGS_VARS(FALSE)

  #
  # Set up coverge testing options
  #

  IF (${PROJECT_NAME}_ENABLE_COVERAGE_TESTING)
    SET(COVERAGE_OPTIONS "-fprofile-arcs -ftest-coverage")
  ELSE()
    SET(COVERAGE_OPTIONS "")
  ENDIF()

  #
  # C compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C CMAKE_C_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_C AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
    TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS(C DEBUG)
    TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS(C RELEASE)
    TRIBITS_SET_LANGUAGE_GENERAL_FLAGS(C)
    TRIBITS_SET_LANGUAGE_COVERAGE_FLAGS(C)
  ENDIF()
  TRIBITS_SET_LANGUAGE_ALL_BUILDTYPES_FLAGS_OVERRIDE(C)

  #
  # C++ compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX CMAKE_CXX_COMPILER_ID)
  IF (${PROJECT_NAME}_ENABLE_CXX AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS(CXX DEBUG)
    TRIBITS_SET_LANGUAGE_BUILDTYPE_FLAGS(CXX RELEASE)
    TRIBITS_SET_LANGUAGE_GENERAL_FLAGS(CXX)
    TRIBITS_SET_LANGUAGE_COVERAGE_FLAGS(CXX)
  ENDIF()
  TRIBITS_SET_LANGUAGE_ALL_BUILDTYPES_FLAGS_OVERRIDE(CXX)

  #
  # Fortran compiler options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_Fortran AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    TRIBITS_SET_LANGUAGE_COVERAGE_FLAGS(Fortran)
  ENDIF()
  TRIBITS_SET_LANGUAGE_ALL_BUILDTYPES_FLAGS_OVERRIDE(Fortran)

  #
  # Linker options
  #

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_COVERAGE_TESTING COVERAGE_OPTIONS)
  IF (${PROJECT_NAME}_ENABLE_COVERAGE_TESTING AND COVERAGE_OPTIONS)
    DUAL_SCOPE_PREPEND_CMNDLINE_ARGS(CMAKE_EXE_LINKER_FLAGS
     "${COVERAGE_OPTIONS} ${CMAKE_EXE_LINKER_FLAGS}")
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Adding coverage linker flags flags \"${COVERAGE_OPTIONS}\"")
      PRINT_VAR(CMAKE_EXE_LINKER_FLAGS)
    ENDIF()
  ENDIF()
  
ENDFUNCTION()
