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

MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/package_arch"
  )

INCLUDE(MessageWrapper)
INCLUDE(TribitsSetupBasicCompileLinkFlags)
INCLUDE(TribitsPackageSetupCompilerFlags)
INCLUDE(UnitTestHelpers)
INCLUDE(GlobalSet)


#####################################################################
#
# Unit tests for setting up compiler options
#
#####################################################################


#
# Set up unit test functions that will be called below to actually run the
# unit tests.
#
# The reason that we use functions is so that we can change varibles just
# inside of the functions that have their own variable scoping.  In that way,
# we can keep variables that are set in one unit test from affecting the
# others.
#


MACRO(TRIBITS_SET_ALL_COMPILER_ID ID)
  SET(CMAKE_C_COMPILER_ID ${ID})
  SET(CMAKE_CXX_COMPILER_ID ${ID})
  SET(CMAKE_Fortran_COMPILER_ID ${ID})
ENDMACRO()


MACRO(TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS)
  IF (NOT PACKAGE_NAME)
    SET(PACKAGE_NAME "DummyPackage")
  ENDIF()
  TRIBITS_SETUP_BASIC_COMPILE_LINK_FLAGS()
  TRIBITS_SETUP_COMPILER_FLAGS(${PACKAGE_NAME})
ENDMACRO()


FUNCTION(UNITEST_GCC_BASE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_SHADOW_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with shadow warnings")
  MESSAGE("***\n")

  SET(CMAKE_BUILD_TYPE DEBUG)
  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(PARSE_ENABLE_SHADOWING_WARNINGS TRUE)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings -Wshadow -Woverloaded-virtual" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_GLOBAL_ENABLE_SHADOW_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler with global enable of shadow warnings")
  MESSAGE("***\n")

  SET(CMAKE_BUILD_TYPE DEBUG)
  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS ON)
  SET(PARSE_ENABLE_SHADOWING_WARNINGS OFF)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings -Wshadow -Woverloaded-virtual" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_GLOBAL_DISABLE_SHADOW_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler with global disable of shadow warnings")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS OFF)
  SET(PARSE_ENABLE_SHADOWING_WARNINGS ON)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_COVERAGE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with coverage flags")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_COVERAGE_TESTING ON)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -fprofile-arcs -ftest-coverage" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings -fprofile-arcs -ftest-coverage" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS
    "-fprofile-arcs -ftest-coverage" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_CHECKED_STL_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with checked STL enables")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_CHECKED_STL ON)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long  -D_GLIBCXX_DEBUG" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings  -D_GLIBCXX_DEBUG" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_CLEANED_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options warnings as errors")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(PARSE_CLEANED TRUE)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "--warnings_as_errors_placeholder -ansi -pedantic -Wall -Wno-long-long" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "--warnings_as_errors_placeholder -ansi -pedantic -Wall -Wno-long-long -Wwrite-strings" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_NO_STRONG_WARNINGS_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC without strong warnings")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS FALSE)
  SET(${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS FALSE)
  SET(${PROJECT_NAME}_ENABLE_STRONG_Fortran_COMPILE_WARNINGS FALSE)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_NO_PACKAGE_STRONG_WARNINGS_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC without strong warnings for a package")
  MESSAGE("***\n")

  SET(PACKAGE_NAME "MyPackage")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PACKAGE_NAME}_DISABLE_STRONG_WARNINGS TRUE)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_SHADOW_CLEANED_CHECKED_STL_COVERAGE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with shadow warnings,"
    " warnings as errors, checked stl, and coverage tests")
  MESSAGE("***\n")

  SET(CMAKE_BUILD_TYPE DEBUG)
  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(PARSE_ENABLE_SHADOWING_WARNINGS TRUE)
  SET(${PROJECT_NAME}_ENABLE_COVERAGE_TESTING ON)
  SET(${PROJECT_NAME}_ENABLE_CHECKED_STL ON)
  SET(PARSE_CLEANED TRUE)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "--warnings_as_errors_placeholder -ansi -pedantic -Wall -Wno-long-long -fprofile-arcs -ftest-coverage -D_GLIBCXX_DEBUG" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "--warnings_as_errors_placeholder -ansi -pedantic -Wall -Wno-long-long -Wwrite-strings -Wshadow -Woverloaded-virtual -fprofile-arcs -ftest-coverage -D_GLIBCXX_DEBUG" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "-fprofile-arcs -ftest-coverage" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_ADDITIONAL_USER_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with user amended options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  # These flags stay at the end of CMAKE_<LANG>_FLAGS
  SET(CMAKE_C_FLAGS "--additional-user-c-flags")
  SET(CMAKE_CXX_FLAGS "--additional-user-cxx-flags")
  SET(CMAKE_Fortran_FLAGS "--additional-user-fortran-flags")
  # These flags don't stay in CMAKE_<LANG>_FLAGS_<BUILDTYPE>
  SET(CMAKE_C_FLAGS_DEBUG "--additional-user-c-flags-dbg")
  SET(CMAKE_CXX_FLAGS_DEBUG "--additional-user-cxx-flags-dbg")
  SET(CMAKE_Fortran_FLAGS_DEBUG "--additional-user-fortran-flags-dbg")
  SET(CMAKE_C_FLAGS_RELEASE "--additional-user-c-flags-rel")
  SET(CMAKE_CXX_FLAGS_RELEASE "--additional-user-cxx-flags-rel")
  SET(CMAKE_Fortran_FLAGS_RELEASE "--additional-user-fortran-flags-rel")

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long   --additional-user-c-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings   --additional-user-cxx-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS " --additional-user-fortran-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "--additional-user-fortran-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "--additional-user-fortran-flags-rel" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_USER_DEBUG_OVERRIDE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC with user override of debug options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)

  SET(CMAKE_C_FLAGS_DEBUG_OVERRIDE "--additional-user-c-flags-dbg")
  SET(CMAKE_CXX_FLAGS_DEBUG_OVERRIDE "--additional-user-cxx-flags-dbg")
  SET(CMAKE_Fortran_FLAGS_DEBUG_OVERRIDE "--additional-user-fortran-flags-dbg")

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "--additional-user-c-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "--additional-user-cxx-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "--additional-user-fortran-flags-dbg" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_USER_RELEASE_OVERRIDE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC with user override of release options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)

  SET(CMAKE_C_FLAGS_RELEASE_OVERRIDE "--additional-user-c-flags-rel")
  SET(CMAKE_CXX_FLAGS_RELEASE_OVERRIDE "--additional-user-cxx-flags-rel")
  SET(CMAKE_Fortran_FLAGS_RELEASE_OVERRIDE "--additional-user-fortran-flags-rel")

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "--additional-user-c-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "--additional-user-cxx-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "--additional-user-fortran-flags-rel" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_USER_OVERRIDE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC with fully user overridden options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  # These flags stay at the end of CMAKE_<LANG>_FLAGS
  SET(CMAKE_C_FLAGS "--additional-user-c-flags")
  SET(CMAKE_CXX_FLAGS "--additional-user-cxx-flags")
  SET(CMAKE_Fortran_FLAGS "--additional-user-fortran-flags")
  # These flags don't stay in CMAKE_<LANG>_FLAGS_<BUILDTYPE>
  SET(CMAKE_C_FLAGS_DEBUG "--additional-user-c-flags-dbg")
  SET(CMAKE_CXX_FLAGS_DEBUG "--additional-user-cxx-flags-dbg")
  SET(CMAKE_Fortran_FLAGS_DEBUG "--additional-user-fortran-flags-dbg")
  SET(CMAKE_C_FLAGS_RELEASE "--additional-user-c-flags-rel")
  SET(CMAKE_CXX_FLAGS_RELEASE "--additional-user-cxx-flags-rel")
  SET(CMAKE_Fortran_FLAGS_RELEASE "--additional-user-fortran-flags-rel")
  # Turn off internal options
  #SET(CMAKE_BUILD_TYPE NONE) This just affects what CMake uses
  SET(${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS OFF)
  SET(${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS OFF)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "  --additional-user-c-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "  --additional-user-cxx-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS " --additional-user-fortran-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_NONE "" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_NONE "" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_NONE "" )
  # Since CMAKE_BUILD_TYPE=NONE, these are not actually used
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "--additional-user-fortran-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "--additional-user-fortran-flags-rel" )

ENDFUNCTION()


FUNCTION(UNITEST_GCC_WITH_DEBUG_SYMBOLES_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing GCC base compiler options with debut symboles")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(GNU)
  SET(${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS ON)

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long  -g" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS
    "-ansi -pedantic -Wall -Wno-long-long -Wwrite-strings  -g" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "-O3" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "" )

ENDFUNCTION()


FUNCTION(UNITEST_OTHER_WITH_SHADOW_CLEANED_CHECKED_STL_COVERAGE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing OTHER base compiler options with shadow warnings,"
    " warnings as errors, checked stl, and coverage tests")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(OTHER)
  SET(PARSE_ENABLE_SHADOWING_WARNINGS TRUE)
  SET(${PROJECT_NAME}_ENABLE_COVERAGE_TESTING ON)
  SET(${PROJECT_NAME}_ENABLE_CHECKED_STL ON)
  SET(PARSE_CLEANED TRUE)
  SET(CMAKE_C_FLAGS "--default-c-flags")
  SET(CMAKE_CXX_FLAGS "--default-cxx-flags")
  SET(CMAKE_Fortran_FLAGS "--default-fortran-flags")
  SET(CMAKE_C_FLAGS_DEBUG "--default-c-flags-dbg")
  SET(CMAKE_CXX_FLAGS_DEBUG "--default-cxx-flags-dbg")
  SET(CMAKE_Fortran_FLAGS_DEBUG "--default-fortran-flags-dbg")
  SET(CMAKE_C_FLAGS_RELEASE "--default-c-flags-rel")
  SET(CMAKE_CXX_FLAGS_RELEASE "--default-cxx-flags-rel")
  SET(CMAKE_Fortran_FLAGS_RELEASE "--default-fortran-flags-rel")

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS "--warnings_as_errors_placeholder --default-c-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS "--warnings_as_errors_placeholder --default-cxx-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "--default-fortran-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "--default-c-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "--default-cxx-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "--default-fortran-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "--default-c-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "--default-cxx-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "--default-fortran-flags-rel" )

ENDFUNCTION()


FUNCTION(UNITEST_OTHER_BASE_OPTIONS)

  MESSAGE("\n***")
  MESSAGE("*** Testing OTHER base compiler options")
  MESSAGE("***\n")

  TRIBITS_SET_ALL_COMPILER_ID(OTHER)
  SET(CMAKE_C_FLAGS "--default-c-flags")
  SET(CMAKE_CXX_FLAGS "--default-cxx-flags")
  SET(CMAKE_Fortran_FLAGS "--default-fortran-flags")
  SET(CMAKE_C_FLAGS_DEBUG "--default-c-flags-dbg")
  SET(CMAKE_CXX_FLAGS_DEBUG "--default-cxx-flags-dbg")
  SET(CMAKE_Fortran_FLAGS_DEBUG "--default-fortran-flags-dbg")
  SET(CMAKE_C_FLAGS_RELEASE "--default-c-flags-rel")
  SET(CMAKE_CXX_FLAGS_RELEASE "--default-cxx-flags-rel")
  SET(CMAKE_Fortran_FLAGS_RELEASE "--default-fortran-flags-rel")

  TRIBITS_COMPILE_OPTIONS_COMMON_ACTIONS()

  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS "--default-c-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS "--default-cxx-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS "--default-fortran-flags" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_DEBUG "--default-c-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_DEBUG "--default-cxx-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_DEBUG "--default-fortran-flags-dbg" )
  UNITTEST_COMPARE_CONST( CMAKE_C_FLAGS_RELEASE "--default-c-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_CXX_FLAGS_RELEASE "--default-cxx-flags-rel" )
  UNITTEST_COMPARE_CONST( CMAKE_Fortran_FLAGS_RELEASE "--default-fortran-flags-rel" )

ENDFUNCTION()


# OTHER with shadow, warnings as errors, checked STL, and coverage

# 

# ???


#####################################################################
#
# Execute the unit tests
#
#####################################################################

# Assume that all unit tests will pass by default
GLOBAL_SET(UNITTEST_OVERALL_PASS TRUE)
GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED 0)
GLOBAL_SET(UNITTEST_OVERALL_NUMRUN 0)

# Set common/base options
SET(PROJECT_NAME "DummyProject")
SET(${PROJECT_NAME}_ENABLE_C TRUE)
SET(${PROJECT_NAME}_ENABLE_CXX TRUE)
SET(${PROJECT_NAME}_ENABLE_Fortran TRUE)
SET(${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS TRUE)
SET(${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS TRUE)
SET(${PROJECT_NAME}_ENABLE_Fortran_DEBUG_COMPILE_FLAGS TRUE)
SET(${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS TRUE)
SET(${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS TRUE)
SET(${PROJECT_NAME}_ENABLE_STRONG_Fortran_COMPILE_WARNINGS TRUE)
SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS "--warnings_as_errors_placeholder")
SET(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS "")
SET(${PROJECT_NAME}_ENABLE_COVERAGE_TESTING OFF)
SET(${PROJECT_NAME}_ENABLE_CHECKED_STL OFF)
SET(${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS OFF)
SET(${PROJECT_NAME}_VERBOSE_CONFIGURE TRUE)
SET(PARSE_DISABLE_STRONG_WARNINGS FALSE)
SET(PARSE_ENABLE_SHADOWING_WARNINGS FALSE)
SET(PARSE_CLEANED FALSE)
SET(CMAKE_BUILD_TYPE DEBUG)

UNITEST_GCC_BASE_OPTIONS()
UNITEST_GCC_WITH_SHADOW_OPTIONS()
UNITEST_GCC_GLOBAL_ENABLE_SHADOW_OPTIONS()
UNITEST_GCC_GLOBAL_DISABLE_SHADOW_OPTIONS()
UNITEST_GCC_WITH_COVERAGE_OPTIONS()
UNITEST_GCC_WITH_CHECKED_STL_OPTIONS()
UNITEST_GCC_WITH_CLEANED_OPTIONS()
UNITEST_GCC_NO_STRONG_WARNINGS_OPTIONS()
UNITEST_GCC_NO_PACKAGE_STRONG_WARNINGS_OPTIONS()
UNITEST_GCC_WITH_SHADOW_CLEANED_CHECKED_STL_COVERAGE_OPTIONS()
UNITEST_GCC_ADDITIONAL_USER_OPTIONS()
UNITEST_GCC_USER_DEBUG_OVERRIDE_OPTIONS()
UNITEST_GCC_USER_RELEASE_OVERRIDE_OPTIONS()
UNITEST_GCC_USER_OVERRIDE_OPTIONS()
UNITEST_GCC_WITH_DEBUG_SYMBOLES_OPTIONS()
UNITEST_OTHER_BASE_OPTIONS()
UNITEST_OTHER_WITH_SHADOW_CLEANED_CHECKED_STL_COVERAGE_OPTIONS()

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(144)
