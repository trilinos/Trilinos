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

INCLUDE(CheckCXXSourceRuns)

FUNCTION(TRIBITS_ENABLE_CXX11)
  ##
  ## Only try introspection if we have not already
  ## or if the user has not overridden then
  ##
  IF(NOT TRIBITS_CXX11_FLAGS)
     MESSAGE(STATUS "Search for CXX11 compiler flag.")
     INCLUDE(CheckCXXSourceCompiles)
     ##
     ## List of possible compiler flags to use
     ##
     SET(CXX11_FLAG_OPTIONS
         "/Qstd=c++11"   # intel windows
         "-std=c++11"    # intel/clang linux/mac
         "-std=gnu++11"  # gcc
     )
     ##
     ## Same CXX11 source
     ##
     SET(TSOURCE
        "
        #include <vector>
        int main() {
          // check >> closing brackets
          std::vector<std::vector<float>> vecvecfloat(1);
          vecvecfloat[0].resize(1);
          vecvecfloat[0][0] = 0.0f;
          std::vector<int> vec(10);
          auto b        = vec.begin();
          decltype(b) e = vec.end();
          std::fill(b,e,1);
          // two examples, taken from the wikipedia article on C++0x :)
          std::vector<int> some_list;
          int total = 0;
          int value = 5;
          std::for_each(some_list.begin(), some_list.end(), [&total](int x) {
              total += x;
              });
          std::for_each(some_list.begin(), some_list.end(), [&, value](int x) {
              total += x * value;
              });
          return 0;
        }
        "
     )
     ##
     ## Try to compile with each flag
     ##
     FOREACH( TOPTION ${CXX11_FLAG_OPTIONS})
        IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "Testing CXX11 flag: ${TOPTION}")
        ENDIF()
        SET(CMAKE_REQUIRED_FLAGS "${TOPTION}")
        CHECK_CXX_SOURCE_COMPILES("${TSOURCE}" CXX11_FLAGS_COMPILE_RESULT
           ## Some compilers do not fail with a bad flag
           FAIL_REGEX "unrecognized .*option"                     # GNU
           FAIL_REGEX "ignoring unknown option"                   # MSVC
           FAIL_REGEX "warning D9002"                             # MSVC, any lang
           FAIL_REGEX "[Uu]nknown option"                         # HP
           FAIL_REGEX "[Ww]arning: [Oo]ption"                     # SunPro
           FAIL_REGEX "command option .* is not recognized"       # XL
        )
        IF(CXX11_FLAGS_COMPILE_RESULT)
          # CACHE the successful flags
          #SET(CXX11_FLAGS_COMPILE_RESULT "${TOPTION}" CACHE INTERNAL "cxx11 support flag")
          # Append compiler flag to CMAKE_CXX_FLAGS
          SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TOPTION}" CACHE INTERNAL "")
          IF(${PROJECT_NAME}_VERBOSE_CONFIGURE})
            MESSAGE(STATUS "Successful CXX11 flag: '${TOPTION}'")
          ENDIF()
          SET(TRIBITS_CXX11_FLAGS ${TOPTION} PARENT_SCOPE)
          BREAK()
        ELSE()
          SET(CXX11_FLAGS_COMPILE_RESULT TRUE)
        ENDIF()
     ENDFOREACH()
  ELSE(NOT TRIBITS_CXX11_FLAGS)
     ##
     ## We already detected or set the cxx11 flags
     ##
     MESSAGE(STATUS "CXX11 Flags already set: '${TRIBITS_CXX11_FLAGS}'")
  ENDIF(NOT TRIBITS_CXX11_FLAGS)
ENDFUNCTION(TRIBITS_ENABLE_CXX11)

FUNCTION(TRIBITS_CHECK_CXX11_SUPPORT VARNAME)

  # support for >> in addition to > > when closing double templates
  SET(SOURCE_CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS
  "
#include <vector>
int main() {
  // check >> closing brackets
  std::vector<std::vector<float>> vecvecfloat(1);
  vecvecfloat[0].resize(1);
  vecvecfloat[0][0] = 0.0f;
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS}" CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS)

  # support for auto and typedecl()
  SET(SOURCE_CXX11_AUTOTYPEDVARIABLES
  "
#include <vector>
int main() {
  std::vector<int> vec(10);
  auto b        = vec.begin();
  decltype(b) e = vec.end();
  std::fill(b,e,1);
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_AUTOTYPEDVARIABLES}" CXX11_AUTOTYPEDVARIABLES)

  # support for lambda expressions
  SET(SOURCE_CXX11_LAMBDAS
  "
#include <vector>
#include <algorithm>
int main() {
  // two examples, taken from the wikipedia article on C++0x :)
  std::vector<int> some_list;
  int total = 0;
  int value = 5;
  std::for_each(some_list.begin(), some_list.end(), [&total](int x) {
      total += x;
      });
  std::for_each(some_list.begin(), some_list.end(), [&, value](int x) {
      total += x * value;
      });
  // 
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_LAMBDAS}" CXX11_LAMBDAS)

  IF (NOT CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS OR NOT CXX11_AUTOTYPEDVARIABLES OR NOT CXX11_LAMBDAS)
    SET(${VARNAME} FALSE PARENT_SCOPE)
  ELSE()
    SET(${VARNAME} TRUE PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()
