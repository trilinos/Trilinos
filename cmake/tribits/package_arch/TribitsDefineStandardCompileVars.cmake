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

INCLUDE(CMakeBuildTypesList)

INCLUDE(AdvancedSet)
INCLUDE(AssertDefined)
INCLUDE(MultilineSet)
INCLUDE(DualScopeSet)
INCLUDE(PrependCmndlineArgs)


#
# Macro that just defines the basic flags
#

MACRO(TRIBITS_DEFINE_STANDARD_COMPILE_FLAGS_VARS  ENABLE_SHADOWING_WARNINGS)

  #
  # Setup and general flags
  #

  SET(GENERAL_BUILD_FLAGS) # Applies to all builds, period
 
  IF (${PROJECT_NAME}_ENABLE_CHECKED_STL)
    PREPEND_CMNDLINE_ARGS(GENERAL_BUILD_FLAGS "-D_GLIBCXX_DEBUG")
  ENDIF()

  SET(DEBUG_SYMBOLS_FLAGS "-g")

  SET(GENERAL_DEBUG_FLAGS "${DEBUG_SYMBOLS_FLAGS} -O0")
 
  SET(GENERAL_RELEASE_FLAGS "-O3")

  IF (${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS)
    PREPEND_CMNDLINE_ARGS(GENERAL_BUILD_FLAGS "${DEBUG_SYMBOLS_FLAGS}")
  ENDIF()

  MULTILINE_SET(C_STRONG_COMPILE_WARNING_FLAGS
    "-ansi" # Check for C89 or C++98 standard code
    " -pedantic" # Adds more static checking to remove non-ANSI GNU extensions
    " -Wall" # Enable a bunch of default warnings
    " -Wno-long-long" # Allow long long int since it is used by MPI, SWIG, etc.
    )
  
  MULTILINE_SET(CXX_STRONG_COMPILE_WARNING_FLAGS
    ${C_STRONG_COMPILE_WARNING_FLAGS}
    " -Wwrite-strings" # Checks for non-const char * copy of string constants
    )

  IF (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS ON)
  ELSEIF (${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS STREQUAL "OFF")
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS OFF)
  ELSE()
    SET(LOCAL_ENABLE_SHADOWING_WARNINGS ${ENABLE_SHADOWING_WARNINGS})
  ENDIF()
  
  IF (LOCAL_ENABLE_SHADOWING_WARNINGS)
    MULTILINE_SET(CXX_STRONG_COMPILE_WARNING_FLAGS
      ${CXX_STRONG_COMPILE_WARNING_FLAGS}
      " -Wshadow" # Warn about general shadowing issues
      " -Woverloaded-virtual" # Warn about hiding virtual functions
      )
  ENDIF()

ENDMACRO()
