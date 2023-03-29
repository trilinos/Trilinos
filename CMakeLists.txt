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

#
# A) Define your project name and set up major project options
#

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 3.23.0 FATAL_ERROR)

# Must set the project name as a variable at very beginning before including anything else
# We set the project name in a separate file so CTest scripts can use it.
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ProjectName.cmake)

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

# Set up to use ccache
if($ENV{CCACHE_NODISABLE})
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
      if(CMAKE_GENERATOR STREQUAL "Xcode")
        message(FATAL_ERROR "ccache not supported with Xcode in Trilinos")
      else()
        set(CMAKE_C_COMPILER_LAUNCHER   "${CCACHE_PROGRAM}")
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
      endif()
  else()
      message(FATAL_ERROR "ccache requested but not found")
  endif()
endif()

# Set an env so we know we are in configure
set(ENV{CMAKE_IS_IN_CONFIGURE_MODE} 1)

# Don't define TriBITS override of include_directories()
set(TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE TRUE)

#
# B) Pull in the TriBITS system and execute
#

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/cmake/tribits/TriBITS.cmake)

# Make Trilinos create <Package>Config.cmake files by default
SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT ON)
# Make Trilinos set up CPack support by default
SET(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT ON)
# Don't allow disabled subpackages to be excluded from tarball
SET(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT FALSE)

# Set up C++ language standard selection (rest is handled by TriBITS)
SET(${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT 17)
SET(${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED "(17|20|23)")
ADVANCED_SET(CMAKE_CXX_STANDARD ${${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT}
  CACHE STRING
  "C++ standard number with values ${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED} (default ${${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT})")
IF (NOT CMAKE_CXX_STANDARD MATCHES "^${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED}$")
  MESSAGE(FATAL_ERROR
    "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} is not in the allowed set "
    "${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED}")
ENDIF ()

message(STATUS "The Current CXX Standard is : " ${CMAKE_CXX_STANDARD})

# Force off support for Makefile.export.* files while TriBITS is refactored to
# remove this (#8498)
SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES OFF CACHE BOOL
  "Support being removed from TriBITS (see trilinos/Trilinos#8498)" FORCE)

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

INSTALL_BUILD_STATS_SCRIPTS()

# Install TriBITS so that other projects can use it.
ADVANCED_SET(${PROJECT_NAME}_INSTALL_TriBITS ON CACHE BOOL
  "If ture, install TriBITS into <lib-install-dir>/cmake/tribits/")
IF (${PROJECT_NAME}_INSTALL_TriBITS)
  ASSERT_DEFINED(Trilinos_SOURCE_DIR)
  ASSERT_DEFINED(${PROJECT_NAME}_INSTALL_LIB_DIR)
  INSTALL(
    DIRECTORY "${Trilinos_SOURCE_DIR}/cmake/tribits"
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake"
    PATTERN "*.pyc" EXCLUDE
    )
ENDIF()

IF(${PROJECT_NAME}_ENABLE_YouCompleteMe)
  INCLUDE(CodeCompletion)
ENDIF()

message(STATUS "If publishing results using Trilinos, please cite us: https://trilinos.github.io/cite.html")
