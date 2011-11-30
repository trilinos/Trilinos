# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Test, and Integrate System
#                 Copyright (2011) Sandia Corporation
#
#
# Copyright (2011) Sandia Corporation. Under the terms of Contract
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

# Definition of the CMakeLists.txt body used by the TriBITS driver job.

# The base directory is the parent of the binary directory.  
# FIXME This is duplicated information. The base directory is set in
# tdd_driver.py. We need to go through two separate CTest invocations
# to get to this point from there. Since environment variables are the
# only way to pass data to CTest, we need to pass this directory
# through the environment.
get_filename_component(TD_BASE_DIR ${CMAKE_BINARY_DIR} PATH)

# CMake versions greater than 2.8.4 have CMAKE_CURRENT_LIST_DIR. Make
# sure it's defined.
if( NOT DEFINED CMAKE_CURRENT_LIST_DIR )
  get_filename_component( CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH )
endif()

# Locate the TriBITS dependencies.
get_filename_component(TRIBITS_ROOT
  "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_LIST_DIR}
  ${TRIBITS_ROOT}/utils
  ${TRIBITS_ROOT}/package_arch
  ${TRIBITS_ROOT}/config_tests
  )

set(TRIBITS_PYTHON_DIR "${TRIBITS_ROOT}/python")

include(CTest)
include(TribitsDriverSupport)

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/CTestCustom.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake
  )

# Function to make exactly one add_subdirectory call based on the site
# name of the machine we're running on. By default, the subdirectory
# is taken to be the site name. The environment variable
# TDD_DRIVER_SUBDIRECTORY can be used to override the default value.
FUNCTION(TDD_PROJECT)
  site_name(site)
  set(subdir "${site}")
  message("site='${site}'")

  # But if that directory does not exist as named, and there's a regex match
  # with the name of a subdirectory, use the exact subdirectory name instead:
  #
  file(GLOB filesAndDirs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*")

  set(dirs "")
  foreach(dir ${filesAndDirs})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
      set(dirs ${dirs} "${dir}")
    endif()
  endforeach()

  foreach(dir ${dirs})
    if("${site}" MATCHES "${dir}")
      set(subdir "${dir}")
      message("site='${site}' MATCHES directory name dir='${dir}'")
    endif()
  endforeach()

  # Allow an environment variable to override the test directory.
  set_default_and_from_env(TDD_DRIVER_SUBDIRECTORY "${subdir}")

  # The one add_subdirectory call:
  #
  message("TDD_DRIVER_SUBDIRECTORY='${TDD_DRIVER_SUBDIRECTORY}'")

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}")
    message(FATAL_ERROR "error: there is no subdirectory of ${CMAKE_CURRENT_SOURCE_DIR} matching '${TDD_DRIVER_SUBDIRECTORY}'")
  endif()

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}/CMakeLists.txt")
    message(FATAL_ERROR "error: there is no CMakeLists.txt file in '${CMAKE_CURRENT_SOURCE_DIR}/${TDD_DRIVER_SUBDIRECTORY}'")
  endif()

  add_subdirectory("${TDD_DRIVER_SUBDIRECTORY}")
ENDFUNCTION()