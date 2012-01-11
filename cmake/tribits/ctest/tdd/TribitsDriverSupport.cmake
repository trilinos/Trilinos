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

INCLUDE(GetCurrentListDir)
INCLUDE(ParseVariableArguments)
INCLUDE(SetDefaultAndFromEnv)

# Set if the inner CMake installs are performed or not
SET_DEFAULT_AND_FROM_ENV(TDD_FORCE_INNER_CMAKE_INSTALL 1)

# Get this value outside of any functions so it will be the path to
# *this* file and not the path to the file calling any of these
# functions.
GET_CURRENT_LIST_DIR(THIS_SCRIPT_DIR)

#
# Placeholder for any global setup necessary on all machines...
#
# Call this function first at the top of any drivers/MACHINE/CMakeLists.txt
# file.
#
function(TRIBITS_DRIVER_SETUP)
  #
  # But right now... no global setup required...
  #
endfunction()


#
# Add a test to download and install a single "flavor" of CMake.
#
# Known values of cmake_type are 'min', 'release', 'rc' and 'dev'.
#
# This function should be called automatically by
# TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS after many calls to
# TRIBITS_DRIVER_ADD_DASHBOARD have collected a list of which flavors
# of CMake are required to run all the dashboards on this machine.
#
function(TRIBITS_DRIVER_ADD_TEST_THAT_INSTALLS_CMAKE cmake_type)

  set(known 0)
  if("${cmake_type}" STREQUAL "min" OR
     "${cmake_type}" STREQUAL "release" OR
     "${cmake_type}" STREQUAL "rc" OR
     "${cmake_type}" STREQUAL "dev")
    set(known 1)
  endif()
  if(NOT known)
    message("error: unknown --installer-type '${cmake_type}' -- no test added")
    return()
  endif()

  if(NOT TD_BASE_DIR)
    message(FATAL_ERROR "TD_BASE_DIR must be defined before calling this function")
  endif()

  if(NOT TRIBITS_PYTHON_DIR)
    message(FATAL_ERROR "TRIBITS_PYTHON_DIR must be defined before calling this function")
  endif()

  find_program(PYTHON_EXE python)

  SET_DEFAULT_AND_FROM_ENV( TDD_HTTP_PROXY "" )
  IF (TDD_HTTP_PROXY)
    SET(TDD_HTTP_PROXY_ARG "--http-proxy=${TDD_HTTP_PROXY}")
  ELSE()
    SET(TDD_HTTP_PROXY_ARG "")
  ENDIF()

  add_test(uninstall-cmake-${cmake_type} ${CMAKE_COMMAND}
    -E remove_directory "${TD_BASE_DIR}/tools/cmake-${cmake_type}"
    )

  add_test(install-cmake-${cmake_type} ${PYTHON_EXE}
    "${TRIBITS_PYTHON_DIR}/download-cmake.py"
    "--skip-detect"
    "--install-dir=${TD_BASE_DIR}/tools/cmake-${cmake_type}"
    "--installer-type=${cmake_type}"
    "${TDD_HTTP_PROXY_ARG}"
    )

  set_property(TEST install-cmake-${cmake_type}
    PROPERTY DEPENDS "uninstall-cmake-${cmake_type}")

endfunction()


#
# Add a test that runs a dashboard script using a known flavor of ctest.
#
# Required arguments:
#  testname - name of the test as it appears on the driver dashboard
#
#  scriptname - name of the file in CMAKE_CURRENT_SOURCE_DIR to run as
#               ctest -S script
#
# Optional arguments:
#
#  [CTEST_INSTALLER_TYPE min|release|rc|dev]
#  [ENVIRONMENT var1=value1;var2=value2;var3=value3;...]
#  [PROCESSORS 1]
#  [RUN_SERIAL]
#  [TIMEOUT_MINUTES 180]
#
function(TRIBITS_DRIVER_ADD_DASHBOARD testname scriptname)

  MESSAGE("TRIBITS_DRIVER_ADD_DASHBOARD:  '${testname}'  '${scriptname}' [${ARGN}]")

  # Uncomment this line to see output of below PARSE_ARGUMENTS:
  #
  #set(PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED TRUE)

  PARSE_ARGUMENTS(
    # prefix
    PARSE
    # args
    "CTEST_INSTALLER_TYPE;ENVIRONMENT;TIMEOUT_MINUTES;DEPENDS;REQUIRED_FILES"
    # options
    "RUN_SERIAL"
    # the stuff to parse:
    ${ARGN}
  )

  if(NOT "${PARSE_DEFAULT_ARGS}" STREQUAL "")
    message("warning: unrecognized arguments '${PARSE_DEFAULT_ARGS}' to TRIBITS_DRIVER_ADD_DASHBOARD")
  endif()

  if(PARSE_CTEST_INSTALLER_TYPE)
    set(ctest_type ${PARSE_CTEST_INSTALLER_TYPE})
  else()
    set(ctest_type "release")
  endif()

  if(PARSE_TIMEOUT_MINUTES)
    math(EXPR timeout_seconds "60 * ${PARSE_TIMEOUT_MINUTES}")
  else()
    set(timeout_seconds 10800) # 3 hours...
  endif()

  add_test(${testname} ${CMAKE_COMMAND}
    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
    -D ctest_type=${ctest_type}
    -D scriptname=${scriptname}
    -D TD_BASE_DIR=${TD_BASE_DIR}
    -D testname=${testname}
    -P ${THIS_SCRIPT_DIR}/LocateCTestAndRunScript.cmake
    )

  # This test that runs the dashboard depends on the test that installs its
  # driving ctest:
  #
  set_property(TEST ${testname} PROPERTY DEPENDS "install-cmake-${ctest_type}")
  
  if(PARSE_DEPENDS)
    set_property(TEST ${testname} PROPERTY DEPENDS "${PARSE_DEPENDS}")
  endif()
  
  if(PARSE_REQUIRED_FILES)
    set_property(TEST ${testname} PROPERTY REQUIRED_FILES "${PARSE_REQUIRED_FILES}")
  endif()

  if(PARSE_ENVIRONMENT)
    set_property(TEST ${testname} PROPERTY ENVIRONMENT
      "${PARSE_ENVIRONMENT};TD_CTEST_VERSION_TYPE=${ctest_type}")
  else()
    set_property(TEST ${testname} PROPERTY ENVIRONMENT
      "TD_CTEST_VERSION_TYPE=${ctest_type}")
  endif()

  set_property(TEST ${testname} PROPERTY
    FAIL_REGULAR_EXPRESSION "error: TRIBITS_CTEST_DRIVER_ERROR_QUEUE")

  if(PARSE_PROCESSORS)
    set_property(TEST ${testname} PROPERTY PROCESSORS "${PARSE_PROCESSORS}")
  endif()

  if(PARSE_RUN_SERIAL)
    set_property(TEST ${testname} PROPERTY RUN_SERIAL ON)
  endif()

  set_property(TEST ${testname} PROPERTY TIMEOUT "${timeout_seconds}")

  # Track the required cmake installer types in a global property:
  #
  set_property(GLOBAL APPEND PROPERTY TD_CMAKE_INSTALLER_TYPES "${ctest_type}")

endfunction()


#
# Call this function last at the bottom of any drivers/MACHINE/CMakeLists.txt
# file. It will add tests as needed to download and install all the right
# flavors of ctest required to drive the other tests added with
# TRIBITS_DRIVER_ADD_DASHBOARD.
#
function(TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS)

  IF (TDD_FORCE_INNER_CMAKE_INSTALL STREQUAL 1) 

    get_property(types GLOBAL PROPERTY TD_CMAKE_INSTALLER_TYPES)
  
    if (types)
      list(REMOVE_DUPLICATES types)
      foreach (type ${types})
        MESSAGE(STATUS "Adding CMake install test ${type}")
        TRIBITS_DRIVER_ADD_TEST_THAT_INSTALLS_CMAKE(${type})
      endforeach()
    else()
      message("warning: no cmake install tests are necessary...")
      message("  Put calls to TRIBITS_DRIVER_ADD_DASHBOARD")
      message("  *before* calls to TRIBITS_ADD_REQUIRED_CMAKE_INSTALLS...")
    endif()

  ELSE()

    MESSAGE(STATUS "Skipping CMake install tests because TDD_FORCE_INNER_CMAKE_INSTALL!=1")

  ENDIF()

endfunction()
