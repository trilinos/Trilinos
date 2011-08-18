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


MESSAGE(
 "\n***"
 "\n*** TrilinosDriverDashboard.cmake"
 "\n***\n"
 )

# Get the base diretory for the Trilinos source.  We only assume that the
# CTest script that is being called is under Trilinos/cmake.
STRING(REGEX MATCH "(.+/Trilinos)/cmake" TRILINOS_CMAKE_DIR
  "${CTEST_SCRIPT_DIRECTORY}" )
IF("${TRILINOS_CMAKE_DIR}" STREQUAL "")
  STRING(REGEX MATCH "(.+)/cmake" TRILINOS_CMAKE_DIR
    "${CTEST_SCRIPT_DIRECTORY}" )
ENDIF()
  
MESSAGE("TRILINOS_CMAKE_DIR = ${TRILINOS_CMAKE_DIR}")

SET( CMAKE_MODULE_PATH
   "${TRILINOS_CMAKE_DIR}/utils"
   )

#MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

INCLUDE(PrintVar)
INCLUDE(SetDefaultAndFromEnv)

#
# A) Set up the environment get options
#

set(CTEST_SITE "$ENV{CTEST_SITE}")
if("${CTEST_SITE}" STREQUAL "")
  site_name(CTEST_SITE)
endif()
if("${CTEST_SITE}" STREQUAL "")
  if(WIN32)
    string(TOLOWER "$ENV{COMPUTERNAME}" CTEST_SITE)
  else()
    execute_process(COMMAND uname -n
      OUTPUT_VARIABLE CTEST_SITE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()
endif()

set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME}")
if("${CTEST_BUILD_NAME}" STREQUAL "")
  if(WIN32)
    set(HOST_TYPE $ENV{OS})
  else()
    execute_process(COMMAND uname
      OUTPUT_VARIABLE HOST_TYPE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()
  set(CTEST_BUILD_NAME "${HOST_TYPE}-TDD-${CTEST_SITE}")
endif()

set(CTEST_CMAKE_GENERATOR "$ENV{CTEST_CMAKE_GENERATOR}")
if("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  if(WIN32)
    set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  else()
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  endif()
endif()

set(CTEST_TEST_TIMEOUT "$ENV{CTEST_TEST_TIMEOUT}")
if("${CTEST_TEST_TIMEOUT}" STREQUAL "")
  set(CTEST_TEST_TIMEOUT 7200)
endif()
  
# Submit the results to the dashboard or not
SET_DEFAULT_AND_FROM_ENV( TDD_DO_SUBMIT TRUE )

# Dashboard model : Nightly, Experimental, Continuous
SET_DEFAULT_AND_FROM_ENV( TDD_CTEST_TEST_TYPE Experimental )

# set this to ON if you need to test something before committing.
SET_DEFAULT_AND_FROM_ENV( TDD_IN_TESTING_MODE OFF )

get_filename_component(CTEST_SOURCE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}" ABSOLUTE)

get_filename_component(CTEST_UPDATE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../.." ABSOLUTE)

get_filename_component(CTEST_BINARY_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../../TDD_BUILD" ABSOLUTE)

get_filename_component(CTEST_NOTES_FILES
  "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" ABSOLUTE)
if(NOT "$ENV{TDD_CRON_DRIVER_LOGFILE}" STREQUAL "")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} "$ENV{TDD_CRON_DRIVER_LOGFILE}")
endif()
if(NOT "$ENV{TDD_CRON_DRIVER_SCRIPT}" STREQUAL "")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} "$ENV{TDD_CRON_DRIVER_SCRIPT}")
endif()

set(parallel_level "$ENV{TDD_PARALLEL_LEVEL}")
if("${parallel_level}" STREQUAL "")
  set(parallel_level 1)
endif()

set(git_exe "$ENV{TDD_GIT_EXE}")
if("${git_exe}" STREQUAL "")
  set(git_exe "git_exe-NOTFOUND")
  find_program(git_exe NAMES git.cmd eg git)
endif()
if(git_exe)
  set(CTEST_UPDATE_TYPE "git")
  set(CTEST_UPDATE_COMMAND "${git_exe}")
endif()

#
# Run the outer dashboard
#

message("\nA) Empty out ${CTEST_BINARY_DIRECTORY} ...")
ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")

ctest_start("${TDD_CTEST_TEST_TYPE}")

message("\nB) Update ${CTEST_UPDATE_DIRECTORY} ...")
message("      CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
message("      CTEST_UPDATE_TYPE='${CTEST_UPDATE_TYPE}'")
if(NOT TDD_IN_TESTING_MODE)
  ctest_update(SOURCE "${CTEST_UPDATE_DIRECTORY}")
else()
  message("In testing mode no update will be performed.")
endif()

message("\nC) Configure ${CTEST_BINARY_DIRECTORY} ...")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}")

ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")

message("\nD) Build ${CTEST_BINARY_DIRECTORY} ...")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)

message("\nE) Submitting update configure notes build ...")
if (TDD_DO_SUBMIT)
  if(NOT "$ENV{TDD_CTEST_DROP_SITE}" STREQUAL "")
    set(CTEST_DROP_SITE "$ENV{TDD_CTEST_DROP_SITE}")
  endif()
  if(NOT "$ENV{TDD_CTEST_DROP_LOCATION}" STREQUAL "")
    set(CTEST_DROP_LOCATION "$ENV{TDD_CTEST_DROP_LOCATION}")
  endif()
  ctest_submit(PARTS update configure notes build)
else()
  message("\nSkipping submit!")
endif()

message("\nF) Run tests (which run all everything really): PARALLEL_LEVEL ${parallel_level} from ${CTEST_BINARY_DIRECTORY} ...")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${parallel_level})

message("\nG) Submitting Test ...")
if (TDD_DO_SUBMIT)
  ctest_submit(PARTS Test)
else()
  message("\nSkipping submit!")
endif()

MESSAGE(
 "\n*** Finished TrilinosDriverDashboard.cmake ***\n"
 )
