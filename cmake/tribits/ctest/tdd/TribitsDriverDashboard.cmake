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


MESSAGE(
 "\n***"
 "\n*** TribitsDriverDashboard.cmake"
 "\n***\n"
 )

# Ge the directly where this file lives in the TriBITS tree.  We use this
# to figure out where everything in in the TriBITS directory tree.
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Get the Tribits base directory
SET(TRIBITS_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
GET_FILENAME_COMPONENT(TRIBITS_ROOT "${TRIBITS_ROOT}" ABSOLUTE)
MESSAGE("TRIBITS_ROOT = '${TRIBITS_ROOT}'")

# Get the directory containing the TriBITS CMake utilities using this
# script's location as the reference point.
SET(TRIBITS_CMAKE_UTILS_DIR "${TRIBITS_ROOT}/utils")
MESSAGE("TRIBITS_CMAKE_UTILS_DIR = '${TRIBITS_CMAKE_UTILS_DIR}'")

SET( CMAKE_MODULE_PATH
  "${TRIBITS_CMAKE_UTILS_DIR}"
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

#
# Set CTEST_BUILD_NAME from TDD_BUILD_NAME in env or set default.
#
# NOTE: CTEST_BUILD_NAME is a built-in CTest varaible and therefore it
# should not be set from the environment since it will give crosstalk
# with TribitsCTestDriverCore.cmake.
#
set(CTEST_BUILD_NAME "$ENV{TDD_BUILD_NAME}")
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
  
# Extra directories to pull updates from
SET_DEFAULT_AND_FROM_ENV( TDD_EXTRA_GIT_PULL_DIRS "" )

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

#
# Allow environment variables to override default values for the
# source, update, and binary directories.
#

get_filename_component(CTEST_SOURCE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}" ABSOLUTE)
SET_DEFAULT_AND_FROM_ENV(CTEST_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY})

get_filename_component(CTEST_UPDATE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../.." ABSOLUTE)
SET_DEFAULT_AND_FROM_ENV(CTEST_UPDATE_DIRECTORY ${CTEST_UPDATE_DIRECTORY})

get_filename_component(CTEST_BINARY_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../../../TDD_BUILD" ABSOLUTE)
SET_DEFAULT_AND_FROM_ENV(CTEST_BINARY_DIRECTORY ${CTEST_BINARY_DIRECTORY})

SET(CTEST_NOTES_FILES)
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
  find_program(git_exe NAMES git.cmd git)
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

if (NOT TDD_IN_TESTING_MODE)

  ctest_update(SOURCE "${CTEST_UPDATE_DIRECTORY}")

  foreach(EXTRA_PULL_DIR ${TDD_EXTRA_GIT_PULL_DIRS})
    SET(EXTRA_PULL_DIR_ABS "${CTEST_UPDATE_DIRECTORY}/${EXTRA_PULL_DIR}")
    SET(PULL_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRA_PULL_DIR}.pull.out")
    set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${PULL_OUT_FILE}) 
    MESSAGE("Pull extra updates in '${EXTRA_PULL_DIR_ABS}' ...")
    execute_process(
      COMMAND ${git_exe} pull
      WORKING_DIRECTORY "${EXTRA_PULL_DIR_ABS}"
      TIMEOUT 60 # seconds
      OUTPUT_FILE "${PULL_OUT_FILE}"
      ERROR_FILE "${PULL_OUT_FILE}"
      RESULT_VARIABLE PULL_RESULT
      )
    IF (NOT PULL_RESULT STREQUAL 0)
      MESSAGE(SEND_ERROR
        "The pull of '${EXTRA_PULL_DIR}' failed with error code ${PULL_RESULT}")
    ENDIF()
  endforeach()

else()

  message("\nTesting mode: no updates of outer sources being performed!")

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
  message("\nSkipping submit on request!")
endif()

message("\nF) Run tests (which runs everything really): PARALLEL_LEVEL ${parallel_level} from ${CTEST_BINARY_DIRECTORY} ...")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${parallel_level})

message("\nG) Submitting Test ...")
if (TDD_DO_SUBMIT)
  ctest_submit(PARTS Test)
else()
  message("\nSkipping submit on request!")
endif()

MESSAGE(
 "\n*** Finished TribitsDriverDashboard.cmake ***\n"
 )
