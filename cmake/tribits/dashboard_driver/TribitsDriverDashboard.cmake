# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


message(
 "\n***"
 "\n*** TribitsDriverDashboard.cmake"
 "\n***\n"
 )

# Ge the directly where this file lives in the TriBITS tree.  We use this
# to figure out where everything in in the TriBITS directory tree.
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Get the Tribits base directory
set(TRIBITS_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
get_filename_component(TRIBITS_ROOT "${TRIBITS_ROOT}" ABSOLUTE)
message("TRIBITS_ROOT = '${TRIBITS_ROOT}'")

# Get the directory containing the TriBITS CMake utilities using this
# script's location as the reference point.
set(TRIBITS_CMAKE_UTILS_DIR "${TRIBITS_ROOT}/core/utils")
message("TRIBITS_CMAKE_UTILS_DIR = '${TRIBITS_CMAKE_UTILS_DIR}'")

set( CMAKE_MODULE_PATH
  "${TRIBITS_CMAKE_UTILS_DIR}"
   )

#message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

include(PrintVar)
include(SetDefaultAndFromEnv)

include(${CMAKE_CURRENT_LIST_DIR}/../ctest_driver/TribitsCTestDriverCoreHelpers.cmake)


#
# Override CTEST_SUBMIT to drive multiple submits and to detect failed
# submissions and track them as queued errors.
#

macro(tdd_ctest_submit)

  # Cache the original CTEST_DROP_SITE and CTEST_DROP_LOCATION
  if ("${TDD_CTEST_DROP_SITE_ORIG}" STREQUAL "")
    set(TDD_CTEST_DROP_SITE_ORIG ${CTEST_DROP_SITE})
    if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      print_var(TDD_CTEST_DROP_SITE_ORIG)
    endif()
  endif()
  if ("${TDD_CTEST_DROP_LOCATION_ORIG}" STREQUAL "")
    set(TDD_CTEST_DROP_LOCATION_ORIG ${CTEST_DROP_LOCATION})
    if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
      print_var(TDD_CTEST_DROP_LOCATION_ORIG)
    endif()
  endif()

  # Do the first submit
  set(CTEST_DROP_SITE ${TDD_CTEST_DROP_SITE_ORIG})
  set(CTEST_DROP_LOCATION ${TDD_CTEST_DROP_LOCATION_ORIG})
  if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
    print_var(CTEST_DROP_SITE)
    print_var(CTEST_DROP_LOCATION)
  endif()

  ctest_submit(${ARGN})

  # Do the second submit if requested!
  if (TDD_2ND_CTEST_DROP_SITE OR TDD_2ND_CTEST_DROP_LOCATION)

    message("\nDoing submit to second CDash site ...\n")

    if (NOT "${TDD_2ND_CTEST_DROP_SITE}" STREQUAL "")
      if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        print_var(TDD_2ND_CTEST_DROP_SITE)
      endif()
      set(CTEST_DROP_SITE ${TDD_2ND_CTEST_DROP_SITE})
    endif()

    if (NOT "${TDD_2ND_CTEST_DROP_LOCATION}" STREQUAL "")
      if (TRIBITS_CTEST_SUBMIT_DEBUG_DUMP)
        print_var(TDD_2ND_CTEST_DROP_LOCATION)
      endif()
      set(CTEST_DROP_LOCATION ${TDD_2ND_CTEST_DROP_LOCATION})
    endif()

    ctest_submit(${ARGN})

  endif()

endmacro()

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
# NOTE: CTEST_BUILD_NAME is a built-in CTest variable and therefore it
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
set_default_and_from_env( TDD_EXTRA_GIT_PULL_DIRS "" )

set(CTEST_TEST_TIMEOUT "$ENV{CTEST_TEST_TIMEOUT}")
if("${CTEST_TEST_TIMEOUT}" STREQUAL "")
  set(CTEST_TEST_TIMEOUT 7200)
endif()

# Submit the results to the dashboard or not
set_default_and_from_env( TDD_DO_SUBMIT TRUE )

# Dashboard model : Nightly, Experimental, Continuous
set_default_and_from_env( TDD_CTEST_TEST_TYPE Experimental )

# set this to ON if you need to test something before committing.
set_default_and_from_env( TDD_IN_TESTING_MODE OFF )

#
# Allow environment variables to override default values for the
# source, update, and binary directories.
#

get_filename_component(CTEST_SOURCE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}" ABSOLUTE)
set_default_and_from_env(CTEST_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY})

get_filename_component(CTEST_UPDATE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../.." ABSOLUTE)
set_default_and_from_env(CTEST_UPDATE_DIRECTORY ${CTEST_UPDATE_DIRECTORY})

get_filename_component(CTEST_BINARY_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../../../TDD_BUILD" ABSOLUTE)
set_default_and_from_env(CTEST_BINARY_DIRECTORY ${CTEST_BINARY_DIRECTORY})

set(CTEST_NOTES_FILES)
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
    set(EXTRA_PULL_DIR_ABS "${CTEST_UPDATE_DIRECTORY}/${EXTRA_PULL_DIR}")
    set(PULL_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRA_PULL_DIR}.pull.out")
    set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${PULL_OUT_FILE})
    message("Pull extra updates in '${EXTRA_PULL_DIR_ABS}' ...")
    tribits_update_git_extrarepo("${git_exe}" "${EXTRA_PULL_DIR_ABS}")
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
  tdd_ctest_submit(PARTS update configure notes build)
else()
  message("\nSkipping submit on request!")
endif()

message("\nF) Run tests (which runs everything really): PARALLEL_LEVEL ${parallel_level} from ${CTEST_BINARY_DIRECTORY} ...")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${parallel_level})

message("\nG) Submitting Test ...")
if (TDD_DO_SUBMIT)
  tdd_ctest_submit(PARTS Test)
else()
  message("\nSkipping submit on request!")
endif()

message(
 "\n*** Finished TribitsDriverDashboard.cmake ***\n"
 )
