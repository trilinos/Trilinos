# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# CMake -P script that prints out a python datastructure for the
# info for extra repos read from a file
#

#
# A) Get the input
#

# Echo input commandline
message("*** Generate a Python datastructure containing TriBITS/git repos ...")
if (NOT SUPPRESS_PRINT_VAR_OUTPUT)
  message("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
  message("TRIBITS_BASE_DIR = ${TRIBITS_BASE_DIR}")
  message("EXTRA_REPOS_FILE = ${EXTRA_REPOS_FILE}")
  message("EXTRA_REPOS = ${EXTRA_REPOS}")
  message("EXTRA_REPOS_PYTHON_OUT_FILE = ${EXTRA_REPOS_PYTHON_OUT_FILE}")
  message("ENABLE_KNOWN_EXTERNAL_REPOS_TYPE = ${ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
  message("IGNORE_MISSING_EXTRA_REPOSITORIES = ${IGNORE_MISSING_EXTRA_REPOSITORIES}")
  message("CHECK_EXTRAREPOS_EXIST = ${CHECK_EXTRAREPOS_EXIST}")
endif()

if ("${CHECK_EXTRAREPOS_EXIST}"  STREQUAL "")
  set(CHECK_EXTRAREPOS_EXIST  TRUE)
endif()

# Set up necessary variables
include(${PROJECT_SOURCE_DIR}/ProjectName.cmake)
if (NOT SUPPRESS_PRINT_VAR_OUTPUT)
  message("PROJECT_NAME = ${PROJECT_NAME}")
endif()
set(${PROJECT_NAME}_TRIBITS_DIR ${TRIBITS_BASE_DIR})
set(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE ${ENABLE_KNOWN_EXTERNAL_REPOS_TYPE})
set(${PROJECT_NAME}_EXTRAREPOS_FILE ${EXTRA_REPOS_FILE})
set(${PROJECT_NAME}_EXTRA_REPOSITORIES ${EXTRA_REPOS})
set(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES ${IGNORE_MISSING_EXTRA_REPOSITORIES})
set(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST  ${CHECK_EXTRAREPOS_EXIST})
#message("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

#set(${PROJECT_NAME}_VERBOSE_CONFIGURE TRUE)

#
# B) Include files from TriBITS
#

include("${CMAKE_CURRENT_LIST_DIR}/../core/common/TribitsCMakePolicies.cmake"  NO_POLICY_SCOPE)

set( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  )
include(Split)
include(AppendStringVar)
include(SetDefaultAndFromEnv) # Used in ExtraRepositoriesList.cmake file?
include(TribitsProcessExtraRepositoriesList)

# Need to split this argument
split("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","
  ${PROJECT_NAME}_EXTRA_REPOSITORIES)

if (NOT ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
  set(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE "Continuous")
endif()

#
# C) Read in and process the extra repos list variable and process the list
#

tribits_get_and_process_extra_repositories_lists()

#
# D) Write the python dictionary/string that will be evaluated by checkin-test.py
#

set(EXTRA_REPOS_PYTHON_STRING)

append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
  "[\n")

if ("${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES}"  STREQUAL ""
  AND  ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
  )
  set(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
    ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT})
endif()

set(EXTRAREPO_IDX 0)
foreach(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

  # Extract the data for current extra repo
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
    EXTRAREPO_DIR )
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
    EXTRAREPO_REPOTYPE )
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
    EXTRAREPO_REPOURL )
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
    EXTRAREPO_HASPKGS )
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_IDX}
    EXTRAREPO_PREPOST )
  list(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES ${EXTRAREPO_IDX}
    EXTRAREPO_CATEGORY )

  # Write the dictorary entries for this extra rep
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "{" )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'NAME' : '${EXTRAREPO_NAME}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'DIR' : '${EXTRAREPO_DIR}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'REPOTYPE' : '${EXTRAREPO_REPOTYPE}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'REPOURL' : '${EXTRAREPO_REPOURL}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'HASPKGS' : '${EXTRAREPO_HASPKGS}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'PREPOST' : '${EXTRAREPO_PREPOST}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "'CATEGORY' : '${EXTRAREPO_CATEGORY}', " )
  append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
    "},\n" )

  math(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

endforeach()

append_string_var_ext(EXTRA_REPOS_PYTHON_STRING
  "]\n")

#
# E) Write the generated python list/dictionary string to the output file
#
# NOTE: This must be the only output from this script since it gets evaluated.
#
# NOTE: You could use message(STATUS ...) to print to stdout but then it would
# also print an annoying '-- ' at the beginning which might not be porable to just
# remove it.  CMake is not a great general scripting language :-(
#

message("")
if (EXTRA_REPOS_PYTHON_OUT_FILE)
  message("Writing Python datastructure to ${EXTRA_REPOS_PYTHON_OUT_FILE} ...")
  file(WRITE ${EXTRA_REPOS_PYTHON_OUT_FILE} "${EXTRA_REPOS_PYTHON_STRING}")
else()
  message("*** Extra Repositories Python Dictionary")
  message("${EXTRA_REPOS_PYTHON_STRING}")
endif()
