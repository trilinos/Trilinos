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

#
# CMake -P script that prints out a python datastructure for the
# info for extra repos read from a file
#

#
# A) Get the input
#

# Echo input commandline
MESSAGE("*** Generate a Python datastructure containing TriBITS/git repos ...")
IF (NOT SUPPRESS_PRINT_VAR_OUTPUT)
  MESSAGE("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
  MESSAGE("TRIBITS_BASE_DIR = ${TRIBITS_BASE_DIR}")
  MESSAGE("EXTRA_REPOS_FILE = ${EXTRA_REPOS_FILE}")
  MESSAGE("EXTRA_REPOS = ${EXTRA_REPOS}")
  MESSAGE("EXTRA_REPOS_PYTHON_OUT_FILE = ${EXTRA_REPOS_PYTHON_OUT_FILE}")
  MESSAGE("ENABLE_KNOWN_EXTERNAL_REPOS_TYPE = ${ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
  MESSAGE("IGNORE_MISSING_EXTRA_REPOSITORIES = ${IGNORE_MISSING_EXTRA_REPOSITORIES}")
  MESSAGE("CHECK_EXTRAREPOS_EXIST = ${CHECK_EXTRAREPOS_EXIST}")
ENDIF()

IF ("${CHECK_EXTRAREPOS_EXIST}"  STREQUAL "")
  SET(CHECK_EXTRAREPOS_EXIST  TRUE)
ENDIF()

# Set up necessary variables
INCLUDE(${PROJECT_SOURCE_DIR}/ProjectName.cmake)
IF (NOT SUPPRESS_PRINT_VAR_OUTPUT)
  MESSAGE("PROJECT_NAME = ${PROJECT_NAME}")
ENDIF()
SET(${PROJECT_NAME}_TRIBITS_DIR ${TRIBITS_BASE_DIR})
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE ${ENABLE_KNOWN_EXTERNAL_REPOS_TYPE})
SET(${PROJECT_NAME}_EXTRAREPOS_FILE ${EXTRA_REPOS_FILE})
SET(${PROJECT_NAME}_EXTRA_REPOSITORIES ${EXTRA_REPOS})
SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES ${IGNORE_MISSING_EXTRA_REPOSITORIES})
SET(${PROJECT_NAME}_CHECK_EXTRAREPOS_EXIST  ${CHECK_EXTRAREPOS_EXIST})
#MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = ${${PROJECT_NAME}_TRIBITS_DIR}")

#SET(${PROJECT_NAME}_VERBOSE_CONFIGURE TRUE)

#
# B) Include files from TriBITS
#

SET( CMAKE_MODULE_PATH
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/utils"
  "${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch"
  )
INCLUDE(TribitsCMakePolicies)
INCLUDE(Split)
INCLUDE(AppendStringVar)
INCLUDE(SetDefaultAndFromEnv) # Used in ExtraRepositoriesList.cmake file?
INCLUDE(TribitsProcessExtraRepositoriesList)

# Need to split this argument
SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  ","
  ${PROJECT_NAME}_EXTRA_REPOSITORIES)

IF (NOT ${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE)
  SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE "Continuous")
ENDIF()

#
# C) Read in and process the extra repos list variable and process the list
#

TRIBITS_GET_AND_PROCESS_EXTRA_REPOSITORIES_LISTS()

#
# D) Write the python dictionary/string that will be evaluated by checkin-test.py
#

SET(EXTRA_REPOS_PYTHON_STRING)

APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
  "[\n")

IF ("${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES}"  STREQUAL ""
  AND  ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS
  )
  SET(${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES
    ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DEFAULT})
ENDIF()

SET(EXTRAREPO_IDX 0)
FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES})

  # Extract the data for current extra repo
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX}
    EXTRAREPO_DIR )
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_VCTYPES ${EXTRAREPO_IDX}
    EXTRAREPO_REPOTYPE )
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX}
    EXTRAREPO_REPOURL )
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_HASPKGS ${EXTRAREPO_IDX}
    EXTRAREPO_HASPKGS )
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_PREPOSTS ${EXTRAREPO_IDX}
    EXTRAREPO_PREPOST )
  LIST(GET ${PROJECT_NAME}_ALL_EXTRA_REPOSITORIES_CATEGORIES ${EXTRAREPO_IDX}
    EXTRAREPO_CATEGORY )

  # Write the dictorary entries for this extra rep
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "{" )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'NAME' : '${EXTRAREPO_NAME}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'DIR' : '${EXTRAREPO_DIR}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'REPOTYPE' : '${EXTRAREPO_REPOTYPE}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'REPOURL' : '${EXTRAREPO_REPOURL}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'HASPKGS' : '${EXTRAREPO_HASPKGS}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'PREPOST' : '${EXTRAREPO_PREPOST}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "'CATEGORY' : '${EXTRAREPO_CATEGORY}', " )
  APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
    "},\n" )

  MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")

ENDFOREACH()

APPEND_STRING_VAR_EXT(EXTRA_REPOS_PYTHON_STRING
  "]\n")

#
# E) Write the generated python list/dictionary string to the output file
#
# NOTE: This must be the only output from this script since it gets evaluated.
#
# NOTE: You could use MESSAGE(STATUS ...) to print to stdout but then it would
# also print an annoying '-- ' at the beginning which might not be porable to just
# remove it.  CMake is not a great general scripting language :-(
#

MESSAGE("")
IF (EXTRA_REPOS_PYTHON_OUT_FILE)
  MESSAGE("Writing Python datastructure to ${EXTRA_REPOS_PYTHON_OUT_FILE} ...")
  FILE(WRITE ${EXTRA_REPOS_PYTHON_OUT_FILE} "${EXTRA_REPOS_PYTHON_STRING}")
ELSE()
  MESSAGE("*** Extra Repositories Python Dictionary")
  MESSAGE("${EXTRA_REPOS_PYTHON_STRING}")
ENDIF()
