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
# Wrapper used for unit testing purposes
#

MACRO(EXTRAREPO_EXECUTE_PROCESS_WRAPPER)
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    EXECUTE_PROCESS(${ARGN}
      RESULT_VARIABLE  EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL)
    IF (NOT EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL STREQUAL "0")
      MESSAGE(SEND_ERROR
        "Error: EXECUTE_PROCESS(${ARGN}) returned"
	" '${EXTRAREPO_EXECUTE_PROCESS_WRAPPER_RTN_VAL}'")
    ENDIF()
  ELSE()
    MESSAGE("EXECUTE_PROCESS(${ARGN})")
  ENDIF()
ENDMACRO()

#
# Function for getting the tracking branch
#
FUNCTION(EXTRAREPO_GET_TRACKING_BRANCH  EXTRAREPO_SRC_DIR  TRACKING_BRANCH_OUT)
  IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    EXECUTE_PROCESS(
      COMMAND "${GIT_EXE}" rev-parse --abbrev-ref --symbolic-full-name @{u}
      WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE  EP_RTN
      OUTPUT_VARIABLE  TRACKING_BRANCH
      )
    IF (NOT EP_RTN STREQUAL "0")
      MESSAGE(SEND_ERROR "Error: obtaining tracking branch for repo"
        " '${EXTRAREPO_SRC_DIR}' failed!" )
    ENDIF()
  ELSE()
    # For unit testing purpose, just return generic tacking branch string
    SET(TRACKING_BRANCH "tracking/branch")
  ENDIF()
  SET(${TRACKING_BRANCH_OUT}  ${TRACKING_BRANCH}  PARENT_SCOPE)
ENDFUNCTION()


#
# Update an existing git repo
#
FUNCTION(EXTRAREPO_CLEAN_FETCH_RESET  GIT_EXE  EXTRAREPO_SRC_DIR)

  SET(EXTRAREPO_CLEAN_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.clean.out")
  SET(EXTRAREPO_FETCH_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.fetch.out")
  SET(EXTRAREPO_RESET_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.reset.out")

  EXTRAREPO_GET_TRACKING_BRANCH("${EXTRAREPO_SRC_DIR}"
    EXTRAREPO_TRACKING_BRANCH)
  #PRINT_VAR(EXTRAREPO_TRACKING_BRANCH)
  SET(CLEAN_CMND_ARGS
    COMMAND "${GIT_EXE}" clean -fdx
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_CLEAN_OUT_FILE}" )
  SET(FETCH_CMND_ARGS
    COMMAND "${GIT_EXE}" fetch
    TIMEOUT 60 # seconds
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_FETCH_OUT_FILE}" )
  SET(RESET_CMND_ARGS
    COMMAND "${GIT_EXE}" reset --hard "${EXTRAREPO_TRACKING_BRANCH}"
    WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
    OUTPUT_FILE "${EXTRAREPO_RESET_OUT_FILE}" )

  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${CLEAN_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${FETCH_CMND_ARGS})
  EXTRAREPO_EXECUTE_PROCESS_WRAPPER(${RESET_CMND_ARGS})

ENDFUNCTION()
