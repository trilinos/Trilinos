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

INCLUDE(PrintVar)
INCLUDE(AppendStringVar)
INCLUDE(Join)


FUNCTION(DELETE_CREATE_WORKING_DIRECTORY  WORKING_DIR_IN)
  IF (EXISTS "${WORKING_DIR_IN}")
    MESSAGE("Removing existing working directory"
      " '${WORKING_DIR_IN}'\n")
    IF (NOT SHOW_COMMANDS_ONLY)
      FILE(REMOVE_RECURSE "${WORKING_DIR_IN}")
    ENDIF()
  ENDIF()
  MESSAGE("Creating new working directory"
    " '${WORKING_DIR_IN}'\n")
  IF (NOT SHOW_COMMANDS_ONLY)
    FILE(MAKE_DIRECTORY "${WORKING_DIR_IN}")
  ENDIF()
ENDFUNCTION()


FUNCTION(DRIVE_ADVANCED_TEST)

  SET(ADVANDED_TEST_SEP
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

  SET(TEST_SEP
    "================================================================================")

  SET(OUTPUT_SEP
    "--------------------------------------------------------------------------------")
  
  MATH(EXPR LAST_CMND_IDX ${NUM_CMNDS}-1)

  MESSAGE("\n${ADVANDED_TEST_SEP}\n")
  MESSAGE("Advanced Test: ${TEST_NAME}\n")

  IF (OVERALL_WORKING_DIRECTORY)
    IF (NOT IS_ABSOLUTE OVERALL_WORKING_DIRECTORY)
      SET(OVERALL_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${OVERALL_WORKING_DIRECTORY})
    ENDIF()
    DELETE_CREATE_WORKING_DIRECTORY("${OVERALL_WORKING_DIRECTORY}")
    SET(BASE_WORKING_DIRECTORY "${OVERALL_WORKING_DIRECTORY}")
  ELSE()
    SET(BASE_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  ENDIF()

  FOREACH ( CMND_IDX RANGE ${LAST_CMND_IDX} )
    IF (CMND_IDX EQUAL 0)
      SET(TEST_NAMES_STR "TEST_0")
    ELSE()
      APPEND_STRING_VAR( TEST_NAMES_STR ", TEST_${CMND_IDX}" )
    ENDIF()
  ENDFOREACH()
  MESSAGE("Running test commands: ${TEST_NAMES_STR}") 

  SET(OVERALL_TEST_PASSED TRUE)
  
  FOREACH ( CMND_IDX RANGE ${LAST_CMND_IDX} )

    MESSAGE("\n${TEST_SEP}\n")
    MESSAGE("TEST_${CMND_IDX}\n")

    IF (TEST_${CMND_IDX}_MESSAGE)
      MESSAGE("${TEST_${CMND_IDX}_MESSAGE}\n")
    ENDIF()

    IF (TEST_${CMND_IDX}_WORKING_DIRECTORY)
      IF (NOT IS_ABSOLUTE TEST_${CMND_IDX}_WORKING_DIRECTORY)
        SET(TEST_${CMND_IDX}_WORKING_DIRECTORY ${BASE_WORKING_DIRECTORY}/${TEST_${CMND_IDX}_WORKING_DIRECTORY})
      ENDIF()
      DELETE_CREATE_WORKING_DIRECTORY("${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
    ENDIF()

    JOIN( TEST_CMND_STR " " TRUE ${TEST_${CMND_IDX}_CMND} )
    MESSAGE("Running: ${TEST_CMND_STR}\n")
    SET(EXEC_CMND COMMAND ${TEST_${CMND_IDX}_CMND})

    SET(WORKING_DIR_SET)
    IF (TEST_${CMND_IDX}_WORKING_DIRECTORY)
      SET(WORKING_DIR_SET "${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
    ELSEIF(OVERALL_WORKING_DIRECTORY)
      SET(WORKING_DIR_SET "${OVERALL_WORKING_DIRECTORY}")
    ENDIF()

    IF (WORKING_DIR_SET)
      MESSAGE("  Running in working directory \"${WORKING_DIR_SET}\"\n")
      SET(EXEC_CMND ${EXEC_CMND}
        WORKING_DIRECTORY "${WORKING_DIR_SET}"
        )
      SET(WORKING_DIR "${WORKING_DIR_SET}")
    ELSE()
      SET(WORKING_DIR "${CMAKE_CURRENT_BINARY_DIR}")
    ENDIF()

    IF (TEST_${CMND_IDX}_OUTPUT_FILE)
      IF (NOT IS_ABSOLUTE ${TEST_${CMND_IDX}_OUTPUT_FILE})
        SET(OUTPUT_FILE_USED "${WORKING_DIR}/${TEST_${CMND_IDX}_OUTPUT_FILE}")
      ELSE()
        SET(OUTPUT_FILE_USED "${TEST_${CMND_IDX}_OUTPUT_FILE}")
      ENDIF()
      MESSAGE("  Writing output to file \"${OUTPUT_FILE_USED}\"\n")
      SET(EXEC_CMND ${EXEC_CMND}
        OUTPUT_FILE "${OUTPUT_FILE_USED}"
        )
    ELSE()
      SET(EXEC_CMND ${EXEC_CMND}
        OUTPUT_VARIABLE TEST_CMND_OUT
        ERROR_VARIABLE TEST_CMND_OUT
        )
    ENDIF()

    IF (NOT SHOW_COMMANDS_ONLY)
      
      #PRINT_VAR(EXEC_CMND)

      # Provide the test configuration in an environment variable.
      IF(TEST_CONFIG)
        SET(ENV{TEST_CONFIG} "${TEST_CONFIG}")
      ENDIF(TEST_CONFIG)

      EXECUTE_PROCESS(${EXEC_CMND} RESULT_VARIABLE EXEC_RESULT)

      IF (TEST_${CMND_IDX}_OUTPUT_FILE)
        FILE(READ "${OUTPUT_FILE_USED}" TEST_CMND_OUT)
      ENDIF()

      MESSAGE("${OUTPUT_SEP}\n")

      IF (NOT TEST_${CMND_IDX}_NO_ECHO_OUTPUT)
        MESSAGE("${TEST_CMND_OUT}")
      ELSE()
        MESSAGE("NO_ECHO_OUTPUT\n")
      ENDIF()

    ELSE()

      MESSAGE("\n*** Not running command on request ***")

    ENDIF()

    MESSAGE("${OUTPUT_SEP}\n")

    IF (NOT SHOW_COMMANDS_ONLY)

      MESSAGE("TEST_${CMND_IDX}: Return code = ${EXEC_RESULT}")

      SET(TEST_PASSED FALSE)
      IF (TEST_${CMND_IDX}_PASS_ANY)
        MESSAGE("TEST_${CMND_IDX}: Pass criteria = Pass Any")
        SET(TEST_PASSED TRUE)
      ELSEIF (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION)
        MESSAGE("TEST_${CMND_IDX}: Pass criteria ="
          " Match REGEX {${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION}}")
        STRING(REGEX MATCH "${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION}"
          MATCH_STR "${TEST_CMND_OUT}" )
        IF (MATCH_STR)
          SET(TEST_PASSED TRUE)
        ELSE()
          SET(TEST_PASSED FALSE)
        ENDIF()
      ELSEIF (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL)
        SET(TEST_PASSED TRUE)
        FOREACH(REGEX_STR ${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL})
          STRING(REGEX MATCH "${REGEX_STR}" MATCH_STR "${TEST_CMND_OUT}" )
          SET(MSG_BEGIN "TEST_${CMND_IDX}: Pass criteria = Match REGEX {${REGEX_STR}}")
          IF (MATCH_STR)
            MESSAGE("${MSG_BEGIN} [PASSED]")
          ELSE()
            MESSAGE("${MSG_BEGIN} [FAILED]")
            SET(TEST_PASSED FALSE)
          ENDIF()
        ENDFOREACH()
      ELSE()
        MESSAGE("TEST_${CMND_IDX}: Pass criteria = Return code")
        IF (EXEC_RESULT EQUAL 0)
          SET(TEST_PASSED TRUE)
        ELSE()
          SET(TEST_PASSED FALSE)
        ENDIF()

      ENDIF()

      IF (TEST_PASSED)
        MESSAGE("TEST_${CMND_IDX}: Result = PASSED")
      ELSE()
        MESSAGE("TEST_${CMND_IDX}: Result = FAILED")
        SET(OVERALL_TEST_PASSED FALSE)
        IF (FAIL_FAST)
          MESSAGE("TEST_${CMND_IDX}: FAIL FAST, SKIPPING REST OF TEST CASES!")
          BREAK()
        ENDIF()
      ENDIF()
    ENDIF()
  
  ENDFOREACH()
        
  MESSAGE("\n${TEST_SEP}\n")

  IF (NOT SHOW_COMMANDS_ONLY)
    IF (OVERALL_TEST_PASSED)
      MESSAGE("OVERALL FINAL RESULT: TEST PASSED")
    ELSE()
      MESSAGE("OVERALL FINAL RESULT: TEST FAILED")
    ENDIF()
  ELSE()
    MESSAGE("OVERALL FINAL RESULT: DID NOT RUN COMMANDS")
  ENDIF()

  MESSAGE("\n${ADVANDED_TEST_SEP}\n")

ENDFUNCTION()
