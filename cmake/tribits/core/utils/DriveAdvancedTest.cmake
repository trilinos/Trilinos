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

INCLUDE(PrintVar)
INCLUDE(AppendStringVar)
INCLUDE(Join)
INCLUDE(TimingUtils)
INCLUDE(TribitsGetCategoriesString)


FUNCTION(PRINT_CURRENT_DATE_TIME  PREFIX_STR)
  EXECUTE_PROCESS( COMMAND  date  OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  DATE_TIME )
  MESSAGE("${PREFIX_STR} ${DATE_TIME}\n")
ENDFUNCTION()


FUNCTION(PRINT_UPTIME  PREFIX_STR)
  EXECUTE_PROCESS( COMMAND  uptime  OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  MACHINE_LOAD )
  MESSAGE("${PREFIX_STR} ${MACHINE_LOAD}")
ENDFUNCTION()


FUNCTION(PRINT_SINGLE_CHECK_RESULT  MSG_BEGIN  TEST_CASE_PASSED_IN)
  IF (TEST_CASE_PASSED_IN)
    MESSAGE("${MSG_BEGIN} [PASSED]")
  ELSE()
    MESSAGE("${MSG_BEGIN} [FAILED]")
  ENDIF()
ENDFUNCTION()


FUNCTION(DELETE_CREATE_WORKING_DIRECTORY  WORKING_DIR_IN   SKIP_CLEAN)
  IF (EXISTS "${WORKING_DIR_IN}" AND NOT SKIP_CLEAN)
    MESSAGE("Removing existing working directory"
      " '${WORKING_DIR_IN}'\n")
    IF (NOT SHOW_COMMANDS_ONLY)
      FILE(REMOVE_RECURSE "${WORKING_DIR_IN}")
    ENDIF()
  ENDIF()
  IF (NOT EXISTS "${WORKING_DIR_IN}")
    MESSAGE("Creating new working directory"
      " '${WORKING_DIR_IN}'\n")
    IF (NOT SHOW_COMMANDS_ONLY)
      FILE(MAKE_DIRECTORY "${WORKING_DIR_IN}")
    ENDIF()
  ENDIF()
ENDFUNCTION()


MACRO(SETUP_AND_RUN_TEST_IDX_COPY_FILES_BLOCK)

  MESSAGE(
    "Copy files from:\n"
    "    ${TEST_${CMND_IDX}_SOURCE_DIR}/\n"
    "  to:\n"
    "    ${TEST_${CMND_IDX}_DEST_DIR}/\n")

  SPLIT("${TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR}"
    "," COPY_FILES_TO_TEST_DIR_ARRAY)

  SET(TEST_CASE_PASSED TRUE)

  MESSAGE("${OUTPUT_SEP}\n")

  # Make the dest directory if not exists (cmake -E copy will not create it)
  IF (NOT EXISTS "${TEST_${CMND_IDX}_DEST_DIR}")
    MESSAGE("Creating dest directory ${TEST_${CMND_IDX}_DEST_DIR}/ ...")
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_${CMND_IDX}_DEST_DIR}"
      RESULT_VARIABLE MKDIR_COMMAND_RTN
      )
    IF (NOT MKDIR_COMMAND_RTN EQUAL 0)
      SET(TEST_CASE_PASSED FALSE)
    ENDIF()
    # NOTE: Above, it would have been great to be able use FILE(MAKE_DIRECTORY
    # ...)  but if that fails it will just abortS the cmake -P script.  There
    # is no way to handle that gracefully.  Therefore, we have to use
    # EXECUTE_PROCESS(cmake -E make_directory ...).
  ENDIF()

  # Copy the full list of files to the dest dir
  FOREACH(FILENAME ${COPY_FILES_TO_TEST_DIR_ARRAY})
    MESSAGE("Copy file ${FILENAME} ...") 
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E copy
        "${TEST_${CMND_IDX}_SOURCE_DIR}/${FILENAME}"
        "${TEST_${CMND_IDX}_DEST_DIR}/"
      RESULT_VARIABLE COPY_COMMAND_RTN
      )
    IF (NOT COPY_COMMAND_RTN EQUAL 0)
      SET(TEST_CASE_PASSED FALSE)
   ENDIF()
  ENDFOREACH()

  # NOTE: Above, it would have been great to be able use
  # CONFIGURE_FILE(... COPYONLY) but if that fails it will just abortS the
  # cmake -P script.  There is no way to handle that gracefully.  Therefore,
  # we have to use EXECUTE_PROCESS(cmake -E copy ...).  Also, it would have
  # been great to just copy the files all at once with `cmake -E copy
  # <srcDir>/<file0> <srDir>/<file1> ... <distDir>/` but older versions of
  # CMake don't support that mode.  They only support copying the files one at
  # a time :-(

  MESSAGE("\n${OUTPUT_SEP}\n")

  MESSAGE("") # Need a vertical space for readability of the output

ENDMACRO()


MACRO(SETUP_AND_RUN_TEST_IDX_CMND_BLOCK)

  # Address working directory for this TEST_<IDX> block if set
  IF (TEST_${CMND_IDX}_WORKING_DIRECTORY)
    IF (NOT  IS_ABSOLUTE  "${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
      SET(TEST_${CMND_IDX}_WORKING_DIRECTORY
        ${BASE_WORKING_DIRECTORY}/${TEST_${CMND_IDX}_WORKING_DIRECTORY})
    ENDIF()
    DELETE_CREATE_WORKING_DIRECTORY("${TEST_${CMND_IDX}_WORKING_DIRECTORY}"
      ${TEST_${CMND_IDX}_SKIP_CLEAN_WORKING_DIRECTORY})
  ENDIF()

  # Set up the TEST_<IDX> command block
  JOIN( TEST_CMND_STR " " TRUE ${TEST_${CMND_IDX}_CMND} )
  MESSAGE("Running: ${TEST_CMND_STR}\n")
  SET(EXEC_CMND COMMAND ${TEST_${CMND_IDX}_CMND})

  # Set up the workig directory that this TEST_<IDX> CMND block will run in

  SET(WORKING_DIR_SET)
  IF (TEST_${CMND_IDX}_WORKING_DIRECTORY)
    SET(WORKING_DIR_SET "${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
  ELSEIF(OVERALL_WORKING_DIRECTORY)
    SET(WORKING_DIR_SET "${OVERALL_WORKING_DIRECTORY}")
  ENDIF()

  IF (WORKING_DIR_SET)
    MESSAGE("  Running in working directory \"${WORKING_DIR_SET}\"\n")
    SET(WORKING_DIR "${WORKING_DIR_SET}")
  ELSE()
    SET(WORKING_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  ENDIF()

  # Set the actual command that will be run with EXECUTE_PROCES()

  SET(EXEC_CMND ${EXEC_CMND}
    WORKING_DIRECTORY "${WORKING_DIR}"
    )

  # Set up the optional output file that the EXECUTE_PROCESS() command will write to

  IF (TEST_${CMND_IDX}_OUTPUT_FILE)
    IF (NOT  IS_ABSOLUTE  "${TEST_${CMND_IDX}_OUTPUT_FILE}")
      SET(OUTPUT_FILE_USED "${WORKING_DIR}/${TEST_${CMND_IDX}_OUTPUT_FILE}")
    ELSE()
      SET(OUTPUT_FILE_USED "${TEST_${CMND_IDX}_OUTPUT_FILE}")
    ENDIF()
    MESSAGE("  Writing output to file \"${OUTPUT_FILE_USED}\"\n")
  ENDIF()

  # Run the actual command with EXECUTTE_PROCESS() (or just print what would run) ...

  IF (NOT SHOW_COMMANDS_ONLY)

    # Provide the test configuration in an environment variable.
    IF(TEST_CONFIG)
      SET(ENV{TEST_CONFIG} "${TEST_CONFIG}")
    ENDIF(TEST_CONFIG)

    EXECUTE_PROCESS(
      ${EXEC_CMND}
      OUTPUT_VARIABLE TEST_CMND_OUT
      ERROR_VARIABLE TEST_CMND_OUT
      RESULT_VARIABLE EXEC_RESULT
      )

    IF (TEST_${CMND_IDX}_OUTPUT_FILE)
      FILE(WRITE "${OUTPUT_FILE_USED}" "${TEST_CMND_OUT}")
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

ENDMACRO()


MACRO(DETERMINE_TEST_IDX_CMND_BLOCK_PASS_FAIL)

  MESSAGE("TEST_${CMND_IDX}: Return code = ${EXEC_RESULT}")

  # A) Apply first set of pass/fail logic
  SET(TEST_CASE_PASSED FALSE)
  IF (TEST_${CMND_IDX}_PASS_ANY)
    SET(TEST_CASE_PASSED TRUE)
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = Pass Any"
      ${TEST_CASE_PASSED} )
  ELSEIF (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION)
    STRING(REGEX MATCH "${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION}"
      MATCH_STR "${TEST_CMND_OUT}" )
    IF (MATCH_STR)
      SET(TEST_CASE_PASSED TRUE)
    ELSE()
      SET(TEST_CASE_PASSED FALSE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = Match REGEX {${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION}}"
      ${TEST_CASE_PASSED})
  ELSEIF (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL)
    SET(TEST_CASE_PASSED TRUE)
    FOREACH(REGEX_STR ${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL})
      STRING(REGEX MATCH "${REGEX_STR}" MATCH_STR "${TEST_CMND_OUT}" )
      IF (NOT "${MATCH_STR}" STREQUAL "")
        SET(THIS_REGEX_MATCHED  TRUE)
      ELSE()
        SET(THIS_REGEX_MATCHED  FALSE)
      ENDIF()
      IF (NOT  THIS_REGEX_MATCHED)
        SET(TEST_CASE_PASSED FALSE)
      ENDIF()
      PRINT_SINGLE_CHECK_RESULT(
        "TEST_${CMND_IDX}: Pass criteria = Match REGEX {${REGEX_STR}}"
        ${THIS_REGEX_MATCHED} )
    ENDFOREACH()
  ELSE()
    IF (EXEC_RESULT EQUAL 0)
      SET(TEST_CASE_PASSED TRUE)
    ELSE()
      SET(TEST_CASE_PASSED FALSE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = Zero return code"
      ${TEST_CASE_PASSED} )
  ENDIF()

  # B) Check for failing regex matching?
  IF (TEST_${CMND_IDX}_FAIL_REGULAR_EXPRESSION)
    STRING(REGEX MATCH "${TEST_${CMND_IDX}_FAIL_REGULAR_EXPRESSION}"
      MATCH_STR "${TEST_CMND_OUT}" )
    IF (MATCH_STR)
      SET(TEST_CASE_PASSED FALSE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = Not match REGEX {${TEST_${CMND_IDX}_FAIL_REGULAR_EXPRESSION}}"
     ${TEST_CASE_PASSED} )
  ENDIF()

  # C) Check for return code always 0?
  IF (TEST_${CMND_IDX}_ALWAYS_FAIL_ON_NONZERO_RETURN)
    IF (NOT EXEC_RESULT EQUAL 0)
      SET(ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED FALSE)
      SET(TEST_CASE_PASSED FALSE)
    ELSE()
      SET(ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED TRUE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = ALWAYS_FAIL_ON_NONZERO_RETURN"
      ${ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED} )
  ELSEIF (TEST_${CMND_IDX}_ALWAYS_FAIL_ON_ZERO_RETURN)
    IF (EXEC_RESULT EQUAL 0)
      SET(ALWAYS_FAIL_ON_ZERO_RETURN_PASSED FALSE)
      SET(TEST_CASE_PASSED FALSE)
    ELSE()
      SET(ALWAYS_FAIL_ON_ZERO_RETURN_PASSED TRUE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = ALWAYS_FAIL_ON_ZERO_RETURN"
      ${ALWAYS_FAIL_ON_ZERO_RETURN_PASSED} )
  ENDIF()

  # D) Invert pass/fail result?
  IF (TEST_${CMND_IDX}_WILL_FAIL)
    IF (TEST_CASE_PASSED)
      SET(TEST_CASE_PASSED FALSE)
    ELSE()
      SET(TEST_CASE_PASSED TRUE)
    ENDIF()
    PRINT_SINGLE_CHECK_RESULT(
      "TEST_${CMND_IDX}: Pass criteria = WILL_FAIL (invert the above 'Pass critera')"
      ${TEST_CASE_PASSED} )
  ENDIF()

ENDMACRO()


FUNCTION(DRIVE_ADVANCED_TEST)

  #
  # A) Print the header for the advanced test
  #

  SET(ADVANDED_TEST_SEP
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

  SET(TEST_SEP
    "================================================================================")

  SET(OUTPUT_SEP
    "--------------------------------------------------------------------------------")

  MATH(EXPR LAST_CMND_IDX ${NUM_CMNDS}-1)

  MESSAGE("\n${ADVANDED_TEST_SEP}\n")
  MESSAGE("Advanced Test: ${TEST_NAME}\n")

  MESSAGE("Selected Test/CTest Propeties:")
  TRIBITS_GET_CATEGORIES_STRING("${CATEGORIES}" CATEGORIES_IN_COMMAS)
  MESSAGE("  CATEGORIES = ${CATEGORIES_IN_COMMAS}")
  MESSAGE("  PROCESSORS = ${PROCESSORS}")
  IF (TIMEOUT)
    MESSAGE("  TIMEOUT    = ${TIMEOUT}\n")
  ELSE()
    MESSAGE("  TIMEOUT    = DEFAULT\n")
  ENDIF()

  IF (SHOW_MACHINE_LOAD  AND  NOT  SHOW_COMMANDS_ONLY)
    PRINT_UPTIME("Starting Uptime:")
  ENDIF()

  IF (SHOW_START_END_DATE_TIME  AND  NOT  SHOW_COMMANDS_ONLY)
    PRINT_CURRENT_DATE_TIME("Starting at:")
  ENDIF()

  IF (OVERALL_WORKING_DIRECTORY)
    IF (NOT  IS_ABSOLUTE  "${OVERALL_WORKING_DIRECTORY}")
      SET(OVERALL_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${OVERALL_WORKING_DIRECTORY})
    ENDIF()
    DELETE_CREATE_WORKING_DIRECTORY("${OVERALL_WORKING_DIRECTORY}"
      ${SKIP_CLEAN_OVERALL_WORKING_DIRECTORY})
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

  IF (SHOW_START_END_DATE_TIME AND NOT SHOW_COMMANDS_ONLY)
   TIMER_GET_RAW_SECONDS(TEST_CMND_START)
   SET(TEST_OVERALL_START ${TEST_CMND_START})
  ENDIF()

  #
  # B) Loop over and run the TEST_<IDX> blocks one at a time
  #

  SET(OVERALL_TEST_PASSED TRUE)

  FOREACH ( CMND_IDX RANGE ${LAST_CMND_IDX} )
    MESSAGE("\n${TEST_SEP}\n")
    MESSAGE("TEST_${CMND_IDX}\n")

    # Print the message for this TEST_<IDX> block if set
    IF (TEST_${CMND_IDX}_MESSAGE)
      MESSAGE("${TEST_${CMND_IDX}_MESSAGE}\n")
    ENDIF()

    # Run the TEST_<IDX> block (Copy files or CMND)
    IF (TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR)
      SETUP_AND_RUN_TEST_IDX_COPY_FILES_BLOCK()
    ELSE()
      SETUP_AND_RUN_TEST_IDX_CMND_BLOCK()
    ENDIF()

    # Print the load and/or timing info for TEST_<IDX> block

    IF (SHOW_MACHINE_LOAD   AND  NOT  SHOW_COMMANDS_ONLY)
      PRINT_UPTIME("TEST_${CMND_IDX}: Uptime:")
    ENDIF()

    IF (SHOW_START_END_DATE_TIME  AND  NOT  SHOW_COMMANDS_ONLY)
      TIMER_GET_RAW_SECONDS(TEST_CMND_END)
      IF (TEST_CMND_START AND TEST_CMND_END)
        TIMER_PRINT_REL_TIME(${TEST_CMND_START} ${TEST_CMND_END}
           "TEST_${CMND_IDX}: Time")
      ELSE()
        MESSAGE("ERROR: Not able to return test times! Is 'date' in your path?")
      ENDIF()
      SET(TEST_CMND_START ${TEST_CMND_END})
    ENDIF()

    # Determine pass/fail for the TEST_<IDX> CMND block

    IF (NOT SHOW_COMMANDS_ONLY)

      # Determine pass/fail for TEST_<IDX> copy files or CMND

      IF (TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR)
        # Pass/fail already determined by setting TEST_CASE_PASSED in
        # SETUP_AND_RUN_TEST_IDX_COPY_FILES_BLOCK above.
      ELSE()
        DETERMINE_TEST_IDX_CMND_BLOCK_PASS_FAIL()
      ENDIF()

      # Print final pass/fail for the TEST_<IDX> block

      IF (TEST_CASE_PASSED)
        MESSAGE("TEST_${CMND_IDX}: Result = PASSED")
      ELSE()
        MESSAGE("TEST_${CMND_IDX}: Result = FAILED")
        SET(OVERALL_TEST_PASSED FALSE)
        IF (FAIL_FAST)
          MESSAGE("TEST_${CMND_IDX}: FAIL FAST, SKIPPING REST OF TEST CASES!")
          BREAK()
        ENDIF()
      ENDIF()

    ENDIF(NOT SHOW_COMMANDS_ONLY)

  ENDFOREACH()

  MESSAGE("\n${TEST_SEP}\n")

  #
  # C) Print the final test data and pass/fail
  #

  IF (NOT SHOW_COMMANDS_ONLY)

    IF (SHOW_START_END_DATE_TIME)
      PRINT_CURRENT_DATE_TIME("Ending at:")
      IF (TEST_OVERALL_START AND TEST_CMND_END)
        TIMER_PRINT_REL_TIME(${TEST_OVERALL_START} ${TEST_CMND_END}
          "OVERALL TEST TIME")
      ELSE()
        MESSAGE("ERROR: Not able to return test times! Is 'date' in your path?")
      ENDIF()
      MESSAGE("")
    ENDIF()

    IF (OVERALL_TEST_PASSED)
      MESSAGE("OVERALL FINAL RESULT: TEST PASSED (${TEST_NAME})")
    ELSE()
      MESSAGE("OVERALL FINAL RESULT: TEST FAILED (${TEST_NAME})")
    ENDIF()
  ELSE()
    MESSAGE("OVERALL FINAL RESULT: DID NOT RUN COMMANDS (${TEST_NAME})")
  ENDIF()

  MESSAGE("\n${ADVANDED_TEST_SEP}\n")

ENDFUNCTION()
