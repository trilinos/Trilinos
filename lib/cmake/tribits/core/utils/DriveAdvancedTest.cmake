# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/PrintVar.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/AppendStringVar.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/Join.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/TimingUtils.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/TribitsGetCategoriesString.cmake")


function(print_current_date_time  PREFIX_STR)
  execute_process( COMMAND  date  OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  DATE_TIME )
  message("${PREFIX_STR} ${DATE_TIME}\n")
endfunction()


function(print_uptime  PREFIX_STR)
  execute_process( COMMAND  uptime  OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE  MACHINE_LOAD )
  message("${PREFIX_STR} ${MACHINE_LOAD}")
endfunction()


function(print_single_check_result  msgBegin  TEST_CASE_PASSED_IN)
  if (TEST_CASE_PASSED_IN)
    message("${msgBegin} [PASSED]")
  else()
    message("${msgBegin} [FAILED]")
  endif()
endfunction()


function(print_any_regex_pass_match  msgBegin  regexMatch)
  if (regexMatch)
    message("${msgBegin} [PASSED]")
  else()
    message("${msgBegin} [not matched]")
  endif()
endfunction()


function(print_any_regex_fail_match  msgBegin  regexMatch)
  if (regexMatch)
    message("${msgBegin} [FAILED]")
  else()
    message("${msgBegin} [does not match]")
  endif()
endfunction()


function(delete_create_working_directory  WORKING_DIR_IN   SKIP_CLEAN)
  if (EXISTS "${WORKING_DIR_IN}" AND NOT SKIP_CLEAN)
    message("Removing existing working directory"
      " '${WORKING_DIR_IN}'\n")
    if (NOT SHOW_COMMANDS_ONLY)
      file(REMOVE_RECURSE "${WORKING_DIR_IN}")
    endif()
  endif()
  if (NOT EXISTS "${WORKING_DIR_IN}")
    message("Creating new working directory"
      " '${WORKING_DIR_IN}'\n")
    if (NOT SHOW_COMMANDS_ONLY)
      file(MAKE_DIRECTORY "${WORKING_DIR_IN}")
    endif()
  endif()
endfunction()


macro(setup_and_run_test_idx_copy_files_block)

  message(
    "Copy files from:\n"
    "    ${TEST_${CMND_IDX}_SOURCE_DIR}/\n"
    "  to:\n"
    "    ${TEST_${CMND_IDX}_DEST_DIR}/\n")

  split("${TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR}"
    "," COPY_FILES_TO_TEST_DIR_ARRAY)

  set(TEST_CASE_PASSED TRUE)

  message("${OUTPUT_SEP}\n")

  # Make the dest directory if not exists (cmake -E copy will not create it)
  if (NOT EXISTS "${TEST_${CMND_IDX}_DEST_DIR}")
    message("Creating dest directory ${TEST_${CMND_IDX}_DEST_DIR}/ ...")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_${CMND_IDX}_DEST_DIR}"
      RESULT_VARIABLE MKDIR_COMMAND_RTN
      )
    if (NOT MKDIR_COMMAND_RTN EQUAL 0)
      set(TEST_CASE_PASSED FALSE)
    endif()
    # NOTE: Above, it would have been great to be able use file(MAKE_DIRECTORY
    # ...)  but if that fails it will just abortS the cmake -P script.  There
    # is no way to handle that gracefully.  Therefore, we have to use
    # execute_process(cmake -E make_directory ...).
  endif()

  # Copy the full list of files to the dest dir
  foreach(FILENAME ${COPY_FILES_TO_TEST_DIR_ARRAY})
    message("Copy file ${FILENAME} ...") 
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy
        "${TEST_${CMND_IDX}_SOURCE_DIR}/${FILENAME}"
        "${TEST_${CMND_IDX}_DEST_DIR}/"
      RESULT_VARIABLE COPY_COMMAND_RTN
      )
    if (NOT COPY_COMMAND_RTN EQUAL 0)
      set(TEST_CASE_PASSED FALSE)
   endif()
  endforeach()

  # NOTE: Above, it would have been great to be able use
  # configure_file(... COPYONLY) but if that fails it will just abortS the
  # cmake -P script.  There is no way to handle that gracefully.  Therefore,
  # we have to use execute_process(cmake -E copy ...).  Also, it would have
  # been great to just copy the files all at once with `cmake -E copy
  # <srcDir>/<file0> <srDir>/<file1> ... <distDir>/` but older versions of
  # CMake don't support that mode.  They only support copying the files one at
  # a time :-(

  message("\n${OUTPUT_SEP}\n")

  message("") # Need a vertical space for readability of the output

endmacro()


macro(setup_and_run_test_idx_cmnd_block)

  # Address working directory for this TEST_<IDX> block if set
  if (TEST_${CMND_IDX}_WORKING_DIRECTORY)
    if (NOT  IS_ABSOLUTE  "${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
      set(TEST_${CMND_IDX}_WORKING_DIRECTORY
        ${BASE_WORKING_DIRECTORY}/${TEST_${CMND_IDX}_WORKING_DIRECTORY})
    endif()
    delete_create_working_directory("${TEST_${CMND_IDX}_WORKING_DIRECTORY}"
      ${TEST_${CMND_IDX}_SKIP_CLEAN_WORKING_DIRECTORY})
  endif()

  # Set up the TEST_<IDX> command block
  join( TEST_CMND_STR " " TRUE ${TEST_${CMND_IDX}_CMND} )
  string(REPLACE "${LIST_SEPARATOR}" ";" TEST_CMND_STR "${TEST_CMND_STR}")
  message("Running: ${TEST_CMND_STR}\n")
  set(EXEC_CMND COMMAND ${TEST_${CMND_IDX}_CMND})
  string(REPLACE "${LIST_SEPARATOR}" "\\\;" EXEC_CMND "${EXEC_CMND}")

  # Set up the working directory that this TEST_<IDX> CMND block will run in

  set(WORKING_DIR_SET)
  if (TEST_${CMND_IDX}_WORKING_DIRECTORY)
    set(WORKING_DIR_SET "${TEST_${CMND_IDX}_WORKING_DIRECTORY}")
  elseif(OVERALL_WORKING_DIRECTORY)
    set(WORKING_DIR_SET "${OVERALL_WORKING_DIRECTORY}")
  endif()

  if (WORKING_DIR_SET)
    message("  Running in working directory \"${WORKING_DIR_SET}\"\n")
    set(WORKING_DIR "${WORKING_DIR_SET}")
  else()
    set(WORKING_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  # Set the actual command that will be run with execute_proces()

  set(EXEC_CMND ${EXEC_CMND}
    WORKING_DIRECTORY "${WORKING_DIR}"
    )

  # Set up the optional output file that the execute_process() command will write to

  if (TEST_${CMND_IDX}_OUTPUT_FILE)
    if (NOT  IS_ABSOLUTE  "${TEST_${CMND_IDX}_OUTPUT_FILE}")
      set(OUTPUT_FILE_USED "${WORKING_DIR}/${TEST_${CMND_IDX}_OUTPUT_FILE}")
    else()
      set(OUTPUT_FILE_USED "${TEST_${CMND_IDX}_OUTPUT_FILE}")
    endif()
    message("  Writing output to file \"${OUTPUT_FILE_USED}\"\n")
  endif()

  # Run the actual command with executte_process() (or just print what would run) ...

  if (NOT SHOW_COMMANDS_ONLY)

    # Provide the test configuration in an environment variable.
    if(TEST_CONFIG)
      set(ENV{TEST_CONFIG} "${TEST_CONFIG}")
    endif(TEST_CONFIG)

    execute_process(
      ${EXEC_CMND}
      OUTPUT_VARIABLE TEST_CMND_OUT
      ERROR_VARIABLE TEST_CMND_OUT
      RESULT_VARIABLE EXEC_RESULT
      )

    if (TEST_${CMND_IDX}_OUTPUT_FILE)
      file(WRITE "${OUTPUT_FILE_USED}" "${TEST_CMND_OUT}")
    endif()

    message("${OUTPUT_SEP}\n")

    if (NOT TEST_${CMND_IDX}_NO_ECHO_OUTPUT)
      message("${TEST_CMND_OUT}")
    else()
      message("NO_ECHO_OUTPUT\n")
    endif()

  else()

    message("\n*** Not running command on request ***")

  endif()

  message("${OUTPUT_SEP}\n")

endmacro()


macro(determine_test_idx_cmnd_block_pass_fail)

  message("TEST_${CMND_IDX}: Return code = ${EXEC_RESULT}")

  # A) Apply first set of pass/fail logic
  set(TEST_CASE_PASSED FALSE)
  if (TEST_${CMND_IDX}_PASS_ANY)
    set(TEST_CASE_PASSED TRUE)
    print_single_check_result(
      "TEST_${CMND_IDX}: Pass criteria = Pass Any"
      ${TEST_CASE_PASSED} )
  elseif (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION)
    set(TEST_CASE_PASSED FALSE)
    foreach(REGEX_STR ${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION})
      string(REGEX MATCH "${REGEX_STR}" MATCH_STR "${TEST_CMND_OUT}")
      if (MATCH_STR)
        set(TEST_CASE_PASSED TRUE)
      endif()
      print_any_regex_pass_match(
        "TEST_${CMND_IDX}: Pass criteria = Match any REGEX {${REGEX_STR}}"
        "${MATCH_STR}")
    endforeach()
  elseif (TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL)
    set(TEST_CASE_PASSED TRUE)
    foreach(REGEX_STR ${TEST_${CMND_IDX}_PASS_REGULAR_EXPRESSION_ALL})
      string(REGEX MATCH "${REGEX_STR}" MATCH_STR "${TEST_CMND_OUT}" )
      if (NOT "${MATCH_STR}" STREQUAL "")
        set(THIS_REGEX_MATCHED  TRUE)
      else()
        set(THIS_REGEX_MATCHED  FALSE)
      endif()
      if (NOT  THIS_REGEX_MATCHED)
        set(TEST_CASE_PASSED FALSE)
      endif()
      print_single_check_result(
        "TEST_${CMND_IDX}: Pass criteria = Match all REGEX {${REGEX_STR}}"
        ${THIS_REGEX_MATCHED} )
    endforeach()
  else()
    if (EXEC_RESULT EQUAL 0)
      set(TEST_CASE_PASSED TRUE)
    else()
      set(TEST_CASE_PASSED FALSE)
    endif()
    print_single_check_result(
      "TEST_${CMND_IDX}: Pass criteria = Zero return code"
      ${TEST_CASE_PASSED} )
  endif()

  # B) Check for failing regex matching?
  if (TEST_${CMND_IDX}_FAIL_REGULAR_EXPRESSION)
    foreach(REGEX_STR ${TEST_${CMND_IDX}_FAIL_REGULAR_EXPRESSION})
      string(REGEX MATCH "${REGEX_STR}" MATCH_STR "${TEST_CMND_OUT}")
      if (MATCH_STR)
        set(TEST_CASE_PASSED FALSE)
      endif()
      print_any_regex_fail_match(
        "TEST_${CMND_IDX}: Pass criteria = Not match REGEX {${REGEX_STR}}"
        "${MATCH_STR}")
    endforeach()
  endif()

  # C) Check for return code always 0?
  if (TEST_${CMND_IDX}_ALWAYS_FAIL_ON_NONZERO_RETURN)
    if (NOT EXEC_RESULT EQUAL 0)
      set(ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED FALSE)
      set(TEST_CASE_PASSED FALSE)
    else()
      set(ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED TRUE)
    endif()
    print_single_check_result(
      "TEST_${CMND_IDX}: Pass criteria = ALWAYS_FAIL_ON_NONZERO_RETURN"
      ${ALWAYS_FAIL_ON_NONZERO_RETURN_PASSED} )
  elseif (TEST_${CMND_IDX}_ALWAYS_FAIL_ON_ZERO_RETURN)
    if (EXEC_RESULT EQUAL 0)
      set(ALWAYS_FAIL_ON_ZERO_RETURN_PASSED FALSE)
      set(TEST_CASE_PASSED FALSE)
    else()
      set(ALWAYS_FAIL_ON_ZERO_RETURN_PASSED TRUE)
    endif()
    print_single_check_result(
      "TEST_${CMND_IDX}: Pass criteria = ALWAYS_FAIL_ON_ZERO_RETURN"
      ${ALWAYS_FAIL_ON_ZERO_RETURN_PASSED} )
  endif()

  # D) Invert pass/fail result?
  if (TEST_${CMND_IDX}_WILL_FAIL)
    if (TEST_CASE_PASSED)
      set(TEST_CASE_PASSED FALSE)
    else()
      set(TEST_CASE_PASSED TRUE)
    endif()
    print_single_check_result(
      "TEST_${CMND_IDX}: Pass criteria = WILL_FAIL (invert the above 'Pass critera')"
      ${TEST_CASE_PASSED} )
  endif()

endmacro()


function(drive_advanced_test)

  #
  # A) Print the header for the advanced test
  #

  set(ADVANDED_TEST_SEP
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

  set(TEST_SEP
    "================================================================================")

  set(OUTPUT_SEP
    "--------------------------------------------------------------------------------")

  math(EXPR LAST_CMND_IDX ${NUM_CMNDS}-1)

  message("\n${ADVANDED_TEST_SEP}\n")
  message("Advanced Test: ${TEST_NAME}\n")

  message("Selected Test/CTest Properties:")
  tribits_get_categories_string("${CATEGORIES}" CATEGORIES_IN_COMMAS)
  message("  CATEGORIES = ${CATEGORIES_IN_COMMAS}")
  message("  PROCESSORS = ${PROCESSORS}")
  if (TIMEOUT)
    message("  TIMEOUT    = ${TIMEOUT}\n")
  else()
    message("  TIMEOUT    = DEFAULT\n")
  endif()

  if (SHOW_MACHINE_LOAD  AND  NOT  SHOW_COMMANDS_ONLY)
    print_uptime("Starting Uptime:")
  endif()

  if (SHOW_START_END_DATE_TIME  AND  NOT  SHOW_COMMANDS_ONLY)
    print_current_date_time("Starting at:")
  endif()

  if (OVERALL_WORKING_DIRECTORY)
    if (NOT  IS_ABSOLUTE  "${OVERALL_WORKING_DIRECTORY}")
      set(OVERALL_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${OVERALL_WORKING_DIRECTORY})
    endif()
    delete_create_working_directory("${OVERALL_WORKING_DIRECTORY}"
      ${SKIP_CLEAN_OVERALL_WORKING_DIRECTORY})
    set(BASE_WORKING_DIRECTORY "${OVERALL_WORKING_DIRECTORY}")
  else()
    set(BASE_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  foreach ( CMND_IDX RANGE ${LAST_CMND_IDX} )
    if (CMND_IDX EQUAL 0)
      set(TEST_NAMES_STR "TEST_0")
    else()
      string(APPEND  TEST_NAMES_STR ", TEST_${CMND_IDX}" )
    endif()
  endforeach()
  message("Running test commands: ${TEST_NAMES_STR}")

  if (SHOW_START_END_DATE_TIME AND NOT SHOW_COMMANDS_ONLY)
   timer_get_raw_seconds(TEST_CMND_START)
   set(TEST_OVERALL_START ${TEST_CMND_START})
  endif()

  #
  # B) Loop over and run the TEST_<IDX> blocks one at a time
  #

  set(OVERALL_TEST_PASSED TRUE)

  foreach ( CMND_IDX RANGE ${LAST_CMND_IDX} )
    message("\n${TEST_SEP}\n")
    message("TEST_${CMND_IDX}\n")

    # Print the message for this TEST_<IDX> block if set
    if (TEST_${CMND_IDX}_MESSAGE)
      message("${TEST_${CMND_IDX}_MESSAGE}\n")
    endif()

    # Run the TEST_<IDX> block (Copy files or CMND)
    if (TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR)
      setup_and_run_test_idx_copy_files_block()
    else()
      setup_and_run_test_idx_cmnd_block()
    endif()

    # Print the load and/or timing info for TEST_<IDX> block

    if (SHOW_MACHINE_LOAD   AND  NOT  SHOW_COMMANDS_ONLY)
      print_uptime("TEST_${CMND_IDX}: Uptime:")
    endif()

    if (SHOW_START_END_DATE_TIME  AND  NOT  SHOW_COMMANDS_ONLY)
      timer_get_raw_seconds(TEST_CMND_END)
      if (TEST_CMND_START AND TEST_CMND_END)
        timer_print_rel_time(${TEST_CMND_START} ${TEST_CMND_END}
           "TEST_${CMND_IDX}: Time")
      else()
        message("ERROR: Not able to return test times! Is 'date' in your path?")
      endif()
      set(TEST_CMND_START ${TEST_CMND_END})
    endif()

    # Determine pass/fail for the TEST_<IDX> CMND block

    if (NOT SHOW_COMMANDS_ONLY)

      # Determine pass/fail for TEST_<IDX> copy files or CMND

      if (TEST_${CMND_IDX}_COPY_FILES_TO_TEST_DIR)
        # Pass/fail already determined by setting TEST_CASE_PASSED in
        # SETUP_AND_RUN_TEST_IDX_COPY_FILES_BLOCK above.
      else()
        determine_test_idx_cmnd_block_pass_fail()
      endif()

      # Print final pass/fail for the TEST_<IDX> block

      if (TEST_CASE_PASSED)
        message("TEST_${CMND_IDX}: Result = PASSED")
      else()
        message("TEST_${CMND_IDX}: Result = FAILED")
        set(OVERALL_TEST_PASSED FALSE)
        if (FAIL_FAST)
          message("TEST_${CMND_IDX}: FAIL FAST, SKIPPING REST OF TEST CASES!")
          break()
        endif()
      endif()

    endif(NOT SHOW_COMMANDS_ONLY)

  endforeach()

  message("\n${TEST_SEP}\n")

  #
  # C) Print the final test data and pass/fail
  #

  if (NOT SHOW_COMMANDS_ONLY)

    if (SHOW_START_END_DATE_TIME)
      print_current_date_time("Ending at:")
      if (TEST_OVERALL_START AND TEST_CMND_END)
        timer_print_rel_time(${TEST_OVERALL_START} ${TEST_CMND_END}
          "OVERALL TEST TIME")
      else()
        message("ERROR: Not able to return test times! Is 'date' in your path?")
      endif()
      message("")
    endif()

    if (OVERALL_TEST_PASSED)
      message("OVERALL FINAL RESULT: TEST PASSED (${TEST_NAME})")
    else()
      message("OVERALL FINAL RESULT: TEST FAILED (${TEST_NAME})")
    endif()
  else()
    message("OVERALL FINAL RESULT: DID NOT RUN COMMANDS (${TEST_NAME})")
  endif()

  message("\n${ADVANDED_TEST_SEP}\n")

endfunction()
