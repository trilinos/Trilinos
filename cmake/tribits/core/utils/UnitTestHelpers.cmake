# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(CMakeParseArguments)
include(GlobalSet)

# Set policy here instead of including TribitCTestPolicis.cmake since we want
# this to be a standalone module
cmake_policy(SET CMP0007 NEW)  # Don't ignore empty list items


# @MACRO: unittest_initialize_vars()
#
# Call to initialize the unit test variables before running unit tests.
#
# Usage::
#
#   unittest_initialize_vars()
#
macro(unittest_initialize_vars)
  # Assume that all unit tests will pass by default
  global_set(UNITTEST_OVERALL_PASS TRUE)
  global_set(UNITTEST_OVERALL_NUMPASSED 0)
  global_set(UNITTEST_OVERALL_NUMRUN 0)
endmacro()


# @FUNCTION: unittest_compare_const()
#
# Perform a single unit test equality check and update overall test statistics
#
# Usage::
#
#   unittest_compare_const(<varName> <expectedValue>)
#
# If ``${<varName>} == <expectedValue>``, then the check passes, otherwise it
# fails.  This prints the variable name and values and shows the test result.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_compare_const VAR_NAME CONST_VAL)

  math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  message(
    "\nCheck:\n"
    "    ${VAR_NAME} =\n"
    "    [${${VAR_NAME}}]\n"
    "  EQUALS:\n"
    "    [${CONST_VAL}]"
    )

  if ("${${VAR_NAME}}" STREQUAL "${CONST_VAL}")
    message("  [PASSED]\n")
    math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  else()
    message("  [FAILED]\n")
    global_set(UNITTEST_OVERALL_PASS FALSE)
    message(WARNING "Stack trace for failed unit test")
  endif()

endfunction()


# @FUNCTION: unittest_compare_list_ele_const()
#
# Perform a single unit test equality check for a single list element
#
# Usage::
#
#   unittest_compare_list_ele_const(<listName> <idx> <expectedConstValue>)
#
# If ``${<listName>[<idx>]} == <expectedValue>``, then the check passes, otherwise it
# fails.  This prints the variable name and values and shows the test result.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_compare_list_ele_const  listName  idx  expectedConstValue)

  math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  list(GET "${listName}" ${idx} listEleIdx)
  set(listNameAndIdx "${listName}[${idx}]")

  message(
    "\nCheck:\n"
    "    ${listNameAndIdx} =\n"
    "    [${listEleIdx}]\n"
    "  EQUALS:\n"
    "    [${expectedConstValue}]"
    )

  if ("${listEleIdx}" STREQUAL "${expectedConstValue}")
    message("  [PASSED]\n")
    math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  else()
    message("  [FAILED]\n")
    global_set(UNITTEST_OVERALL_PASS FALSE)
    message(WARNING "Stack trace for failed unit test")
  endif()

endfunction()


# @FUNCTION: unittest_string_block_compare()
#
# Compare two string blocks (with multiple newlines '\n') line-by-line and
# print the first line that fails the comparison.
#
# Usage::
#
#   unittest_string_block_compare(
#     <stringVar> "<stringExpected>"
#     )
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_string_block_compare  stringVar  stringExpected)

  message("\nCheck: ${stringVar} equals expected string:")

  math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  string(REPLACE "\n" ";" stringList "${${stringVar}}")
  string(REPLACE "\n" ";" stringExpectedList "${stringExpected}")

  list(LENGTH  stringList  stringLen)
  list(LENGTH  stringExpectedList  stringExpectedLen)

  # minLen = min(stringLen, stringExpectedLen)
  set(minLen ${stringLen})
  if (stringExpectedLen LESS minLen)
    set(minLen ${stringExpectedLen})
  endif()

  set(allMatched TRUE)
  set(idx 0)
  while (idx LESS minLen)
    list(GET stringList ${idx} stringEle)
    list(GET stringExpectedList ${idx} stringExpectedEle)
    if (NOT stringEle STREQUAL stringExpectedEle)
      message(
        "  Error: Line ${idx} of ${stringVar}:\n"
        "    '${stringEle}'\n"
        "    !=\n"
        "    '${stringExpectedEle}'\n"
        "  [FAILED]\n"
        )
      global_set(UNITTEST_OVERALL_PASS FALSE)
      message(WARNING "Stack trace for failed unit test")
      set(allMatched FALSED)
      break()
    endif()
    math(EXPR idx "${idx}+1")
  endwhile()

  if (NOT allMatched)
    # The error handling was already handled above
  elseif (NOT stringLen EQUAL stringExpectedLen)
    # The lines of the strings matched but one of the strings had more lines
    # than the other
    if (stringLen GREATER stringExpectedLen)
      list(GET stringList ${stringExpectedLen} nextLine)
      message(
        "  Error: ${stringVar} has ${stringLen} lines where expected string has"
        " only ${stringExpectedLen} lines and the next extra line is:\n"
        "    '${nextLine}'\n"
        "  [FAILED]\n"
        )
    elseif (stringExpectedLen GREATER stringLen)
      list(GET stringExpectedList ${stringLen} nextLine)
      message(
        "  Error: Expected string has ${stringExpectedLen} lines where ${stringVar} has"
        "only ${stringLen} lines and the next extra line is:\n"
        "    '${nextLine}'\n"
        "  [FAILED]\n"
        )
    endif()
    global_set(UNITTEST_OVERALL_PASS FALSE)
    message(WARNING "Stack trace for failed unit test")
  else()
    # The strings matched exactly!
    message("  [PASSED]\n")
    math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  endif()

endfunction()


# @FUNCTION: unittest_string_regex()
#
# Perform a series of regexes on a given string and update overall test
# statistics.
#
# Usage::
#
#   unittest_string_regex(
#     "<inputString>"
#     REGEX_STRINGS "<str0>" "<str1>" ...
#     )
#
# If the ``<inputString>`` matches all of the of the regexs ``"<str0>"``,
# ``"<str1>"``, ..., then the test passes.  Otherwise it fails.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_string_regex INPUT_STRING)

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     ""
     #one_value_keywords
     ""
     #multi_value_keywords
     "REGEX_STRINGS"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  foreach(REGEX ${PARSE_REGEX_STRINGS})

    math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
    global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

    string(REGEX MATCH "${REGEX}" REGEX_MATCH_RESULT "${INPUT_STRING}")

    message("  Searching string:")
    message("     '${INPUT_STRING}'")
    if (REGEX_MATCH_RESULT)
      message("  for REGEX {${REGEX}}: [PASSED]\n")
      math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
      global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
    else()
      message("  for REGEX {${REGEX}}: [FAILED]\n")
      global_set(UNITTEST_OVERALL_PASS FALSE)
      message(WARNING "Stack trace for failed unit test")
    endif()

  endforeach()

endfunction()


# @FUNCTION: unittest_string_var_regex()
#
# Perform a series of regexes on a given string variable and update overall
# test statistics.
#
# Usage::
#
#   unittest_string_var_regex(
#     <inputStringVar>
#     REGEX_STRINGS "<str0>" "<str1>" ...
#     )
#
# If the ``"${<inputStringVar>}"`` matches all of the of the regexs
# ``"<str0>"``, ``"<str1>"``, ..., then the test passes.  Otherwise it fails.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_string_var_regex  inputStringVar)

  cmake_parse_arguments(PARSE_ARGV 1
     PARSE "" "" # prefix, options, one_value_keywords
     "REGEX_STRINGS" #multi_value_keywords
     )
  tribits_check_for_unparsed_arguments(PARSE)

  foreach(REGEX ${PARSE_REGEX_STRINGS})

    math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
    global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

    string(REGEX MATCH "${REGEX}" REGEX_MATCH_RESULT "${${inputStringVar}}")

    message("Searching string variable value '${inputStringVar}':")
    message("     '${${inputStringVar}}'")
    if (REGEX_MATCH_RESULT)
      message("  for REGEX {${REGEX}}: [PASSED]\n")
      math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
      global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
    else()
      message("  for REGEX {${REGEX}}: [FAILED]\n")
      global_set(UNITTEST_OVERALL_PASS FALSE)
      message(WARNING "Stack trace for failed unit test")
    endif()

  endforeach()

endfunction()


# @FUNCTION: unittest_has_substr_const()
#
# Check that a given string var contains the given substring and update
# overall test statistics
#
# Usage::
#
#   unittest_has_substr_const(<varName> <substr>)
#
# If ``${<varName>}`` contains the substring ``<substr>``, then the check
# passes, otherwise it fails.  This prints the variable name and values and
# shows the test result.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_has_substr_const VAR_NAME SUBSTR_VAL)

  math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  message(
    "\nCheck:\n"
    "    ${VAR_NAME} =\n"
    "    [${${VAR_NAME}}]\n"
    "  Contains:\n"
    "    [${SUBSTR_VAL}]"
    )

  string(FIND "${${VAR_NAME}}" "${SUBSTR_VAL}" SUBSTR_START_IDX)
  #print_var(SUBSTR_START_IDX)

  if (${SUBSTR_START_IDX} GREATER -1)
    message("  [PASSED]\n")
    math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  else()
    message("  [FAILED]\n")
    global_set(UNITTEST_OVERALL_PASS FALSE)
    message(WARNING "Stack trace for failed unit test")
  endif()

endfunction()


# @FUNCTION: unittest_not_has_substr_const()
#
# Check that a given string var does **NOT** contains the given substring and
# update overall test statistics
#
# Usage::
#
#   unittest_not_has_substr_const(<varName> <substr>)
#
# If ``${<varName>}`` contains the substring ``<substr>``, then the check
# failed, otherwise it passes.  This prints the variable name and values and
# shows the test result.
#
# This updates the global variables ``UNITTEST_OVERALL_NUMRUN``,
# ``UNITTEST_OVERALL_NUMPASSED``, and ``UNITTEST_OVERALL_PASS`` which are used
# by the unit test harness system to assess overall pass/fail.
#
function(unittest_not_has_substr_const VAR_NAME SUBSTR_VAL)

  math( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  global_set(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  message(
    "\nCheck:\n"
    "    ${VAR_NAME} =\n"
    "    [${${VAR_NAME}}]\n"
    "  Does NOT contain:\n"
    "    [${SUBSTR_VAL}]"
    )

  string(FIND "${${VAR_NAME}}" "${SUBSTR_VAL}" SUBSTR_START_IDX)
  #print_var(SUBSTR_START_IDX)

  if (${SUBSTR_START_IDX} GREATER -1)
    message("  [FAILED]\n")
    global_set(UNITTEST_OVERALL_PASS FALSE)
    message(WARNING "Stack trace for failed unit test")
  else()
    message("  [PASSED]\n")
    math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  endif()

endfunction()


# @FUNCTION: unittest_file_regex()
#
# Perform a series regexes of given strings and update overall test statistics.
#
# Usage::
#
#   unittest_file_regex(
#     <inputFileName>
#     REGEX_STRINGS "<str1>" "<str2>" ...
#     )
#
# The contents of ``<inputFileName>`` are read into a string and then passed
# to `unittest_string_regex()`_ to assess pass/fail.
#
function(unittest_file_regex  INPUT_FILE)
  message("\nRegexing for strings in the file '${INPUT_FILE}':\n")
  file(READ "${INPUT_FILE}" INPUT_FILE_STRING)
  unittest_string_regex("${INPUT_FILE_STRING}" ${ARGN})
endfunction()


# @FUNCTION: unittest_final_result()
#
# Print final statistics from all tests and assert final pass/fail
#
# Usage::
#
#   unittest_final_result(<expectedNumPassed>)
#
# If ``${UNITTEST_OVERALL_PASS}==TRUE`` and ``${UNITTEST_OVERALL_NUMPASSED} ==
# <expectedNumPassed>``, then the overall test program is determined to have
# passed and string::
#
#  "Final UnitTests Result: PASSED"
#
# is printed.  Otherwise, the overall test program is determined to have
# failed, the string::
#
#  "Final UnitTests Result: FAILED"
#
# is printed, and ``message(SEND_ERROR "FAIL")`` is called.
#
# The reason that we require passing in the expected number of passed tests is
# as an extra precaution to make sure that important unit tests are not left
# out.  CMake is a very loosely typed language and it pays to be a little
# paranoid.
#
function(unittest_final_result  EXPECTED_NUMPASSED)
   message("\nFinal UnitTests Result: num_run = ${UNITTEST_OVERALL_NUMRUN}\n")
  if (UNITTEST_OVERALL_PASS)
    if (UNITTEST_OVERALL_NUMPASSED EQUAL ${EXPECTED_NUMPASSED})
      message("Final UnitTests Result: PASSED"
        " (num_passed = ${UNITTEST_OVERALL_NUMPASSED})")
    else()
      message("\nError: num_passed = ${UNITTEST_OVERALL_NUMPASSED}"
        " != num_expected = ${EXPECTED_NUMPASSED}")
      message("\nFinal UnitTests Result: FAILED\n")
      message(SEND_ERROR "FAIL")
    endif()
  else()
    message("\nFinal UnitTests Result: FAILED\n")
  endif()
endfunction()

