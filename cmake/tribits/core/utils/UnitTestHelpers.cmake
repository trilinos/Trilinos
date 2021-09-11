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

include(CMakeParseArguments)
include(GlobalSet)


#
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


#
# @FUNCTION: unittest_string_regex()
#
# Perform a series regexes of given strings and update overall test statistics.
#
# Usage::
#
#   unittest_string_regex(
#     <inputString>
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

    if (REGEX_MATCH_RESULT)
      message("  Searching for REGEX {${REGEX}}:  [PASSED]\n")
      math( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
      global_set(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
    else()
      message("  Searching for REGEX {${REGEX}}:  [FAILED]\n")
      global_set(UNITTEST_OVERALL_PASS FALSE)
      message(WARNING "Stack trace for failed unit test")
    endif()

  endforeach()

endfunction()


#
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


#
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


#
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


#
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

