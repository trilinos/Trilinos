# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/MessageWrapper.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/TribitsParseArgumentsHelpers.cmake")


set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_VALUES_THAT_CALL_MESSAGE
  DEPRECATION  AUTHOR_WARNING  SEND_ERROR  FATAL_ERROR )
set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_VALUES_THAT_DONT_CALL_MESSAGE
  IGNORE )
set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_ALL_VALID_VALUES
  ${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_VALUES_THAT_CALL_MESSAGE}
  ${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_VALUES_THAT_DONT_CALL_MESSAGE} )


# @FUNCTION: tribits_deprecated()
#
# Notify the user that some TriBITS functionality is deprecated.
#
# Usage::
#
#   tribits_deprecated(<message>)
#
# Depending on the value of the cache variable
# `TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE`_, this can do one of several
# things:
#
# - ``DEPRECATION`` or empty string (or variable not defined): Issue a CMake
#    ``DEPRECATION`` message and continue.
# - ``AUTHOR_WARNING``: Issue a CMake ``AUTHOR_WARNING`` message and continue.
# - ``SEND_ERROR``: Issue a CMake ``SEND_ERROR`` message and continue.
# - ``FATAL_ERROR``: Issue a CMake ``FATAL_ERROR`` message and exit.
# - ``IGNORE``: Issue no message and continue.
#
function(tribits_deprecated)
  cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")

  if ("${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE}" STREQUAL "")
    set(TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE DEPRECATION)
  endif()

  if ("${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE}"  IN_LIST  TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_VALUES_THAT_CALL_MESSAGE)
    message_wrapper("${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE}"
      ${FWD_UNPARSED_ARGUMENTS}
      "\n\nNOTE: To Make these warnings go away, set -D"
      " TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE=IGNORE (see the build reference guide).")
  elseif (NOT "${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE}"  IN_LIST
      TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE_ALL_VALID_VALUES
    )
    message_wrapper(FATAL_ERROR "Invalid value for"
      " TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE="
      "'${TRIBITS_HANDLE_TRIBITS_DEPRECATED_CODE}'")
  endif()
endfunction()


# @FUNCTION: tribits_deprecated_command()
#
# Notify the user that a TriBITS function or macro is deprecated. This should
# be the first command called at the top of any deprecated function or macro.
#
# Usage::
#
#   tribits_deprecated_command(<name>
#     [MESSAGE <message>]
#     )
#
function(tribits_deprecated_command  name)
  # Parse input arguments
  set(argMultiValArgKeywords  MESSAGE)
  cmake_parse_arguments(PARSE_ARGV  1  PREFIX
    ""   # options
    ""   # one_value_keywords
    "${argMultiValArgKeywords}"   # multi_value_keywords
    )
  tribits_check_for_unparsed_arguments(PREFIX)

  set(deprecationMessage "TriBITS command '${name}' is deprecated.")
  if (NOT "${PREFIX_MESSAGE}" STREQUAL "")
    string(APPEND deprecationMessage "\n\n${PREFIX_MESSAGE}")
  endif()

  tribits_deprecated("${deprecationMessage}")
endfunction()
