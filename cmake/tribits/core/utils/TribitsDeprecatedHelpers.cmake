# @HEADER
# ************************************************************************
#
# TriBITS: Tribal Build, Integrate, and Test System
# Copyright 2013 Sandia Corporation
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
