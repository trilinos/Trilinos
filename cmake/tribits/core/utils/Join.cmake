# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: join()
#
# Join a set of strings into a single string using a join string.
#
# Usage::
#
#   join(<outputStrVar> "<sepStr>" <quoteElements>
#     "<string0>" "<string1>" ...)
#
# Arguments:
#
#   ``<outputStrVar>``
#
#     The name of a variable that will hold the output string.
#
#   ``"<sepStr>"``
#
#     A string to use to join the list of strings.
#
#   ``<quoteElements>``
#
#     If ``TRUE``, then each ``<stringi>`` is quoted using an escaped quote
#      char ``\"``.  If ``FALSE`` then no escaped quote is used.
#
#   ``"<string0>" "<string1>" ...``
#
#     Zero or more string arguments to be joined.
#
# On output, the variable ``<outputStrVar>`` is set to::
#
#   "<string0><sepStr><string1><sepStr>..."
#
# If ``<quoteElements>=TRUE``, then ``<outputStrVar>`` is set to::
#
#   "\"<string0>\"<sepStr>\"<string1>\"<sepStr>..."
#
# For example, the latter can be used to set up a set of command-line
# arguments given a CMake array like::
#
#   join(CMND_LINE_ARGS " " TRUE ${CMND_LINE_ARRAY})
#
# WARNING: Be careful to quote string arguments that have spaces because CMake
# interprets those as array boundaries.
#
function(join  OUTPUT_STRING_VAR  SEP_STR  QUOTE_ELEMENTS)
  set(QUOTE_CHAR)
  if (QUOTE_ELEMENTS)
    set(QUOTE_CHAR "\"")
  endif()
  set(OUTPUT_STRING "")
  foreach(STRING_VAL ${ARGN})
    if (OUTPUT_STRING STREQUAL "")
      set(OUTPUT_STRING "${QUOTE_CHAR}${STRING_VAL}${QUOTE_CHAR}")
    else()
      set(OUTPUT_STRING "${OUTPUT_STRING}${SEP_STR}${QUOTE_CHAR}${STRING_VAL}${QUOTE_CHAR}")
    endif()
  endforeach()
  set(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
endfunction()
