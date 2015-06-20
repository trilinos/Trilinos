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
# @FUNCTION: JOIN()
#
# Join a set of strings into a single string using a join string.
#
# Usage::
#
#   JOIN(<outputStrVar> "<sepStr>" <quoteElements>
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
#   JOIN(CMND_LINE_ARGS " " TRUE ${CMND_LINE_ARRAY})
#
# WARNING: Be careful to quote string arguments that have spaces because CMake
# interprets those as array boundaries.
#
FUNCTION(JOIN  OUTPUT_STRING_VAR  SEP_STR  QUOTE_ELEMENTS)
  SET(QUOTE_CHAR)
  IF (QUOTE_ELEMENTS)
    SET(QUOTE_CHAR "\"")
  ENDIF()
  SET(OUTPUT_STRING "")
  FOREACH(STRING_VAL ${ARGN})
    IF (OUTPUT_STRING STREQUAL "")
      SET(OUTPUT_STRING "${QUOTE_CHAR}${STRING_VAL}${QUOTE_CHAR}")
    ELSE()
      SET(OUTPUT_STRING "${OUTPUT_STRING}${SEP_STR}${QUOTE_CHAR}${STRING_VAL}${QUOTE_CHAR}")
    ENDIF()
  ENDFOREACH()
  SET(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
ENDFUNCTION()
