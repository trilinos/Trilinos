# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

# @FUNCTION: print_var()
#
# Unconditionally print a variable giving its name then value.
#
# Usage::
#
#   print_var(<varName>)
#
# This prints::
#
#   message("-- " "${VARIBLE_NAME}='${${VARIBLE_NAME}}'")
#
# The variable ``<varName>`` can be defined or undefined or empty.  This uses
# an explicit "-- " line prefix so that it prints nice even on Windows CMake.
#
function(print_var VARIBLE_NAME)
  if (MESSAGE_WRAPPER_UNIT_TEST_MODE)
    global_set(MESSAGE_WRAPPER_INPUT "${MESSAGE_WRAPPER_INPUT}"
      "-- " "${VARIBLE_NAME}='${${VARIBLE_NAME}}'" "${PRINT_MSG}")
  else()
    message("-- " "${VARIBLE_NAME}='${${VARIBLE_NAME}}'")
  endif()
endfunction()

# NOTE: Above, I was not able to call message_wrapper() directly because it
# was removing the ';' in array arguments.  This broke a bunch of unit tests.
# Therefore, I have to duplicate code and call it in two separate places.  I
# have to admit that CMake behavior surprises me many times.  This is not a
# great programming language.
