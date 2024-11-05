# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AssertDefined)
include(PrintVar)


# @FUNCTION: print_nonempty_var()
#
# Print a defined variable giving its name then value only if it is not empty.
#
# Usage::
#
#    print_nonempty_var(<varName>)
#
# Calls ``print_var(<varName>)`` if ``${<varName>}`` is not empty.
#
function(print_nonempty_var VARIABLE_NAME)
  assert_defined(VARIABLE_NAME)
  if (NOT "${${VARIABLE_NAME}}" STREQUAL "")
    print_var(${VARIABLE_NAME})
  endif()
endfunction()
