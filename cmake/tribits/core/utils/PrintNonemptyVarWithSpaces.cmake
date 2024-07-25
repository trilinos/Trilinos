# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AssertDefined)
include(AppendStringVarWithSep)


# @FUNCTION: print_nonempty_var_with_spaces()
#
# Print a list variable giving its name then value printed with spaces instead
# of ``';'``, but only if the list is non-empty.
#
# Usage::
#
#    print_nonempty_var_with_spaces(<varName>  <printedVarOut>)
#
# Prints the variable as::
#
#    <varName>: <ele0> <ele1> ...
#
# If ``<printedVarOut>`` is ``TRUE`` on input, then the variable is not
# touched. If however, the variable ``<printedVarOut>`` is not ``TRUE`` and
# the list ``<varName>`` in non-empty, then ``<printedVarOut>`` is set to
# ``TRUE`` on output.
#
function(print_nonempty_var_with_spaces  variableName  printedVarOut)
  if (NOT "${${variableName}}" STREQUAL "")
    string(REPLACE  ";"  " "  OUTSTR  "${${variableName}}")
    message("-- ${variableName}: ${OUTSTR}")
    set(${printedVarOut}  TRUE  PARENT_SCOPE)
  endif()
endfunction()
