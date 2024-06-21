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


# @FUNCTION: print_var_with_spaces()
#
# Print a defined variable giving its name then value printed with spaces
# instead of ``';'``.
#
# Usage::
#
#    print_var_with_spaces(<varName>  <printedVarInOut>)
#
# Prints the variable as::
#
#    <varName>: <ele0> <ele1> ...
#
# If ``$<printedVarInOut>`` is ``TRUE`` on input, then the variable is not
# touched. If however, the variable ``$<printedVarInOut>`` is not ``TRUE`` on
# input, then it is set to ``TRUE`` on output.
#
function(print_var_with_spaces  VARIBLE_NAME  PRINTED_VAR_OUT)
  assert_defined(VARIBLE_NAME)
  string(REPLACE ";" " " OUTSTR "${${VARIBLE_NAME}}")
  message("-- ${VARIBLE_NAME}: ${OUTSTR}")
  set(${PRINTED_VAR_OUT} TRUE PARENT_SCOPE)
endfunction()
