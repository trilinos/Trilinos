# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(PrependCmndlineArgs)
include(DualScopeSet)


# @MACRO: dual_scope_prepend_cmndline_args()
#
# Utility function that prepends command-line arguments to a variable of
# command-line arguments and sets the result in current scope and parent
# scope.
#
# Usage::
#
#   dual_scope_prepend_cmndline_args(<var> "<extraArgs>")
#
# Just calls `prepend_cmndline_args()`_ and then ``set(<var> ${<var>} PARENT_SCOPE)``.
#
macro(dual_scope_prepend_cmndline_args  CMNDLINE_VAR_NAME  EXTRAARGS)
  prepend_cmndline_args(${CMNDLINE_VAR_NAME} "${EXTRAARGS}")
  set(${CMNDLINE_VAR_NAME} "${${CMNDLINE_VAR_NAME}}" PARENT_SCOPE)
endmacro()
