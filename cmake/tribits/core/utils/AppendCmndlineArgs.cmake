# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: append_cmndline_args()
#
# Utility function that appends command-line arguments to a variable of
# command-line arguments.
#
# Usage::
#
#   append_cmndline_args(<var> "<extraArgs>")
#
# This function just appends the command-line arguments in the string
# ``"<extraArgs>"`` but does not add an extra space if ``<var>`` is empty on
# input.  This just makes the formatting of command-line arguments easier.
#
function(append_cmndline_args  CMNDLINE_VAR_NAME  EXTRAARGS)
  if (${CMNDLINE_VAR_NAME})
    set(${CMNDLINE_VAR_NAME} "${${CMNDLINE_VAR_NAME}} ${EXTRAARGS}" PARENT_SCOPE)
  else()
    set(${CMNDLINE_VAR_NAME} "${EXTRAARGS}" PARENT_SCOPE)
  endif()
#  print_var(${CMNDLINE_VAR_NAME})
endfunction()
