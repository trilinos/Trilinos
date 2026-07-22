# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: prepend_cmndline_args()
#
# Utility function that prepends command-line arguments to a variable of
# command-line arguments.
#
# Usage::
#
#   prepend_cmndline_args(<var> "<extraArgs>")
#
# This function just prepends the command-line arguments in the string
# ``"<extraArgs>"`` but does not add an extra space if ``<var>`` is empty on
# input.
#
function(prepend_cmndline_args  CMNDLINE_VAR_NAME  EXTRAARGS)
  if (${CMNDLINE_VAR_NAME})
    set(${CMNDLINE_VAR_NAME} "${EXTRAARGS} ${${CMNDLINE_VAR_NAME}}" PARENT_SCOPE)
  else()
    set(${CMNDLINE_VAR_NAME} "${EXTRAARGS}" PARENT_SCOPE)
  endif()
  #print_var(${CMNDLINE_VAR_NAME})
endfunction()
