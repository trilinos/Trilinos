# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: append_set()
#
# Utility function to append elements to a variable (reduces boiler-plate
# code).
#
# Usage::
#
#   append_set(<varName> <arg0> <arg1> ...)
#
# This just calls::
#
#   list(APPEND <varName> <arg0> <arg1> ...)
#
# There is better error reporting if one misspells ``APPEND_SET`` than if one
# misspells ``APPEND``.
#
macro(append_set VARNAME)
  list(APPEND ${VARNAME} ${ARGN})
endmacro()
