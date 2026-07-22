# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: prepend_set()
#
# Utility macro to prepend elements to a variable (reduces boiler-plate code).
#
# Usage::
#
#   prepend_set(<varName> <arg0> <arg1> ...)
#
# Just calls::
#
#   set(<varName> <arg0> <arg1> ... ${<varName>})
#
# NOTE: Prepending is not as efficient as appending so prefer `append_set()`_
# or just ``list(APPEND ...)``.
#
macro(prepend_set VARNAME)
  set(${VARNAME} ${ARGN} ${${VARNAME}})
endmacro()
