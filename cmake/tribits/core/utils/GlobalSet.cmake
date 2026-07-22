# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: global_set()
#
# Set a variable as an internal global (cache) variable (removes boiler-plate
# code).
#
# Usage::
#
#   global_set(<varName> [other args])
#
# This just calls::
#
#   set(<varName> [other args] CACHE INTERNAL "")
#
# This avoid misspelling ``CACHE``.
#
macro(global_set VARNAME)
  set(${VARNAME} ${ARGN} CACHE INTERNAL "")
endmacro()
