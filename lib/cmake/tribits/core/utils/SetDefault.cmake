# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: set_default()
#
# Give a local variable a default value if a non-empty value is not already
# set.
#
# Usage::
#
#   set_default(<varName> <arg0> <arg1> ...)
#
# If on input ``"${<varName>}"==""``, then ``<varName>`` is set to the given
# default ``<arg0> <arg1> ...``.  Otherwise, the existing non-empty value is
# preserved.
#
macro(set_default VAR)
  if ("${${VAR}}" STREQUAL "")
    set(${VAR} ${ARGN})
  endif()
endmacro()
