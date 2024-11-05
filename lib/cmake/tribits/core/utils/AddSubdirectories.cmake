# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: add_subdirectories()
#
# Macro that adds a list of subdirectories all at once (removes boiler-plate
# code).
#
# Usage::
#
#   add_subdirectories(<dir1> <dir2> ...)
#
# instead of::
#
#   add_subdirectory(<dir1>)
#   add_subdirectory(<dir2>)
#   ...
#
macro(add_subdirectories)
  foreach(DIR ${ARGV})
    add_subdirectory(${DIR})
  endforeach()
endmacro()
