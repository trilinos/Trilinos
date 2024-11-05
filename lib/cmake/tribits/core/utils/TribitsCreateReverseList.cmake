# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: tribits_create_reverse_list()
#
# Create a reverse list var in one shot.
#
# Usage::
#
#   tribits_create_reverse_list(<oldListName> <newListName>)
#
macro(tribits_create_reverse_list  oldListName  newListName)
  set(${newListName} ${${oldListName}})
  list(REVERSE ${newListName})
endmacro()
