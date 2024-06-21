# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

# @FUNCTION: tribits_print_list()
#
# Unconditionally the name and values of a list var.
#
# Usage::
#
#   tribits_print_list(<listName>)
#
# This prints::
#
#   -- <listName> (size=<len>):
#   --   <val0>
#   --   <val1>
#   ...
#   --   <valN>
#
# The variable ``<listName>`` can be defined or undefined or empty.  This uses
# an explicit "-- " line prefix so that it prints nice even on Windows CMake.
#
function(tribits_print_list listName)
  list(LENGTH ${listName} len)
  message("-- " "${listName} (size=${len})")
  foreach(ele IN LISTS ${listName})
    message("-- " "  '${ele}'")
  endforeach()
endfunction()

# NOTE: Above, I was not able to call message_wrapper() directly because it
# was removing the ';' in array arguments.  This broke a bunch of unit tests.
# Therefore, I have to duplicate code and call it in two separate places.  I
# have to admit that CMake behavior surprises me many times.  This is not a
# great programming language.
