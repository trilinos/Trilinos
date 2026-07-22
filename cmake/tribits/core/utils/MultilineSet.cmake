# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: multiline_set()
#
# Function to set a single string by concatenating a list of separate strings
#
# Usage::
#
#   multiline_set(<outputStrVar>
#     "<string0>"
#     "<string1>"
#     ...
#     )
#
# On output, the local variables ``<outputStrVar>`` is set to::
#
#   "<string0><string1>..."
#
# The purpose of this is function to make it easier to set longer strings over
# multiple lines.
#
# This function is exactly the same as `concat_strings()`_ and should not even
# exist :-(
#
function(multiline_set VARAIBLE_NAME)

  set(MULTILINE_SET_LOCAL_STR "")

  foreach(LINE_STR ${ARGN})
    set(MULTILINE_SET_LOCAL_STR "${MULTILINE_SET_LOCAL_STR}${LINE_STR}")
  endforeach()

  set(${VARAIBLE_NAME} "${MULTILINE_SET_LOCAL_STR}" PARENT_SCOPE)

endfunction()
