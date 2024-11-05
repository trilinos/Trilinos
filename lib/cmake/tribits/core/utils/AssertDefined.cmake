# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: assert_defined()
#
# Assert that a variable is defined and if not call ``message(SEND_ERROR
# ...)``.
#
# Usage::
#
#   assert_defined(<varName>)
#
# This is used to get around the problem of CMake not asserting the
# dereferencing of undefined variables.  For example, how does one know if one
# did not misspell the name of a variable in an if statement like::
#
#   if (SOME_VARBLE)
#     ...
#   endif()
#
# ?
#
#  If one misspelled the variable ``SOME_VARBLE`` (which is likely in this
#  case), then the if statement will always be false!  To avoid this problem
#  when one always expects that a variable is explicitly set, instead do::
#
#   assert_defined(SOME_VARBLE)
#   if (SOME_VARBLE)
#     ...
#   endif()
#
# Now if one misspells this variable, then CMake will asset and stop
# processing.  This is not a perfect solution since one can misspell the
# variable name in the following if statement but typically one would always
# just copy and paste between the two statements so these names are always the
# same.  This is the best that can be done in CMake unfortunately to catch
# usage of misspelled undefined variables.
#
function(assert_defined VARS)
  foreach(VAR ${VARS})
    if(NOT DEFINED ${VAR})
      message(SEND_ERROR "Error, the variable ${VAR} is not defined!")
    endif()
  endforeach()
endfunction()

# ToDo: The VARS arg This really needs to be replaced with ${ARGV}.  I fear
# that only the first arg passed in is asserted.  However, to change this now
# is breaking backward compatibility.
