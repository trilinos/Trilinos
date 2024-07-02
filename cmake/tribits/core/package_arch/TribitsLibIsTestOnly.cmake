# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()


# @FUNCTION: tribits_set_lib_is_testonly()
#
# See if a library is a TESTONLY library
#
# Usage::
#
#   tribits_set_lib_is_testonly(<libName>)
#
# This sets the ``TRIBITS_TESTONLY_LIB`` on the library target ``<libName>``.
#
function(tribits_set_lib_is_testonly  libName)
  set_target_properties(${libName}  PROPERTIES  TRIBITS_TESTONLY_LIB  TRUE)
endfunction()


# @FUNCTION: tribits_lib_is_testonly()
#
# See if a library is a TESTONLY library
#
# Usage::
#
#   tribits_lib_is_testonly(<libName> <libIsTestOnlyOut>)
#
# This will only return ``TRUE`` in `` <libIsTestOnlyOut>`` if ``<libName>``
# is a target and the target property ``TRIBITS_TESTONLY_LIB`` is set to
# ``TRUE``.
#
function(tribits_lib_is_testonly  libName  libIsTestOnlyOut)
  if (TARGET ${libName})
    get_target_property(libIsTestOnly ${libName} TRIBITS_TESTONLY_LIB)
  else()
    set(libIsTestOnly FALSE)
  endif()
  set(${libIsTestOnlyOut} ${libIsTestOnly} PARENT_SCOPE)
endfunction()
