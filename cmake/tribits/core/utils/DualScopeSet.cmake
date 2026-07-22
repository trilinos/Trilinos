# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: dual_scope_set()
#
# Macro that sets a variable name both in the current scope and the
# parent scope.
#
# Usage::
#
#    dual_scope_set(<varName> [other args])
#
# It turns out that when one calls ``add_subdirectory(<someDir>)`` or enters a
# ``FUNCTION`` that CMake actually creates a copy of all of the regular
# non-cache variables in the current scope in order to create a new set of
# variables for the ``CMakeLists.txt`` file in ``<someDir>``.  This means that
# if you call ``set(SOMEVAR Blah PARENT_SCOPE)`` that it will not affect the
# value of ``SOMEVAR`` in the current scope!  This macro therefore is designed
# to set the value of the variable in the current scope and the parent scope
# in one shot to avoid confusion.
#
# Global variables are different.  When one moves to a subordinate
# ``CMakeLists.txt`` file or enters a ``FUNCTION``, then a local copy of the
# variable is *not* created.  If one sets the variable locally, it will shadow
# the global variable.  However, if one sets the global cache value with
# ``set(SOMEVAR someValue CACHE INTERNAL "")``, then the value will get
# changed in the current subordinate scope and in all parent scopes all in one
# shot!
#
macro(dual_scope_set VARNAME)
  set(${VARNAME} ${ARGN} PARENT_SCOPE)
  set(${VARNAME} ${ARGN})
endmacro()
