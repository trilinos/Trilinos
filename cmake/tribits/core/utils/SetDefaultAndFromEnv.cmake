# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(SetDefault)
include(PrintVar)


# @MACRO: set_default_and_from_env()
#
# Set a default value for a local variable and override from an environment
# variable of the same name if it is set.
#
# Usage::
#
#   set_default_and_from_env(<varName> <defaultVal>)
#
# First calls ``set_default(<varName> <defaultVal>)`` and then looks for an
# environment variable named ``<varName>``, and if non-empty then overrides
# the value of the local variable ``<varName>``.
#
# This macro is primarily used in CTest code to provide a way to pass in the
# value of CMake variables.  Older versions of ``ctest`` did not support the
# option ``-D <var>:<type>=<value>`` to allow variables to be set through the
# command-line like ``cmake`` always allowed.
#
macro(set_default_and_from_env  VAR  DEFAULT_VAL)

  set_default(${VAR} "${DEFAULT_VAL}")

  set(ENV_${VAR} $ENV{${VAR}})
  if (NOT "${ENV_${VAR}}" STREQUAL "")
    print_var(ENV_${VAR})
    set(${VAR} ${ENV_${VAR}})
  endif()

  print_var(${VAR})

endmacro()
