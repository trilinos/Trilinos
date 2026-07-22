# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: split()
#
# Split a string variable into a string array/list variable.
#
# Usage::
#
#   split("<inputStr>" "<sepStr>" <outputStrListVar>)
#
# The ``<sepStr>`` string is used with ``string(REGEX ...)`` to replace all
# occurrences of ``<sepStr>`` in ``<inputStr>`` with ``";"`` and writing into
# ``<outputStrListVar>``.
#
# WARNING: ``<sepStr>`` is interpreted as a regular expression (regex) so keep
# that in mind when considering special regex chars like ``'*'``, ``'.'``,
# etc!
#
function(split  INPUT_STRING  SEP_STR  OUTPUT_STRING_VAR)
  string(REGEX REPLACE "${SEP_STR}" ";" OUTPUT_STRING "${INPUT_STRING}")
  set(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
endfunction()
