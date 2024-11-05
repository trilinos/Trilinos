# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AppendSet)


# @MACRO: append_glob()
#
# Utility macro that does a ``file(GLOB ...)`` and appends to an existing list
# (removes boiler-plate code).
#
# Usage::
#
#   append_glob(<fileListVar> <glob0> <glob1> ...)
#
# On output, ``<fileListVar>`` will have the list of glob files appended.
#
macro(append_glob VAR)
  file(GLOB LOCAL_TMP_VAR ${ARGN})
  append_set(${VAR} ${LOCAL_TMP_VAR})
endmacro()
