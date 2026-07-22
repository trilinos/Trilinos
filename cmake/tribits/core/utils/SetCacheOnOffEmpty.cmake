# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: set_cache_on_off_empty()
#
# Usage::
#
#   set_cache_on_off_empty(<varName> <initialVal> "<docString>" [FORCE])
#
# Sets a special string cache variable with possible values "", "ON", or
# "OFF".  This results in a nice drop-down box in the CMake cache manipulation
# GUIs.
#
function(set_cache_on_off_empty VAR INITIAL_VALUE DOCSTR)
  set(FORCE_ARG)
  foreach(ARG ${ARGN})
    if (ARG STREQUAL FORCE)
      set(FORCE_ARG FORCE)
    else()
      message(FATAL_ERROR "set_cache_on_off_empty(...): Error, last arg '${ARG}' is"
        "invalid!  Must be 'FORCE' or nothing." )
    endif()
  endforeach()
  set( ${VAR} "${INITIAL_VALUE}" CACHE STRING "${DOCSTR}" ${FORCE_ARG})
  set_property(CACHE ${VAR} PROPERTY STRINGS "" "ON" "OFF")
endfunction()
