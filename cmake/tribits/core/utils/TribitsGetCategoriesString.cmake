# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


function(tribits_get_categories_string  CATEGORIES_IN  CATEGORIES_STR_OUT)
  if (CATEGORIES_IN)
    string(REPLACE ";" ", " CATEGORIES_STR "${CATEGORIES_IN}")
  else()
    set(CATEGORIES_STR "BASIC")
  endif()
  set(${CATEGORIES_STR_OUT} "${CATEGORIES_STR}" PARENT_SCOPE)
endfunction()