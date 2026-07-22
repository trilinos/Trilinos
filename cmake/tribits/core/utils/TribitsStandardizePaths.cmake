# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: tribits_standardize_abs_paths()
#
# Function uses get_filename_component() to standardize a list of paths to be
# absolute paths.
#
# Usage::
#
#   tribits_standardize_abs_paths(<pathsListvar> <path0> <path1> ...)
#
# On output, ``<pathsListLvar>`` will be set to the list of paths
#
function(tribits_standardize_abs_paths  PATHS_LIST_VAR_OUT)
  set(PATHS_LIST)
  foreach(PATH_I ${ARGN})
    #print_var(PATH_I)
    get_filename_component(STD_ABS_PATH_I "${PATH_I}" ABSOLUTE)
    #print_var(STD_ABS_PATH_I)
    list(APPEND PATHS_LIST "${STD_ABS_PATH_I}")
  endforeach()
  set(${PATHS_LIST_VAR_OUT} ${PATHS_LIST} PARENT_SCOPE)
endfunction()

