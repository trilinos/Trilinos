# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


function( get_last_dir_name  INPUT_DIR  OUTPUT_DIR_VAR )
  file(TO_CMAKE_PATH "${INPUT_DIR}" STANDARD_INPUT_DIR)
  string(REGEX REPLACE "/.+/(.+)" "\\1" LOCAL_OUTPUT_DIR "${STANDARD_INPUT_DIR}")
  set(${OUTPUT_DIR_VAR} "${LOCAL_OUTPUT_DIR}" PARENT_SCOPE)
endfunction()
