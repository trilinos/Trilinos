# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(AppendStringVar)


function(tribits_strip_comments_from_cmake_cache_file  INPUT_FILE  OUTPUT_FILE)
  execute_process(
    COMMAND cat "${INPUT_FILE}"
    COMMAND grep -v "^#"
    COMMAND grep -v "^//"
    COMMAND grep -v "^$"
    OUTPUT_FILE "${OUTPUT_FILE}"
    )
endfunction()
