# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsTplFindIncludeDirsAndLibraries)
include(TribitsDeprecatedHelpers)

function(tribits_tpl_declare_libraries TPL_NAME)
  tribits_deprecated_command(tribits_tpl_declare_libraries
    MESSAGE
    "Use tribits_tpl_find_include_dirs_and_libraries() instead."
    "  Make this change in the file:\n"
    "  ${${TPL_NAME}_FINDMOD}\n"
    "which is pointed to by the file:\n"
    "  ${${TPL_NAME}_TPLS_LIST_FILE}\n"
    )
  tribits_tpl_find_include_dirs_and_libraries(${TPL_NAME} ${ARGN})
endfunction()
