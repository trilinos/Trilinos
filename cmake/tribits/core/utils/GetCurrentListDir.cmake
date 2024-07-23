# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Get the directory containing the file currently being processed.
#
# This is equivalent to getting the value of CMAKE_CURRENT_LIST_DIR in
# cmake versions greater than 2.8.4, but we provide this wrapper macro
# for compatibility. If called outside of a function or macro, this
# will return the directory of the calling file. If called within a
# function or macro, this will return the directory containing the
# caller.
function(get_current_list_dir output_variable)
  if(NOT DEFINED CMAKE_CURRENT_LIST_DIR)
    get_filename_component(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
  endif()
  set(${output_variable} "${CMAKE_CURRENT_LIST_DIR}" PARENT_SCOPE)
endfunction()