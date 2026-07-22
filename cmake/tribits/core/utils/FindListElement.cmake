# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# Function that searches for an element in a list and if found returns
# true.  Otherwise returns false.
#

function(find_list_element LIST_NAME ELEMENT_VAL ELEMENT_FOUND_OUT)
  #message("FIND_LIST_ELEMENT: ${LIST_NAME} ${ELEMENT_VAL}")
  #print_var(${LIST_NAME})
  list(FIND ${LIST_NAME} ${ELEMENT_VAL} ELEMENT_IDX)
  if (ELEMENT_IDX EQUAL -1)
    set(ELEMENT_FOUND FALSE)
  else()
    set(ELEMENT_FOUND TRUE)
  endif()
  #print_var(ELEMENT_FOUND)
  set(${ELEMENT_FOUND_OUT} ${ELEMENT_FOUND} PARENT_SCOPE)
endfunction()
