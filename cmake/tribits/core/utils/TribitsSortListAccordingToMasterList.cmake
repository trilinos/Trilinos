# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include("${CMAKE_CURRENT_LIST_DIR}/../utils/PrintVar.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/AppendSet.cmake")


# Do an in-place sort of a list of items according to the ordering in a master
# list.
#
# NOTE: This function has worst-case complexity N*n where N is the number of
# elements in the ``<masterList>`` and n is the number of elements in the
# ``<listVarInout>`` list.
#
function(tribits_sort_list_according_to_master_list  masterList  listVarInOut)

  set(sortedList)

  foreach(item ${masterList})
    list(FIND ${listVarInOut} ${item} itemIdx)
     if (NOT itemIdx EQUAL -1)
      list(APPEND sortedList ${item})
    endif()
  endforeach()

  set(${listVarInOut} ${sortedList} PARENT_SCOPE)

endfunction()
