# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AssertDefined)
include(GlobalSet)


# @FUNCTION: remove_global_duplicates()
#
# Remove duplicate elements from a global list variable (removes boiler-plate
# code and errors).
#
# Usage::
#
#   remove_global_duplicates(<globalVarName>)
#
# This function is necessary in order to preserve the "global" nature of the
# variable.  If one just calls ``list(REMOVE_DUPLICATES ...)`` it will
# actually create a local variable of the same name and shadow the global
# variable!  That is a fun bug to track down!  The variable
# ``<globalVarName>`` must be defined before this function is called.  If
# ``<globalVarName>`` is actually not a global cache variable before this
# function is called it will be after it completes.
#
function(remove_global_duplicates VARNAME)
  assert_defined(${VARNAME})
  if (${VARNAME})
    set(TMP ${${VARNAME}})
    list(REMOVE_DUPLICATES TMP)
    global_set(${VARNAME} ${TMP})
  endif()
endfunction()
