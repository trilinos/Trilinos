# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

include(TribitsPkgExportCacheVars)
include(GlobalSet)


# @FUNCTION: tribits_add_option_and_define()
#
# Add an option and an optional macro define variable in one shot.
#
# Usage::
#
#  tribits_add_option_and_define( <userOptionName>  <macroDefineName>
#    "<docStr>"  <defaultValue> [NONCACHE])
#
# This macro sets the user cache ``BOOL`` variable ``<userOptionName>`` and if
# it is true, then sets the global (internal cache) macro define variable
# ``<macroDefineName>`` to ``ON``, and otherwise sets it to ``OFF``.  If
# ``NONCACHE`` is passed in, then ``<macroDefineName>`` is set as a non-cache
# local variable instead of a cache variable.
#
# This is designed to make it easy to add a user-enabled option to a
# configured header file and have the define set in one shot.  This would
# require that the package's configure file (see `tribits_configure_file()`_)
# have the line::
#
#   #cmakedefine <macroDefineName>
#
# NOTE: This also calls `tribits_pkg_export_cache_var()`_ to export the
# variables ``<userOptionName>`` and ``<macroDefineName>`` (when ``NONCACHE``
# is **not** passed).  This also requires that local variables with the same
# names of these cache variables not be assigned with a different value from
# these cache variables.  If they are, then an error will occur later when
# these variables are read.
#
# NOTE: The define var name ``<macroDefineName>`` can be empty "" in which
# case all logic related to ``<macroDefineName>`` is skipped.  (But in this
# case, it would be better to just call::
#
#   set(<userOptionName> <defaultValue> CACHE BOOL "<docStr>")
#
function(tribits_add_option_and_define  USER_OPTION_NAME  MACRO_DEFINE_NAME
    DOCSTRING  DEFAULT_VALUE
  )
  # Assert inputs
  set(set_MACRO_DEFINE_NAME_noncache OFF)
  if ("${ARGN}" STREQUAL "NONCACHE")
    set(set_MACRO_DEFINE_NAME_noncache ON)
  elseif (NOT "${ARGN}" STREQUAL "")
    message(SEND_ERROR "Error, optional argument"
      " '${ARGN}' can only be 'NONCACHE' or empty!")
  endif()
  # Assert these are not already defined local vars
  if (DEFINED ${USER_OPTION_NAME})
    if (NOT DEFINED CACHE{${USER_OPTION_NAME}})
      message(SEND_ERROR "ERROR: The option"
	" ${USER_OPTION_NAME}='${${USER_OPTION_NAME}}'"
        " is already defined as a non-cache variable!")
    endif()
  endif()
  if (DEFINED ${MACRO_DEFINE_NAME})
    if (NOT DEFINED CACHE{${MACRO_DEFINE_NAME}})
      message(SEND_ERROR "ERROR: The macro define name"
	" ${MACRO_DEFINE_NAME}='${${MACRO_DEFINE_NAME}}'"
        " is already defined as a non-cache variable!")
    endif()
  endif()
  # Set ${USER_OPTION_NAME} as exported cache var
  set( ${USER_OPTION_NAME} "${DEFAULT_VALUE}" CACHE BOOL "${DOCSTRING}" )
  tribits_pkg_export_cache_var(${USER_OPTION_NAME})
  # Set ${MACRO_DEFINE_NAME} as local or exported cache var
  if(NOT "${MACRO_DEFINE_NAME}" STREQUAL "")
    if(${USER_OPTION_NAME})
      set(macroDefineValue ON)
    else()
      set(macroDefineValue OFF)
    endif()
    if (set_MACRO_DEFINE_NAME_noncache)
      set(${MACRO_DEFINE_NAME} ${macroDefineValue} PARENT_SCOPE)
    else()
      global_set(${MACRO_DEFINE_NAME} ${macroDefineValue})
      tribits_pkg_export_cache_var(${MACRO_DEFINE_NAME})
    endif()
  endif()
endfunction()
