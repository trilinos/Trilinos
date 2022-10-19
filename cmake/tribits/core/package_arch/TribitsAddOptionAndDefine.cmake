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


# @MACRO: tribits_add_option_and_define()
#
# Add an option and an optional macro define variable in one shot.
#
# Usage::
#
#  tribits_add_option_and_define( <userOptionName>  <macroDefineName>
#    "<docStr>"  <defaultValue> )
#
# This macro sets the user cache ``BOOL`` variable ``<userOptionName>`` and if
# it is true, then sets the global (internal cache) macro define variable
# ``<macroDefineName>`` to ``ON``, and otherwise sets it to ``OFF``.  This is
# designed to make it easy to add a user-enabled option to a configured header
# file and have the define set in one shot.  This would require that the
# package's configure file (see `tribits_configure_file()`_) have the line::
#
#   #cmakedefine <macroDefineName>
#
# NOTE: This also calls `tribits_pkg_export_cache_var()`_ to export the
# variables ``<userOptionName>`` and ``<macroDefineName>``.  This also
# requires that local variables with the same names of these cache variables
# not be assigned with a different value from these cache variables.  If they
# are, then an error will occur later when these variables are read.
#
# NOTE: The define var name ``<macroDefineName>`` can be empty "" in which
# case all logic related to ``<macroDefineName>`` is skipped.  (But in this
# case, it would be better to just call::
#
#   set(<userOptionName> <defaultValue> CACHE BOOL "<docStr>")
#
macro(tribits_add_option_and_define  USER_OPTION_NAME  MACRO_DEFINE_NAME
  DOCSTRING  DEFAULT_VALUE
  )
  #message("TRIBITS_ADD_OPTION_AND_DEFINE: '${USER_OPTION_NAME}' '${MACRO_DEFINE_NAME}' '${DEFAULT_VALUE}'")
  set( ${USER_OPTION_NAME} "${DEFAULT_VALUE}" CACHE BOOL "${DOCSTRING}" )
  if(NOT "${MACRO_DEFINE_NAME}" STREQUAL "")
    if(${USER_OPTION_NAME})
      global_set(${MACRO_DEFINE_NAME} ON)
    else()
      global_set(${MACRO_DEFINE_NAME} OFF)
    endif()
  endif()
  tribits_pkg_export_cache_var(${USER_OPTION_NAME})
  if(NOT "${MACRO_DEFINE_NAME}" STREQUAL "")
    tribits_pkg_export_cache_var(${MACRO_DEFINE_NAME})
  endif()
endmacro()

# 2008/10/05: rabartl: ToDo: Add an option to automatically add the macro
# define to any XXX_config.h file that gets configured by
# tribits_configure_file(...).  This will help to eliminate
# duplication.
