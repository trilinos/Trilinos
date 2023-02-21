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


# @MACRO: tribits_set_and_inc_dirs()
#
# Set a variable to an include directory and call
# `tribits_include_directories()`_ (removes boiler-plate code).
#
# Usage:
#
#   tribits_set_and_inc_dirs(<dirVarName> <includeDir>)
#
# On output, this sets ``<dirVarName>`` to ``<includeDir>`` in the local scope
# and calls ``tribits_include_directories(<includeDir>)``.
#
macro(tribits_set_and_inc_dirs  dirVarName  includeDir)
  set(${dirVarName} ${includeDir})
  tribits_include_directories(${${dirVarName}})
endmacro()


# Deprecated!  Use tribits_set_and_inc_dirs() instead!
#
macro(set_and_inc_dirs  DIR_VAR_NAME  INCLUDE_DIR)
  tribits_deprecated_command(set_and_inc_dirs
    MESSAGE "Use tribits_set_and_inc_dirs() instead." )
  set(${DIR_VAR_NAME} ${INCLUDE_DIR})
  tribits_include_directories(${${DIR_VAR_NAME}})
endmacro()
