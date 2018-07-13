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

#
# @FUNCTION: TRIBITS_STANDARDIZE_ABS_PATHS()
#
# Function uses GET_FILENAME_COMPONENT() to standardize a list of paths to be
# absolute paths.
#
# Usage::
#
#   TRIBITS_STANDARDIZE_ABS_PATHS(<pathsListvar> <path0> <path1> ...)
#
# On output, ``<pathsListLvar>`` will be set to the list of paths
#
FUNCTION(TRIBITS_STANDARDIZE_ABS_PATHS  PATHS_LIST_VAR_OUT)
  SET(PATHS_LIST)
  FOREACH(PATH_I ${ARGN})
    #PRINT_VAR(PATH_I)
    GET_FILENAME_COMPONENT(STD_ABS_PATH_I "${PATH_I}" ABSOLUTE)
    #PRINT_VAR(STD_ABS_PATH_I)
    LIST(APPEND PATHS_LIST "${STD_ABS_PATH_I}")
  ENDFOREACH()
  SET(${PATHS_LIST_VAR_OUT} ${PATHS_LIST} PARENT_SCOPE)
ENDFUNCTION()

