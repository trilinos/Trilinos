# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

INCLUDE(GlobalSet)

MACRO(TRIBITS_ADD_OPTION_AND_DEFINE USER_OPTION_NAME MACRO_DEFINE_NAME DOCSTRING DEFAULT_VALUE)
  #MESSAGE("TRIBITS_ADD_OPTION_AND_DEFINE: '${USER_OPTION_NAME}' '${MACRO_DEFINE_NAME}' '${DEFAULT_VALUE}'")
  SET( ${USER_OPTION_NAME} "${DEFAULT_VALUE}" CACHE BOOL "${DOCSTRING}" )
  IF(NOT ${MACRO_DEFINE_NAME} STREQUAL "")
    IF(${USER_OPTION_NAME})
      GLOBAL_SET(${MACRO_DEFINE_NAME} ON)
    ELSE()
      GLOBAL_SET(${MACRO_DEFINE_NAME} OFF)
    ENDIF()
  ENDIF()
ENDMACRO()

# 2008/10/05: rabartl: ToDo: Add an option to automatically add the macro
# define to any XXX_config.h file that gets configured by
# TRIBITS_CONFIGURE_FILE(...).  This will help to eliminate
# duplication.
