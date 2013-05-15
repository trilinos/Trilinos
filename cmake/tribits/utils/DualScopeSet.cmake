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


#
# Macro that sets a variable name both in the current scope and the
# parent scope.
#
# It turns out that when you call ADD_SUBDIRECTORY(someDir) that CMake
# actaully creates a copy of all of the regular non-cache varaibles in
# the current scope in order to create a new set of variables for the
# CMakeLists.txt file in 'someDir'.  This means that if you call
# SET(SOMEVAR Blah PARENT_SCOPE) that it will not affect the value of
# SOMEVAR in the current scope.  This macro therefore is designed to
# set the value of the variable in the current scope and the parent
# scope in one shot.
#
# Global variables are different.  When you move to a subordinate
# CMakeLists.txt file, a local copy of the variable is *not* created.
# If you set the value name locally, it will shadow the global
# variable.  However, if you set the globlal value with SET(SOMEVAR
# someValue CACHE INTERNAL ""), then the value will get changed in the
# current subordinate scope and in all parent scopes.
#

MACRO(DUAL_SCOPE_SET VARNAME)
  SET(${VARNAME} ${ARGN} PARENT_SCOPE)
  SET(${VARNAME} ${ARGN})
ENDMACRO()
