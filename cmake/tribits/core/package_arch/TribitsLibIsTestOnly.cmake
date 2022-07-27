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

include_guard()


# @FUNCTION: tribits_set_lib_is_testonly()
#
# See if a library is a TESTONLY library
#
# Usage::
#
#   tribits_set_lib_is_testonly(<libName>)
#
# This sets the ``TRIBITS_TESTONLY_LIB`` on the library target ``<libName>``.
#
function(tribits_set_lib_is_testonly  libName)
  set_target_properties(${libName}  PROPERTIES  TRIBITS_TESTONLY_LIB  TRUE)
endfunction()


# @FUNCTION: tribits_lib_is_testonly()
#
# See if a library is a TESTONLY library
#
# Usage::
#
#   tribits_lib_is_testonly(<libName> <libIsTestOnlyOut>)
#
# This will only return ``TRUE`` in `` <libIsTestOnlyOut>`` if ``<libName>``
# is a target and the target property ``TRIBITS_TESTONLY_LIB`` is set to
# ``TRUE``.
#
function(tribits_lib_is_testonly  libName  libIsTestOnlyOut)
  if (TARGET ${libName})
    get_target_property(libIsTestOnly ${libName} TRIBITS_TESTONLY_LIB)
  else()
    set(libIsTestOnly FALSE)
  endif()
  set(${libIsTestOnlyOut} ${libIsTestOnly} PARENT_SCOPE)
endfunction()
