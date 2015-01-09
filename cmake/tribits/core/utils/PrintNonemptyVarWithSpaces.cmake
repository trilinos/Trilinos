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

INCLUDE(AssertDefined)
INCLUDE(AppendStringVarWithSep)

#
# @FUNCTION: PRINT_NONEMPTY_VAR_WITH_SPACES()
#
# Print a defined variable giving its name then value printed with spaces
# instead of ``';'`` but only if it is not empty.
#
# Usage::
#
#    PRINT_NONEMPTY_VAR_WITH_SPACES(<varName>  <printedVarInOut>)
#
# Prints the variable as::
#
#    <varName>: <ele0> <ele1> ...
#
# If ``$<printedVarInOut>`` is ``TRUE`` on input, then the varible is not
# touched. If however, the varible ``$<printedVarInOut>`` is not ``TRUE`` on
# input, then it is set to ``TRUE`` on output.
#
FUNCTION(PRINT_NONEMPTY_VAR_WITH_SPACES  VARIBLE_NAME  PRINTED_VAR_OUT)
  ASSERT_DEFINED(VARIBLE_NAME)
  IF (NOT "${${VARIBLE_NAME}}" STREQUAL "")
    STRING(REPLACE ";" " " OUTSTR "${${VARIBLE_NAME}}")
    MESSAGE("-- ${VARIBLE_NAME}: ${OUTSTR}")
    SET(${PRINTED_VAR_OUT} TRUE PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()
