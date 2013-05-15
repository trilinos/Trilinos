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

INCLUDE(ParseVariableArguments)

#This function is to override the standard behavior of include_directories.
#  We are overriding the default behavior for installation testing, this allows
#us to ensure that include directories will not be inadvertently added to the
#build lines for tests during installation testing. Normally we want the include
#directories to be handled as cmake usually does.However during installation
#testing we do not want most of the include directories to be used as the majority
#of the files should come from the installation we are building against. There is
#an exception to this and that is when there are test only headers that are needed.
#For that case we allow people to set "REQUIRED_DURING_INSTALLATION_TESTING" to
#tell us that this include directory does need to be set for instaltion testing.
FUNCTION(INCLUDE_DIRECTORIES)
  PARSE_ARGUMENTS(
    PARSE #prefix
    "" # Lists
    "REQUIRED_DURING_INSTALLATION_TESTING" #Options
    ${ARGN}
    )

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING)
    _INCLUDE_DIRECTORIES(${PARSE_DEFAULT_ARGS})
  ENDIF()
ENDFUNCTION()
