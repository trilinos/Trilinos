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


INCLUDE(TribitsAddTestHelpers)


# Allow for a maximum of 20 (0 through 19) test commands
SET(PACKAGE_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX 19)


#
# Join arguments together to add to a SET(...) statement to passed to
# EXECUTE_PROCESS(...)
#

FUNCTION(TRIBITS_JOIN_EXEC_PROCESS_SET_ARGS  OUTPUT_STRING_VAR)
  SET(OUTPUT_STRING "")
  FOREACH(STRING_VAL_RAW ${ARGN})
    # Remove quotes around arguments because CTest does not need them
    STRING(REGEX REPLACE "\"" "" STRING_VAL "${STRING_VAL_RAW}")
    IF (OUTPUT_STRING STREQUAL "")
      SET(OUTPUT_STRING "\"${STRING_VAL}\"")
    ELSE()
      SET(OUTPUT_STRING "${OUTPUT_STRING} \"${STRING_VAL}\"")
    ENDIF()
  ENDFOREACH()
  SET(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
ENDFUNCTION()


#
# Unit test helper function for TRIBITS_ADD_ADVANCED_TEST(...) that resets
# state before calling TRIBITS_ADD_ADVANCED_TEST(...) in unit test mode.
#

FUNCTION(TRIBITS_ADD_ADVANCED_TEST_UNITTEST_RESET)

  GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_UNITTEST TRUE)

  GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS)

  FOREACH( TEST_CMND_IDX RANGE ${PACKAGE_ADD_ADVANCED_TEST_MAX_NUM_TEST_CMND_IDX})
    GLOBAL_SET(PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_${TEST_CMND_IDX})
  ENDFOREACH()

ENDFUNCTION()
