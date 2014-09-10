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

FUNCTION(TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG  TARGET_NAME_IN  LINKER_LANGUAGE_IN)

  IF( TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP)
    MESSAGE("TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG( '${TARGET_NAME_IN}' '${LINKER_LANGUAGE_IN}' )")
  ENDIF()

  IF(LINKER_LANGUAGE_IN)
    SET(LINKER_LANGUAGE ${LINKER_LANGUAGE_IN})
  ELSEIF (${PROJECT_NAME}_ENABLE_CXX)
    SET(LINKER_LANGUAGE CXX)
  ELSEIF(${PROJECT_NAME}_ENABLE_C)
    SET(LINKER_LANGUAGE C)
  ELSE()
    SET(LINKER_LANGUAGE)
  ENDIF()

  IF (LINKER_LANGUAGE)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE OR TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP)
      MESSAGE("-- Setting linker language for target '${TARGET_NAME_IN}' to '${LINKER_LANGUAGE}'") 
    ENDIF()

    SET_PROPERTY(
      TARGET ${TARGET_NAME_IN}
      APPEND PROPERTY LINKER_LANGUAGE ${LINKER_LANGUAGE}
      )
  ENDIF()

ENDFUNCTION()


