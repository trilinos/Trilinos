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

function(tribits_set_linker_language_from_arg  TARGET_NAME_IN  LINKER_LANGUAGE_IN)

  if( TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP)
    message("tribits_set_linker_language_from_arg( '${TARGET_NAME_IN}' '${LINKER_LANGUAGE_IN}' )")
  endif()

  if(LINKER_LANGUAGE_IN)
    set(LINKER_LANGUAGE ${LINKER_LANGUAGE_IN})
  elseif (${PROJECT_NAME}_ENABLE_CXX)
    set(LINKER_LANGUAGE CXX)
  elseif(${PROJECT_NAME}_ENABLE_C)
    set(LINKER_LANGUAGE C)
  else()
    set(LINKER_LANGUAGE)
  endif()

  if (LINKER_LANGUAGE)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE
        OR TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP
      )
      message("-- Setting linker language for target '${TARGET_NAME_IN}' to '${LINKER_LANGUAGE}'")
    endif()

    set_property(
      TARGET ${TARGET_NAME_IN}
      APPEND PROPERTY LINKER_LANGUAGE ${LINKER_LANGUAGE}
      )

  endif()

endfunction()


