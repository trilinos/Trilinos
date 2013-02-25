
# @HEADER
#  ************************************************************************
#
#                    Trios: Trilinos I/O Support
#                  Copyright 2011 Sandia Corporation
#
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#
#  1. Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
#  3. Neither the name of the Corporation nor the names of the
#  contributors may be used to endorse or promote products derived from
#  this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
#  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
#
#  *************************************************************************
# @HEADER

INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(programs)

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


  TRIBITS_ADD_EXECUTABLE(
    kill-service
    SOURCES nssi-kill.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )

# Gemini means Cray XE6.
# On Cray XE6, libtopomap requires the gpcd lib in addition to METIS and ParMETIS.
IF (TPL_ENABLE_Gemini   AND
    TPL_ENABLE_gpcd     AND
    TPL_ENABLE_METIS    AND
    TPL_ENABLE_ParMETIS AND
    TPL_ENABLE_LibTopoMap)
  TRIBITS_ADD_EXECUTABLE(
    nssi-placement
    SOURCES nssi-placement.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )
ENDIF()

# SPI means BG/P.
# On BG/P, libtopomap only requires METIS and ParMETIS.
IF (TPL_ENABLE_SPI      AND
    TPL_ENABLE_METIS    AND
    TPL_ENABLE_ParMETIS AND
    TPL_ENABLE_LibTopoMap)
  TRIBITS_ADD_EXECUTABLE(
    nssi-placement
    SOURCES nssi-placement.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )
ENDIF()

  TRIBITS_ADD_EXECUTABLE(
    ping-service
    SOURCES nssi-ping.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )

TRIBITS_SUBPACKAGE_POSTPROCESS()
