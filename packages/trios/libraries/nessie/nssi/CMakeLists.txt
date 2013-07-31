
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
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TriosProcessXDR)

SET(NSSI_HEADERS "")
SET(NSSI_NONINSTALLHEADERS "")
SET(NSSI_SOURCES "")


############# COMMON ##########################

TRIBITS_SUBPACKAGE(nssi)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


############# NSSI ##########################

TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/Trios_nssi_types_xdr.x)
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/nssi_service_args.x)

APPEND_SET(NSSI_NOINSTHEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.h
  nssi_opcodes.h
  nssi_trace.h
)

APPEND_SET(NSSI_HEADERS
  Trios_nssi_rpc.h
  Trios_nssi_types.h
  Trios_nssi_debug.h
  Trios_nssi_fprint_types.h
  Trios_nssi_xdr.h
  Trios_nssi_TypeManip.h
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.h
)

APPEND_SET(NSSI_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.c
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.c
  fprint_types.cpp
  nssi_common.cpp
  nssi_debug.c
  nssi_xdr.cpp
)

############# CLIENT ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../client)

APPEND_SET(NSSI_HEADERS
  ../client/Trios_nssi_client.h
  ../client/Trios_nssi_request.h
)

APPEND_SET(NSSI_SOURCES
  ../client/nssi_client.cpp
)

############# SERVER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../server)

APPEND_SET(NSSI_HEADERS
  ../server/Trios_nssi_server.h
)

APPEND_SET(NSSI_SOURCES
  ../server/nssi_server.cpp
)

############# libtrios_nssi.a ##########################

TRIBITS_ADD_LIBRARY(
   trios_nssi
   HEADERS ${NSSI_HEADERS}
   NOINSTALLHEADERS ${NSSI_NOINSTALLHEADERS}
   SOURCES ${NSSI_SOURCES}
   DEPLIBS ${DEPLIBS}
)


TRIBITS_SUBPACKAGE_POSTPROCESS()
