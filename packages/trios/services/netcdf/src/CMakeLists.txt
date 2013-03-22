
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


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)


SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")

# Need to include the source dir so we can find aggregation.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# The rpcgen files should be treated as system files to avoid annoying warnings
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/netcdf_args.x)


APPEND_SET(NOINSTHEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.h
  aggregation.h
  create_subchunks.h
  io_timer.h
  NcAttInfo.h
  NcDimInfo.h
  NcFileInfo.h
  NcGroupInfo.h
  NcVarInfo.h
  netcdf_client_private.h
  netcdf_debug.h
)

APPEND_SET(LIB_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.c
    netcdf_debug.c
    NcFileInfo.cpp
    NcAttInfo.cpp
    NcVarInfo.cpp
    NcDimInfo.cpp
    NcGroupInfo.cpp
    netcdf_client.cpp
    ncmpi_client.cpp
    create_subchunks.cpp
)

############# libtrios_netcdf_client ##########################

TRIBITS_ADD_LIBRARY(
   trios_netcdf_client
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${LIB_SOURCES}
   TESTONLY
)

# trios_netcdf_client is created as a TESTONLY library which means it won't
# get installed.  forced the install here.
INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libtrios_netcdf_client.a
  DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
)
GLOBAL_SET(${PACKAGE_NAME}_NETCDF_LIBRARIES "trios_netcdf_client")

######################  PROGRAMS  ##################################################

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

APPEND_SET(SVC_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.c
    netcdf_server.cpp
    netcdf_debug.c
    aggregation.cpp
)

TRIBITS_ADD_EXECUTABLE(
  netcdf-server
  SOURCES ${SVC_SOURCES}
  DEPLIBS ${Trios_LIBRARIES}
  NOEXEPREFIX
)
