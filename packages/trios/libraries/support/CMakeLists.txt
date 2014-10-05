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

TRIBITS_SUBPACKAGE(support)

# A) Package-specific configuration options
#

# Initializations
SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")
SET(MODULES "")


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

#INCLUDE_DIRECTORIES(BEFORE SYSTEM ${${PARENT_PACKAGE}_BINARY_DIR})

#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

message(STATUS "${${PARENT_PACKAGE_NAME}_BINARY_DIR}/Trios_config.h")

APPEND_SET(HEADERS
  ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/Trios_config.h
)

############# LOGGER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/logger)

APPEND_SET(HEADERS
  logger/Trios_logger.h
)

APPEND_SET(SOURCES
  logger/logger.cpp
)


############# THREADS ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/threads)

APPEND_SET(HEADERS
  threads/Trios_threads.h
  threads/Trios_threads_types.h
)

APPEND_SET(SOURCES
  threads/trios_threads.cpp
)

############# TIMER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/timer)

APPEND_SET(HEADERS
  timer/Trios_timer.h
)

APPEND_SET(SOURCES
  timer/timer.c
  timer/Trios_Timer.cpp
)


############# SIGNAL ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/signal)

APPEND_SET(HEADERS
  signal/Trios_signal.h
)

APPEND_SET(SOURCES
  signal/signal.cpp
)


########### Trace  ###############



INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/trace)

APPEND_SET(HEADERS
  trace/Trios_trace.h
)

APPEND_SET(NOINSTALLHEADERS
  trace/SDDF.h
  trace/TraceFile.h
  trace/TraceGroup.h
)

APPEND_SET(SOURCES
  trace/trace.cc
)

# See if PABLO is enabled
IF (${PACKAGE_NAME}_ENABLE_Pablo)
  APPEND_SET(SOURCES trace/SDDF.cc)
ENDIF (${PACKAGE_NAME}_ENABLE_Pablo)


############# XDR ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/xdr)

APPEND_SET(HEADERS
  xdr/Trios_xdr.h
)

########################### TRIOS SUPPORT LIBRARY ###############

TRIBITS_ADD_LIBRARY(
   trios_support
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${SOURCES}
   DEPLIBS ${DEPLIBS}
   )


########### install targets ###############

# INSTALL(FILES ${MODULES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

TRIBITS_SUBPACKAGE_POSTPROCESS()
