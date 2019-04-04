# Copyright (c) 2013, Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of Sandia Corporation nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 


#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS_IMPL "")
SET(SOURCES_IMPL "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

FILE(GLOB HEADERS_INTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/internal/*.hpp)
FILE(GLOB SOURCES_INTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/internal/*.cpp)

FILE(GLOB HEADERS_SEARCH_TOLERANCE ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance/*.hpp)

FILE(GLOB HEADERS_SEARCH_TOLERANCE_ALGS ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance_algs/*.hpp)
FILE(GLOB SOURCES_SEARCH_TOLERANCE_ALGS ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance_algs/*.cpp)

FILE(GLOB SOURCES_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/balance_main/*.cpp)

FILE(GLOB SOURCES_M2N_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/m2n_main/*.cpp)

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_balance_lib
  NOINSTALLHEADERS ${HEADERS} ${HEADERS_INTERNAL} ${HEADERS_SEARCH_TOLERANCE} ${HEADERS_SEARCH_TOLERANCE_ALGS}
  SOURCES ${SOURCES} ${SOURCES_INTERNAL} ${SOURCES_SEARCH_TOLERANCE_ALGS}
  DEPLIBS
  )

# TRIBITS_ADD_LIBRARY(
#   stk_balance_search_tolerance_algs
#   NOINSTALLHEADERS ${HEADERS_SEARCH_TOLERANCE_ALGS} 
#   SOURCES ${SOURCES_SEARCH_TOLERANCE_ALGS}
#   DEPLIBS
#   )

TRIBITS_ADD_EXECUTABLE(
    stk_balance
    SOURCES ${SOURCES_MAIN}
    DEPLIBS
    NOEXEPREFIX INSTALLABLE
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    stk_balance_m2n
    SOURCES ${SOURCES_M2N_MAIN}
    DEPLIBS
    NOEXEPREFIX INSTALLABLE
    )

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance)
INSTALL(FILES ${HEADERS_INTERNAL} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/internal)
INSTALL(FILES ${HEADERS_SEARCH_TOLERANCE} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/search_tolerance)
INSTALL(FILES ${HEADERS_SEARCH_TOLERANCE_ALGS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/search_tolerance_algs)

