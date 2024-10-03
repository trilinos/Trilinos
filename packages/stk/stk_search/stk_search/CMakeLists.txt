# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
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
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
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
SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS_KDTREE  
     ${CMAKE_CURRENT_SOURCE_DIR}/kdtree/*.hpp)
FILE(GLOB HEADERS_MORTON  
     ${CMAKE_CURRENT_SOURCE_DIR}/morton_lbvh/*.hpp)
FILE(GLOB HEADERS_ARBORX  
     ${CMAKE_CURRENT_SOURCE_DIR}/arborx/*.hpp)
FILE(GLOB HEADERS *.hpp)

FILE(GLOB SOURCES_KDTREE
     ${CMAKE_CURRENT_SOURCE_DIR}/kdtree/*.cpp)
FILE(GLOB SOURCES_MORTON
     ${CMAKE_CURRENT_SOURCE_DIR}/morton_lbvh/*.cpp)

FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_LIBRARY(stk_search
    NOINSTALLHEADERS ${HEADERS} ${HEADERS_KDTREE} ${HEADERS_MORTON}
    SOURCES ${SOURCES} ${SOURCES_KDTREE} ${SOURCES_MORTON}
  )
else()
  add_library(stk_search ${SOURCES} ${SOURCES_KDTREE} ${SOURCES_MORTON})

  target_link_libraries(stk_search PUBLIC stk_math)
  target_link_libraries(stk_search PUBLIC stk_util_parallel)

  if(STK_HAS_ARBORX)
    target_link_libraries(stk_search PUBLIC ArborX::ArborX)
    INSTALL(FILES ${HEADERS_ARBORX} DESTINATION
        ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_search/arborx)
  endif()

endif()

target_include_directories(stk_search PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_search>
  $<INSTALL_INTERFACE:include>
)

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_search)
INSTALL(FILES ${HEADERS_KDTREE} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_search/kdtree)
INSTALL(FILES ${HEADERS_MORTON} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_search/morton_lbvh)

if(NOT HAVE_STK_Trilinos)
  INSTALL(TARGETS stk_search EXPORT stkTargets DESTINATION ${STK_INSTALL_LIBDIR})
endif()
