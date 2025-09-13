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

SET(HEADERS_IMPL "")
SET(SOURCES_IMPL "")

FILE(GLOB HEADERS CONFIGURE_DEPENDS *.hpp)
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

FILE(GLOB HEADERS_IMPL CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../baseImpl/*.hpp)
FILE(GLOB SOURCES_IMPL CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../baseImpl/*.cpp)

FILE(GLOB HEADERS_ELEMGRAPH CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../baseImpl/elementGraph/*.hpp)
FILE(GLOB SOURCES_ELEMGRAPH CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../baseImpl/elementGraph/*.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_LIBRARY(stk_mesh_base
    NOINSTALLHEADERS ${HEADERS} ${HEADERS_IMPL} ${HEADERS_ELEMGRAPH}
    SOURCES ${SOURCES} ${SOURCES_IMPL} ${SOURCES_ELEMGRAPH}
  )
else()
  find_package(Shards REQUIRED)
  add_library(stk_mesh_base ${SOURCES} ${SOURCES_IMPL} ${SOURCES_ELEMGRAPH})
  target_link_libraries(stk_mesh_base ${Shards_LIBRARIES})
  target_link_libraries(stk_mesh_base stk_topology)
  target_link_libraries(stk_mesh_base stk_util_diag)
  target_link_libraries(stk_mesh_base stk_util_env)
  target_link_libraries(stk_mesh_base stk_util_parallel)

  find_package(BLAS REQUIRED)
  if(TARGET BLAS::BLAS)
    target_link_libraries(stk_util_util PUBLIC BLAS::BLAS)
  else()
    target_link_libraries(stk_util_util PUBLIC ${BLAS_LIBRARIES})
  endif()

endif()

target_include_directories(stk_mesh_base PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_mesh>
  $<INSTALL_INTERFACE:include>
)

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_mesh/base/)

INSTALL(FILES ${HEADERS_IMPL} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_mesh/baseImpl/)

INSTALL(FILES ${HEADERS_ELEMGRAPH} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_mesh/baseImpl/elementGraph)

if(NOT HAVE_STK_Trilinos)
  INSTALL(TARGETS stk_mesh_base EXPORT stkTargets DESTINATION ${STK_INSTALL_LIBDIR})
endif()
