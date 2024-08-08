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

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)
FILE(GLOB PRED_HEADERS predicates/*.hpp)
FILE(GLOB PRED_SOURCES predicates/*.cpp)

message("STKMiddle_mesh requires CDT. (use -DSTK_ENABLE_STKMiddle_mesh=OFF if you don't need it and you don't have CDT.)")
find_package(CDT REQUIRED)
if(CDT_FOUND)
  message(STATUS "Found CDT library ${CDT_DIR}")
else()
  message(FATAL_ERROR "Failed to find CDT package")
endif()

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_LIBRARY(stk_middle_mesh
    NOINSTALLHEADERS ${HEADERS} ${PRED_HEADERS}
    SOURCES ${SOURCES} ${PRED_SOURCES}
  )
else()
  add_library(stk_middle_mesh ${SOURCES} ${PRED_SOURCES})
  target_link_libraries(stk_middle_mesh PUBLIC stk_util_parallel)
  target_link_libraries(stk_middle_mesh PUBLIC stk_search)
  find_package(CBLAS REQUIRED)
  find_package(Lapack REQUIRED)
  target_link_libraries(stk_middle_mesh PUBLIC lapack)
endif()

target_include_directories(stk_middle_mesh PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_middle_mesh>
  $<INSTALL_INTERFACE:include>
)

target_include_directories(stk_middle_mesh PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_util>
  $<INSTALL_INTERFACE:include>
)

if(CDT_INCLUDE_DIRS)
SET(CDT_INC_DIR "${CDT_INCLUDE_DIRS}")
else()
message(STATUS "CDT_INCLUDE_DIRS not set, trying get_target_property...")
get_target_property(CDT_INC_DIR CDT::CDT INTERFACE_INCLUDE_DIRECTORIES)
endif()

target_include_directories(stk_middle_mesh PUBLIC ${CDT_INC_DIR})

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_middle_mesh/)

INSTALL(FILES ${PRED_HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_middle_mesh/predicates)

if(NOT HAVE_STK_Trilinos)
  INSTALL(TARGETS stk_middle_mesh EXPORT stkTargets DESTINATION ${STK_INSTALL_LIBDIR})
endif()
