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

IF(NOT ${PROJECT_NAME}_ENABLE_STKIO)

LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/GearsFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Hex20Fixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Hex27Fixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/PyramidFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/QuadShellFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Tet10Fixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/TetFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/TriFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/WedgeFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/degenerate_mesh.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/heterogeneous_mesh.cpp")

ENDIF()

if(HAVE_STK_Trilinos)
  TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE}_SOURCE_DIR/stk_mesh)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

  TRIBITS_ADD_LIBRARY(stk_mesh_fixtures
    NOINSTALLHEADERS ${HEADERS}
    SOURCES ${SOURCES}
    ${Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG}
  )
else()
  add_library(stk_mesh_fixtures ${SOURCES})

  target_include_directories(stk_mesh_fixtures PUBLIC
    $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_unit_test_utils>
    $<INSTALL_INTERFACE:include>
  )

  target_link_libraries(stk_mesh_fixtures PUBLIC stk_unit_test_utils)
  if(STK_ENABLE_STKIO)
    target_link_libraries(stk_mesh_fixtures PUBLIC stk_io)
  endif()
  target_link_libraries(stk_mesh_fixtures PUBLIC stk_mesh_base)
endif()

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_unit_test_utils/stk_mesh_fixtures)

if(NOT HAVE_STK_Trilinos)
  INSTALL(TARGETS stk_mesh_fixtures EXPORT stkTargets DESTINATION ${STK_INSTALL_LIBDIR})
endif()
