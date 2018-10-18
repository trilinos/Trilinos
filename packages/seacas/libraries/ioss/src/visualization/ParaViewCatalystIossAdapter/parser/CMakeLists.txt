# Copyright(C) 1999-2017 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
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
#     * Neither the name of NTESS nor the names of its
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

PROJECT(catalyst_parser)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(BISON REQUIRED)
SET(BisonOutput ${CMAKE_CURRENT_BINARY_DIR}/catalystParser.tab.cpp 
                ${CMAKE_CURRENT_BINARY_DIR}/catalystParser.tab.hpp)
IF(BISON_FOUND)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${BisonOutput}
      COMMAND ${BISON_EXECUTABLE}
              -d ${CMAKE_CURRENT_SOURCE_DIR}/catalystParser.ypp
              -p catalyst
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/catalystParser.ypp
      COMMENT "Generating catalystParser.tab.cpp and catalystParser.tab.hpp"
    )
ENDIF()

FIND_PACKAGE(FLEX REQUIRED)
SET(FlexOutput ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.c)
IF(FLEX_FOUND)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${FlexOutput}
      COMMAND ${FLEX_EXECUTABLE}
              --outfile=${FlexOutput}
              ${CMAKE_CURRENT_SOURCE_DIR}/catalystParser.l
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/catalystParser.l
      COMMENT "Generating lex.yy.c"
    )
ENDIF()

SET_SOURCE_FILES_PROPERTIES(${FlexOutput} PROPERTIES GENERATED TRUE)
SET_SOURCE_FILES_PROPERTIES(${BisonOutput} PROPERTIES GENERATED TRUE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(ParaView QUIET)

ADD_LIBRARY(lexer SHARED ${FlexOutput} ${BisonOutput})

IF(ParaView_FOUND)
  INCLUDE(${PARAVIEW_USE_FILE})
  ADD_DEFINITIONS("-DUSE_VTK_JSON")
  ADD_LIBRARY(catalystparser SHARED ${BisonOutput}
              ${CMAKE_CURRENT_SOURCE_DIR}/CatalystParserInterface.cpp)
  TARGET_LINK_LIBRARIES(lexer vtkjsoncpp)
  TARGET_LINK_LIBRARIES(catalystparser vtkjsoncpp)
ELSE()
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp-master/dist/)
  ADD_LIBRARY(catalystparser SHARED ${BisonOutput} 
              ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp-master/dist/jsoncpp.cpp 
              ${CMAKE_CURRENT_SOURCE_DIR}/CatalystParserInterface.cpp)
ENDIF()

TARGET_LINK_LIBRARIES(catalystparser lexer)

ADD_EXECUTABLE(catalyst_parser_driver ${CMAKE_CURRENT_SOURCE_DIR}/driver.cpp)
TARGET_LINK_LIBRARIES(catalyst_parser_driver catalystparser)

ADD_EXECUTABLE(catalyst_parser_tester ${CMAKE_CURRENT_SOURCE_DIR}/driver_test.cpp)
TARGET_LINK_LIBRARIES(catalyst_parser_tester catalystparser)

INSTALL(TARGETS lexer LIBRARY DESTINATION lib)
INSTALL(TARGETS catalystparser LIBRARY DESTINATION lib)
INSTALL(TARGETS catalyst_parser_driver RUNTIME DESTINATION bin)
INSTALL(TARGETS catalyst_parser_tester RUNTIME DESTINATION bin)
