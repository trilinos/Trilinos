# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

FIND_PACKAGE(BISON 3.0 REQUIRED)
SET(BisonOutput ${CMAKE_CURRENT_BINARY_DIR}/phactoriParser.tab.cpp
                ${CMAKE_CURRENT_BINARY_DIR}/phactoriParser.tab.hpp)
IF(BISON_FOUND)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${BisonOutput}
      COMMAND ${BISON_EXECUTABLE}
              -d ${CMAKE_CURRENT_SOURCE_DIR}/phactoriParser.ypp
              -p catalyst
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/phactoriParser.ypp
      COMMENT "Generating phactoriParser.tab.cpp and phactoriParser.tab.hpp"
    )
ENDIF()

FIND_PACKAGE(FLEX 2.5 REQUIRED)
SET(FlexOutput ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.c)
IF(FLEX_FOUND)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${FlexOutput}
      COMMAND ${FLEX_EXECUTABLE}
              --outfile=${FlexOutput}
              ${CMAKE_CURRENT_SOURCE_DIR}/phactoriParser.l
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/phactoriParser.l
      COMMENT "Generating lex.yy.c"
    )
ENDIF()

SET_SOURCE_FILES_PROPERTIES(${FlexOutput} PROPERTIES GENERATED TRUE)
SET_SOURCE_FILES_PROPERTIES(${BisonOutput} PROPERTIES GENERATED TRUE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(ParaView QUIET)

ADD_LIBRARY(phactori_lexer SHARED ${FlexOutput} ${BisonOutput})
TARGET_LINK_LIBRARIES(phactori_lexer
    PRIVATE
      VTK::jsoncpp)

ADD_LIBRARY(phactori_parser SHARED ${BisonOutput}
              ${CMAKE_CURRENT_SOURCE_DIR}/PhactoriParserInterface.cxx)

TARGET_LINK_LIBRARIES(phactori_parser
    PRIVATE
      phactori_lexer
      VTK::jsoncpp)

ADD_EXECUTABLE(parse_phactori ${CMAKE_CURRENT_SOURCE_DIR}/parse_phactori.cxx)
TARGET_LINK_LIBRARIES(parse_phactori
    PRIVATE
      phactori_parser
      VTK::jsoncpp)

INSTALL(TARGETS phactori_lexer LIBRARY DESTINATION lib)
INSTALL(TARGETS phactori_parser LIBRARY DESTINATION lib)
INSTALL(TARGETS parse_phactori RUNTIME DESTINATION bin)
