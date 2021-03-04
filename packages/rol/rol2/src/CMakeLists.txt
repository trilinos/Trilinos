include(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
yellow("ROL source directory:  ${ROL_SRC}")

set(ROL2_SRC   ${CMAKE_CURRENT_SOURCE_DIR})
yellow("ROL2 source directory:  ${ROL2_SRC}")


set(ROL2_ALGORITHM   ${ROL2_SRC}/algorithm   )
set(ROL2_ELEMENTWISE ${ROL2_SRC}/elementwise )
set(ROL2_FUNCTION    ${ROL2_SRC}/function    )
set(ROL2_UTILS       ${ROL2_SRC}/utils       )
set(ROL2_VECTOR      ${ROL2_SRC}/vector      )

#
# B) Define the header and source files (and directories)
#

#
# src
#

set(HEADERS "")
set(SOURCES "")

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

include_directories( ${ROL2_SRC}   )
APPEND_GLOB( HEADERS ${ROL2_SRC}/ROL2.hpp )

include_directories( ${ROL2_UTILS} )
APPEND_GLOB( HEADERS ${ROL2_UTILS}/*.hpp )

include_directories( ${ROL2_ELEMENTWISE} )
APPEND_GLOB( HEADERS ${ROL2_ELEMENTWISE}/*.hpp )

#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL_Utilities.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL_TypeTraits.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL2_FoldExpressions.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL2_Ptr.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL2_ParameterList.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL2_Stream_Decl.hpp )
#APPEND_GLOB( HEADERS ${ROL2_UTILS}/ROL2_Stream_Def.hpp )


bold_blue("List of headers")
foreach(H ${HEADERS})
  blue(${H})
endforeach()

#
# All sources
#
APPEND_SET( SOURCES ${ROL2_UTILS}/ROL2_ParameterList.cpp )

#APPEND_SET(HEADERS
#  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#)

TRIBITS_ADD_LIBRARY( rol
                     HEADERS ${HEADERS}
                     SOURCES ${SOURCES}
                   )        
