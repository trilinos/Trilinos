include(TribitsLibraryMacros)
#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
yellow("ROL source directory:  ${ROL_SRC}")

set(ROL_FUNCTION  ${ROL_SRC}/function       )
set(ROL_OBJECTIVE ${ROL_FUNCTION}/objective )
set(ROL_OPERATOR  ${ROL_FUNCTION}/operator  )




set(ROL2_SRC   ${CMAKE_CURRENT_SOURCE_DIR})
yellow("ROL2 source directory:  ${ROL2_SRC}")

set(ROL2_VECTOR      ${ROL2_SRC}/vector      )
set(ROL2_ALGORITHM   ${ROL2_SRC}/algorithm   )
set(ROL2_ELEMENTWISE ${ROL2_SRC}/elementwise )
set(ROL2_FUNCTION    ${ROL2_SRC}/function    )

set(ROL2_OBJECTIVE   ${ROL2_FUNCTION}/objective )
set(ROL2_OPERATIOR   ${ROL2_FUNCTION}/operator  )

set(ROL2_UTILS       ${ROL2_SRC}/utils       )

set(HEADERS "")
set(SOURCES "")

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


include_directories( ${ROL2_SRC}   )
APPEND_GLOB( HEADERS ${ROL2_SRC}/ROL2.hpp )

include_directories( ${ROL2_UTILS} )
APPEND_GLOB( HEADERS ${ROL2_UTILS}/*.hpp )

include_directories( ${ROL2_ELEMENTWISE} )
APPEND_GLOB( HEADERS ${ROL2_ELEMENTWISE}/*.hpp )

include_directories( ${ROL2_VECTOR} )
APPEND_GLOB( HEADERS ${ROL2_VECTOR}/*.hpp )

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
