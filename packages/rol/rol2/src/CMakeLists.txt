include(TribitsLibraryMacros)
#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
yellow("ROL source directory:  ${ROL_SRC}")

set(ROL_FUNCTION  ${ROL_SRC}/function       )
set(ROL_OBJECTIVE ${ROL_FUNCTION}/objective )
set(ROL_OPERATOR  ${ROL_FUNCTION}/operator  )


set(DIR ${CMAKE_CURRENT_SOURCE_DIR})
yellow("ROL2 source directory: ${DIR}")


set(HEADERS "")
set(SOURCES "")

include_directories( ${DIR}
                     ${DIR}/algorithm
                     ${DIR}/algorithm/TypeU
                     ${DIR}/algorithm/TypeU/linesearch
                     ${DIR}/algorithm/TypeU/linesearch/descent
                     ${DIR}/algorithm/TypeU/trustregion
                     ${DIR}/elementwise
                     ${DIR}/function
                     ${DIR}/function/objective
                     ${DIR}/function/operator
                     ${DIR}/status
                     ${DIR}/step
                     ${DIR}/step/krylov
                     ${DIR}/step/secant
                     ${DIR}/utils
                     ${DIR}/vector
)

APPEND_GLOB( HEADERS ${DIR}/ROL2.hpp )

add_subdirectory( algorithm )
add_subdirectory( elementwise )
add_subdirectory( function )
add_subdirectory( utils )
add_subdirectory( vector )



bold_blue("List of headers")
foreach(H ${HEADERS})
  blue(${H})
endforeach()


#include_directories(${CMAKE_CURRENT_BINARY_DIR})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})


#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()


#include_directories( ${DIR}   )


#
# All sources
#
#APPEND_SET( SOURCES ${ROL2_SRC}/utils/ROL2_ParameterList.cpp )

#APPEND_SET(HEADERS
#  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#)

TRIBITS_ADD_LIBRARY( rol
                     HEADERS ${HEADERS}
                     SOURCES ${SOURCES}
                   )        
