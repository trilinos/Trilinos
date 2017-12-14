
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS_PUBLIC "")
SET(HEADERS_PRIVATE "")
SET(SOURCES "")

FILE(GLOB HEADERS_PUBLIC *.hpp)
LIST( APPEND HEADERS_PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )


# Experiment to try and get both Basker versions compiling - no help
APPEND_SET(HEADERS				
  shylubasker_decl.hpp
  shylubasker_def.hpp
  shylubasker_trilinos_decl.hpp
  )


APPEND_SET(SOURCES
  shylubasker.cpp)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
    shylu_nodebasker
    HEADERS ${HEADERS_PUBLIC}
    NOINSTALLHEADERS ${HEADERS_PRIVATE}
    SOURCES ${SOURCES}
    DEPLIBS
    )
