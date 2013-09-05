
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# SET(HEADERS ${HEADERS}
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#   )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  gmres.h
  gmres_tools.h
  ProjectionPreconditioner.hpp
  )

SET(SOURCES ${SOURCES}
  dummy.cpp
  )

TRIBITS_ADD_LIBRARY(
  iqr
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
