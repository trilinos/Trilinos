

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Sacado_Random.hpp 
  )

SET(SOURCES ${SOURCES}
#  Sacado_Random.cpp
  )

#
# C) Define the targets for package's library(s)
#

#TRIBITS_ADD_LIBRARY(
#  sacadoutils
#  HEADERS ${HEADERS}
#  SOURCES ${SOURCES}
#  TESTONLY
#  )
