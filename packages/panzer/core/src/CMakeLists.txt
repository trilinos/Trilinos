#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PARENT_PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  panzer-core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
