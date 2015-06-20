
#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#


#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#SET(HEADERS ${HEADERS}
#  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
#  )

#APPEND_SET(HEADERS MeshingGenie_2d.h)
APPEND_SET(SOURCES
	MeshingGenie_2d.cpp
	MeshingGenie_3d.cpp
	)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  meshinggenie
#  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
