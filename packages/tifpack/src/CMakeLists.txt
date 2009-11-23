
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Tifpack_ConfigDefs.hpp
  Tifpack_Version.hpp
  Tifpack_IlukGraph.hpp
  Tifpack_Parameters.hpp
  Tifpack_PointRelaxation.hpp
  Tifpack_Preconditioner.hpp
  Tifpack_Version.hpp
  )

SET(SOURCES ${SOURCES}
  Tifpack_Parameters.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tifpack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

