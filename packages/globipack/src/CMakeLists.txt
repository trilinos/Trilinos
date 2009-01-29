
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  ${PACKAGE_NAME}_ArmijoPolyInterpLineSearch_decl.hpp
  ${PACKAGE_NAME}_ArmijoPolyInterpLineSearch_def.hpp
  ${PACKAGE_NAME}_ArmijoPolyInterpLineSearch.hpp
  ${PACKAGE_NAME}_ConfigDefs.hpp
  ${PACKAGE_NAME}_LineSearchBase.hpp
  ${PACKAGE_NAME}_MeritFunc1DBase.hpp
  ${PACKAGE_NAME}_Types.hpp
  ${PACKAGE_NAME}_Version.hpp
  )

SET(SOURCES ${SOURCES}
  ${PACKAGE_NAME}_ArmijoPolyInterpLineSearch.cpp
  ${PACKAGE_NAME}_Version.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  globipack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
