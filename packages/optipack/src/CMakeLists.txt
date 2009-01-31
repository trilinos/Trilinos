
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
  ${PACKAGE_NAME}_ConfigDefs.hpp
  ${PACKAGE_NAME}_NonlinearCG.hpp 
  ${PACKAGE_NAME}_NonlinearCG_decl.hpp 
  ${PACKAGE_NAME}_NonlinearCG_def.hpp 
  ${PACKAGE_NAME}_UnconstrainedOptMeritFunc1D.hpp 
  ${PACKAGE_NAME}_UnconstrainedOptMeritFunc1D_decl.hpp 
  ${PACKAGE_NAME}_UnconstrainedOptMeritFunc1D_def.hpp 
  ${PACKAGE_NAME}_Version.hpp
  )

SET(SOURCES ${SOURCES}
  ${PACKAGE_NAME}_NonlinearCG.cpp 
  ${PACKAGE_NAME}_UnconstrainedOptMeritFunc1D.cpp 
  ${PACKAGE_NAME}_Version.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  optipack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
