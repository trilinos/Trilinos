
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

APPEND_SET(HEADERS 
  Krylov_ConfigDefs.hpp
  Krylov_Leja.hpp
  Krylov_Operator.hpp
  Krylov_OperatorTraits.hpp
  Krylov_PolyBasis.hpp
  Krylov_ProductIsZero.hpp
  Krylov_ShiftsOutOfOrder.hpp
  )

APPEND_SET(SOURCES
  Krylov_Version.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  krylovsupport
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
