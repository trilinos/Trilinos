
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

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Kokkos_ConfigDefs.hpp
  #Kokkos_DefaultSparseMultiply.hpp
  #Kokkos_DefaultSparseSolve.hpp
  #Kokkos_CisMatrix.hpp
  #Kokkos_MultiVector.hpp
  #Kokkos_Vector.hpp
  Kokkos_Version.hpp
  )

APPEND_SET(SOURCES
  Kokkos_Version.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  kokkos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
