
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
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# general headers, sources
APPEND_SET(HEADERS
    # standard
    Kokkos_LinAlgVersion.hpp
    # data structures
    Kokkos_CrsMatrix.hpp
    Kokkos_MultiVector.hpp
    Kokkos_Vector.hpp
    # computation
    Kokkos_DefaultArithmetic.hpp
    Kokkos_DefaultSparseMultiply.hpp
    Kokkos_DefaultSparseSolve.hpp
  )

APPEND_SET(SOURCES
    Kokkos_LinAlgVersion.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  kokkoslinalg
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS kokkos kokkosnodeapi
  )

PACKAGE_ADD_TEST_DIRECTORIES(test)
