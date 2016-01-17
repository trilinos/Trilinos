

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# general headers, sources
APPEND_SET(HEADERS
    # standard
    Kokkos_LinAlgVersion.hpp
    # kernels (NOTE (mfh 16 Jan 2016): Leave in place; Stokhos might still use it)
    Kokkos_MultiVectorKernelOps.hpp
  )

APPEND_SET(SOURCES
    Kokkos_LinAlgVersion.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetraclassiclinalg
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

TRIBITS_ADD_TEST_DIRECTORIES(examples)
TRIBITS_ADD_TEST_DIRECTORIES(test)
