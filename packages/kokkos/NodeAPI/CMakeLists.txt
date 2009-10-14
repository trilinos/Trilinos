
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

# TBB-specific source and headers
APPEND_SET(TBB_HEADERS
    Kokkos_TBBNode.hpp
)
APPEND_SET(TBB_SOURCES
    Kokkos_TBBNode.cpp
)

# TPI-specific source and headers
APPEND_SET(TPI_HEADERS
    Kokkos_TPINode.hpp
)

# general headers, sources
APPEND_SET(HEADERS
    Kokkos_DefaultNode.hpp
    Kokkos_NodeHelpers.hpp
    Kokkos_NodeAPIConfigDefs.hpp
    Kokkos_NodeAPIVersion.hpp
    Kokkos_SerialNode.hpp
    Kokkos_StandardNodeMemoryModel.hpp
  )
APPEND_SET(SOURCES
    Kokkos_DefaultNode.cpp
    Kokkos_NodeAPIVersion.cpp
  )

ASSERT_DEFINED(Kokkos_ENABLE_TBB)
IF (Kokkos_ENABLE_TBB)
  APPEND_SET(HEADERS ${TBB_HEADERS})
  APPEND_SET(SOURCES ${TBB_SOURCES})
ENDIF()

ASSERT_DEFINED(Kokkos_ENABLE_ThreadPool)
IF (Kokkos_ENABLE_ThreadPool)
  APPEND_SET(HEADERS ${TPI_HEADERS})
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  kokkosnodeapi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS kokkos
  )

PACKAGE_ADD_TEST_DIRECTORIES(test)
