
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

# CUDA-specific source and headers
APPEND_SET(CUDA_HEADERS
    Kokkos_CUDA_util_inline_runtime.h
    Kokkos_CUDANode.cuh
    Kokkos_CUDANode.hpp
    Kokkos_CUDANodeImpl.hpp
    Kokkos_CUDA_util_sharedmem.cuh
)
APPEND_SET(CUDA_SOURCES
    Kokkos_CUDANode.cpp
)

# TBB-specific source and headers
APPEND_SET(TBB_HEADERS
    Kokkos_TBBNode.hpp
)
APPEND_SET(TBB_SOURCES
    Kokkos_TBBNode.cpp
)

# general headers, sources
APPEND_SET(HEADERS
    Kokkos_NodeAPIVersion.hpp
    Kokkos_SerialNode.hpp
    Kokkos_StandardNodeMemoryModel.hpp
  )
APPEND_SET(SOURCES
    Kokkos_NodeAPIVersion.cpp
  )

ASSERT_DEFINED(Kokkos_ENABLE_CUDA)
IF (Kokkos_ENABLE_CUDA)
  APPEND_SET(HEADERS ${CUDA_HEADERS})
  APPEND_SET(SOURCES ${CUDA_SOURCES})
ENDIF()

ASSERT_DEFINED(Kokkos_ENABLE_TBB)
IF (Kokkos_ENABLE_TBB)
  APPEND_SET(HEADERS ${TBB_HEADERS})
  APPEND_SET(SOURCES ${TBB_SOURCES})
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

