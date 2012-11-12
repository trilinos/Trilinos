
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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
    Kokkos_CUDANodeMemoryModel.hpp
    Kokkos_CUDANodeMemoryModelImpl.hpp
    Kokkos_ThrustGPUNode.cuh
    Kokkos_ThrustGPUNode.hpp
    Kokkos_CUDANodeUtils.hpp
)
APPEND_SET(CUDA_SOURCES
    Kokkos_CUDANodeMemoryModel.cpp
    Kokkos_ThrustGPUNode.cpp
    Kokkos_CUDANodeUtils.cpp
)

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
APPEND_SET(TPI_SOURCES
    Kokkos_TPINode.cpp
)

# OpenMP-specific source and headers
APPEND_SET(OMP_HEADERS
    Kokkos_OpenMPNode.hpp
)
APPEND_SET(OMP_SOURCES
    Kokkos_OpenMPNode.cpp
)

# general headers, sources
APPEND_SET(HEADERS
    Kokkos_BufferMacros.hpp
    Kokkos_DefaultNode.hpp
    Kokkos_NodeHelpers.hpp
    Kokkos_NodeAPIConfigDefs.hpp
    Kokkos_NodeAPIVersion.hpp
    Kokkos_NodeTrace.hpp
    Kokkos_SerialNode.hpp
    Kokkos_StandardNodeMemoryModel.hpp
  )
APPEND_SET(SOURCES
    Kokkos_DefaultNode.cpp
    Kokkos_NodeAPIVersion.cpp
  )

MESSAGE(STATUS "Enabling node: Serial")

ASSERT_DEFINED(KokkosClassic_ENABLE_Thrust)
IF (KokkosClassic_ENABLE_Thrust)
  MESSAGE(STATUS "Enabling node: CUDA/Thrust")
  APPEND_SET(HEADERS ${CUDA_HEADERS})
  APPEND_SET(SOURCES ${CUDA_SOURCES})
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_TBB)
IF (KokkosClassic_ENABLE_TBB)
  MESSAGE(STATUS "Enabling node: TBB")
  APPEND_SET(HEADERS ${TBB_HEADERS})
  APPEND_SET(SOURCES ${TBB_SOURCES})
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_ThreadPool)
IF (KokkosClassic_ENABLE_ThreadPool)
  MESSAGE(STATUS "Enabling node: ThreadPool")
  APPEND_SET(HEADERS ${TPI_HEADERS})
  APPEND_SET(SOURCES ${TPI_SOURCES})
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_OpenMP)
IF (KokkosClassic_ENABLE_OpenMP)
  MESSAGE(STATUS "Enabling node: OpenMP")
  APPEND_SET(HEADERS ${OMP_HEADERS})
  APPEND_SET(SOURCES ${OMP_SOURCES})
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  kokkosnodeapi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS kokkos
  )

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_TEST_DIRECTORIES(examples)
