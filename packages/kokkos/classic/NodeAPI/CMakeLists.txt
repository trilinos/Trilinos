

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

# CUDA-specific source and headers
APPEND_SET(CUDA_HEADERS
    Kokkos_CUDANodeMemoryModel.hpp
    Kokkos_CUDANodeMemoryModelImpl.hpp
    Kokkos_ThrustGPUNode.cuh
    Kokkos_ThrustGPUNode.hpp
    Kokkos_ThrustGPUWrappers.hpp
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
    Kokkos_SerialNode.cpp
    Kokkos_DefaultNode.cpp
    Kokkos_NodeAPIVersion.cpp
  )

SET(kokkos_enabled_nodes "")

IF (KokkosClassic_ENABLE_Serial)
  MESSAGE(STATUS "Enabling node: Serial")
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::SerialNode")
  set(kokkos_defaultnode_default "KokkosClassic::SerialNode")
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_Thrust)
IF (KokkosClassic_ENABLE_Thrust)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::ThrustGPUNode")
  MESSAGE(STATUS "Enabling node: CUDA/Thrust")
  APPEND_SET(HEADERS ${CUDA_HEADERS})
  APPEND_SET(SOURCES ${CUDA_SOURCES})
ENDIF()

IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  ASSERT_DEFINED(KokkosClassic_ENABLE_ThreadPool)
ENDIF()
IF (KokkosClassic_ENABLE_ThreadPool)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TPINode")
  MESSAGE(STATUS "Enabling node: ThreadPool")
  APPEND_SET(HEADERS ${TPI_HEADERS})
  APPEND_SET(SOURCES ${TPI_SOURCES})
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_OpenMP)
IF (KokkosClassic_ENABLE_OpenMP)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::OpenMPNode")
  MESSAGE(STATUS "Enabling node: OpenMP")
  APPEND_SET(HEADERS ${OMP_HEADERS})
  APPEND_SET(SOURCES ${OMP_SOURCES})
ENDIF()

ASSERT_DEFINED(KokkosClassic_ENABLE_TBB)
IF (KokkosClassic_ENABLE_TBB)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TBBNode")
  MESSAGE(STATUS "Enabling node: TBB")
  APPEND_SET(HEADERS ${TBB_HEADERS})
  APPEND_SET(SOURCES ${TBB_SOURCES})
ENDIF()

# Check whether to add the new Kokkos wrapper Nodes to the list.
# These Nodes live in the KokkosCompat subpackage.
ASSERT_DEFINED(KokkosClassic_ENABLE_KokkosCompat)
IF (KokkosClassic_ENABLE_KokkosCompat)
  ASSERT_DEFINED(Kokkos_ENABLE_OpenMP)
  IF (Kokkos_ENABLE_OpenMP)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosOpenMPWrapperNode")
    MESSAGE(STATUS "Enabling node: OpenMP/Wrapper")
    APPEND_SET(HEADERS ../../compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Pthread)
  IF (Kokkos_ENABLE_Pthread)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosThreadsWrapperNode")
    MESSAGE(STATUS "Enabling node: Threads/Wrapper")
    APPEND_SET(HEADERS ../../compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
  IF (Kokkos_ENABLE_Cuda)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosCudaWrapperNode")
    MESSAGE(STATUS "Enabling node: Cuda/Wrapper")
    APPEND_SET(HEADERS ../../compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  # FIXME (mfh 02 Oct 2014) I can't find a CMake variable that tells me
  # whether the Serial device is enabled.  Is it enabled all the time?
  LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosSerialWrapperNode")
  MESSAGE(STATUS "Enabling node: Serial/Wrapper")
  APPEND_SET(HEADERS ../../compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
ENDIF() # KokkosClassic_ENABLE_KokkosCompat


IF(NOT DEFINED KokkosClassic_DefaultNode)
  SET(KokkosClassic_DefaultNode ${kokkos_defaultnode_default})
ENDIF()
MESSAGE(STATUS "Default node: ${KokkosClassic_DefaultNode}")
IF    (KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TPINode"       AND KokkosClassic_ENABLE_ThreadPool)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TBBNode"       AND KokkosClassic_ENABLE_TBB)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::OpenMPNode"    AND KokkosClassic_ENABLE_OpenMP)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::ThrustGPUNode" AND KokkosClassic_ENABLE_Thrust)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THRUSTGPUNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::SerialNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosOpenMPWrapperNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosThreadsWrapperNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosSerialWrapperNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
ELSE()
  MESSAGE(WARNING "Unrecognized or disabled node type; setting default node to KokkosClassic::SerialNode.")
  SET(KokkosClassic_DefaultNode "KokkosClassic::SerialNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALNODE ON CACHE BOOL "")
ENDIF()
SET(KokkosClassic_DefaultNode "${KokkosClassic_DefaultNode}" CACHE STRING "" FORCE)

SET(KokkosClassic_ENABLED_NODES ${kokkos_enabled_nodes} CACHE STRING "List of nodes enabled in KokkkosClassic.")

CONFIGURE_FILE(cmake/KokkosClassic_DefaultNode_config.h.in  KokkosClassic_DefaultNode_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
APPEND_SET(HEADERS  ${CMAKE_CURRENT_BINARY_DIR}/KokkosClassic_DefaultNode_config.h)

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
