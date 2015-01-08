

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

# mfh 08 Jan 2015: Any header (.hpp) or source (.cpp) file that
# relates to a specific Node type, should have its entire contents
# protected by the macro that controls whether that Node type is
# enabled.  This makes it a bit easier to manage the build system.

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

IF (${PACKAGE_NAME}_ENABLE_Serial)
  MESSAGE(STATUS "Enabling node: Serial")
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::SerialNode")
  set(kokkos_defaultnode_default "KokkosClassic::SerialNode")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thrust)
IF (${PACKAGE_NAME}_ENABLE_Thrust)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::ThrustGPUNode")
  MESSAGE(STATUS "Enabling node: CUDA/Thrust")
  APPEND_SET(HEADERS ${CUDA_HEADERS})
  APPEND_SET(SOURCES ${CUDA_SOURCES})
ENDIF()

IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThreadPool)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_ThreadPool)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TPINode")
  MESSAGE(STATUS "Enabling node: ThreadPool")
  APPEND_SET(HEADERS ${TPI_HEADERS})
  APPEND_SET(SOURCES ${TPI_SOURCES})
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_OpenMP)
IF (${PACKAGE_NAME}_ENABLE_OpenMP)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::OpenMPNode")
  MESSAGE(STATUS "Enabling node: OpenMP")
  APPEND_SET(HEADERS ${OMP_HEADERS})
  APPEND_SET(SOURCES ${OMP_SOURCES})
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TBB)
IF (${PACKAGE_NAME}_ENABLE_TBB)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TBBNode")
  MESSAGE(STATUS "Enabling node: TBB")
  APPEND_SET(HEADERS ${TBB_HEADERS})
  APPEND_SET(SOURCES ${TBB_SOURCES})
ENDIF()

# Check whether to add the new Kokkos wrapper Nodes to the list.
# These Nodes live in the TeuchosKokkosCompat subpackage.
#
# FIXME (mfh 17 Dec 2014) The relative path to
# KokkosCompat_ClassicNodeAPI_Wrapper.hpp is a hack.  Come up with a
# more generalized way to point to that file.
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat)
IF (${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat)
  ASSERT_DEFINED(Kokkos_ENABLE_OpenMP)
  IF (Kokkos_ENABLE_OpenMP)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosOpenMPWrapperNode")
    MESSAGE(STATUS "Enabling node: OpenMP/Wrapper")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Pthread)
  IF (Kokkos_ENABLE_Pthread)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosThreadsWrapperNode")
    MESSAGE(STATUS "Enabling node: Threads/Wrapper")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
  IF (Kokkos_ENABLE_Cuda)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosCudaWrapperNode")
    MESSAGE(STATUS "Enabling node: Cuda/Wrapper")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Serial)
  IF (Kokkos_ENABLE_Serial)
    LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosSerialWrapperNode")
    MESSAGE(STATUS "Enabling node: Serial/Wrapper")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
ENDIF() # ${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat


IF(NOT DEFINED KokkosClassic_DefaultNode)
  SET(KokkosClassic_DefaultNode ${kokkos_defaultnode_default})
ENDIF()
MESSAGE(STATUS "Default node: ${KokkosClassic_DefaultNode}")
IF    (KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TPINode"       AND ${PACKAGE_NAME}_ENABLE_ThreadPool)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TBBNode"       AND ${PACKAGE_NAME}_ENABLE_TBB)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::OpenMPNode"    AND ${PACKAGE_NAME}_ENABLE_OpenMP)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::ThrustGPUNode" AND ${PACKAGE_NAME}_ENABLE_Thrust)
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
  MESSAGE(WARNING "Unrecognized or disabled node type \"${KokkosClassic_DefaultNode}\"; setting default node to KokkosClassic::SerialNode.")
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
  tpetraclassicnodeapi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_TEST_DIRECTORIES(examples)
