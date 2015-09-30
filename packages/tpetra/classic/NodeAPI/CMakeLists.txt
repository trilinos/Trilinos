

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
    Kokkos_BufferMacros.hpp
    Kokkos_DefaultNode.hpp
    Kokkos_NodeHelpers.hpp
    Kokkos_NodeAPIConfigDefs.hpp
    Kokkos_NodeAPIVersion.hpp
    Kokkos_NodeTrace.hpp
    Kokkos_StandardNodeMemoryModel.hpp
  )
APPEND_SET(SOURCES
    Kokkos_DefaultNode.cpp
    Kokkos_NodeAPIVersion.cpp
  )

SET(kokkos_enabled_nodes "")

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_SerialNode)
IF (${PACKAGE_NAME}_ENABLE_SerialNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::SerialNode class was deprecated long ago and is no longer supported.")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TPINode)
IF (${PACKAGE_NAME}_ENABLE_TPINode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::TPINode class was deprecated long ago and has since been removed.")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_OpenMPNode)
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::OpenMPNode class was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set Trilinos_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest.")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TBBNode)
IF (${PACKAGE_NAME}_ENABLE_TBBNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::TBBNode class was deprecated long ago and has since been removed.")
ENDIF()

ASSERT_DEFINED(Kokkos_ENABLE_OpenMP)
ASSERT_DEFINED(Kokkos_ENABLE_Pthread)
ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
ASSERT_DEFINED(Kokkos_ENABLE_Serial)

IF(NOT DEFINED KokkosClassic_DefaultNode)
  SET(KokkosClassic_DefaultNode ${kokkos_defaultnode_default})
ENDIF()
MESSAGE(STATUS "Default node: ${KokkosClassic_DefaultNode}")
IF (KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
ELSEIF (KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosOpenMPWrapperNode")
  SET (HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
ELSEIF (KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosThreadsWrapperNode")
  SET (HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
ELSEIF (KokkosClassic_DefaultNode STREQUAL "Kokkos::Compat::KokkosSerialWrapperNode")
  SET (HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
ELSE()
  IF (Tpetra_INST_CUDA)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_OPENMP)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_PTHREAD)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Kokkos_ENABLE_Cuda)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_OpenMP)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Serial)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosSerialWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosSerialWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Pthread)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
  ELSE ()
    MESSAGE(ERROR "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\".  None of the possibilities Kokkos::Compat::KokkosSerialWrapperNode, Kokkos::Compat::KokkosOpenMPWrapperNode, Kokkos::Compat::KokkosThreadsWrapperNode, or Kokkos::Compat::KokkosCudaWrapperNode are enabled, so we do not know what default Node type to pick.")
  ENDIF ()
ENDIF()
SET(KokkosClassic_DefaultNode "${KokkosClassic_DefaultNode}" CACHE STRING "" FORCE)

SET(KokkosClassic_ENABLED_NODES ${kokkos_enabled_nodes} CACHE STRING "List of nodes enabled in KokkosClassic.")

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
