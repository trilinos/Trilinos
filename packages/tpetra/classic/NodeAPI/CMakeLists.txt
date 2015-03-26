

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_SerialNode)
IF (${PACKAGE_NAME}_ENABLE_SerialNode)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::SerialNode")
  MESSAGE(STATUS "Enabling Node: KokkosClassic::SerialNode (DEPRECATED)")  
  set(kokkos_defaultnode_default "KokkosClassic::SerialNode")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TPINode)
IF (${PACKAGE_NAME}_ENABLE_TPINode)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TPINode")
  MESSAGE(STATUS "Enabling Node: KokkosClassic::TPINode (DEPRECATED)")
  APPEND_SET(HEADERS ${TPI_HEADERS})
  APPEND_SET(SOURCES ${TPI_SOURCES})
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_OpenMPNode)
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::OpenMPNode")
  MESSAGE(STATUS "Enabling Node: KokkosClassic::OpenMPNode (DEPRECATED)")
  APPEND_SET(HEADERS ${OMP_HEADERS})
  APPEND_SET(SOURCES ${OMP_SOURCES})
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TBBNode)
IF (${PACKAGE_NAME}_ENABLE_TBBNode)
  LIST(APPEND kokkos_enabled_nodes "KokkosClassic::TBBNode")
  MESSAGE(STATUS "Enabling Node: KokkosClassic::TBBNode (DEPRECATED)")
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
ASSERT_DEFINED(Tpetra_ENABLE_Kokkos_Refactor)
IF (NOT ${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat is OFF, but Tpetra_ENABLE_Kokkos_Refactor is ON.  In order to build the Kokkos refactor version of Tpetra, the TeuchosKokkosCompat subpackage of Teuchos must be enabled.")
ENDIF ()

IF (${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat)
  ASSERT_DEFINED(Kokkos_ENABLE_OpenMP)
  IF (Kokkos_ENABLE_OpenMP AND Tpetra_ENABLE_Kokkos_Refactor)
    #LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosOpenMPWrapperNode")
    #MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosOpenMPWrapperNode")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Pthread)
  IF (Kokkos_ENABLE_Pthread AND Tpetra_ENABLE_Kokkos_Refactor)
    #LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosThreadsWrapperNode")
    #MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosThreadsWrapperNode")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
  IF (Kokkos_ENABLE_Cuda AND Tpetra_ENABLE_Kokkos_Refactor)
    #LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosCudaWrapperNode")
    #MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosCudaWrapperNode")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
  
  ASSERT_DEFINED(Kokkos_ENABLE_Serial)
  IF (Kokkos_ENABLE_Serial AND Tpetra_ENABLE_Kokkos_Refactor)
    #LIST(APPEND kokkos_enabled_nodes "Kokkos::Compat::KokkosSerialWrapperNode")
    #MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosSerialWrapperNode")
    #APPEND_SET(HEADERS ../../../teuchos/compat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp)
  ENDIF()
ENDIF() # ${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat


IF(NOT DEFINED KokkosClassic_DefaultNode)
  SET(KokkosClassic_DefaultNode ${kokkos_defaultnode_default})
ENDIF()
MESSAGE(STATUS "Default node: ${KokkosClassic_DefaultNode}")
IF    (KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TPINode"       AND ${PACKAGE_NAME}_ENABLE_TPINode)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::TBBNode"       AND ${PACKAGE_NAME}_ENABLE_TBBNode)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::OpenMPNode"    AND ${PACKAGE_NAME}_ENABLE_OpenMPNode)
  SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE ON CACHE BOOL "")
ELSEIF(KokkosClassic_DefaultNode STREQUAL "KokkosClassic::SerialNode"    AND ${PACKAGE_NAME}_ENABLE_SerialNode)
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
  IF (Tpetra_INST_CUDA AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_OPENMP AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_PTHREAD AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Kokkos_ENABLE_Cuda AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_OpenMP AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Serial AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosSerialWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosSerialWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Pthread AND Tpetra_ENABLE_Kokkos_Refactor)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(KokkosClassic_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (${PACKAGE_NAME}_ENABLE_SerialNode)
    MESSAGE(WARNING "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\"; setting default Node to KokkosClassic::SerialNode.")
    SET(KokkosClassic_DefaultNode "KokkosClassic::SerialNode")
    SET(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALNODE ON CACHE BOOL "")
  ELSE ()
    MESSAGE(ERROR "Unrecognized or disabled Node type \"${KokkosClassic_DefaultNode}\".  None of the possibilities Kokkos::Compat::KokkosSerialWrapperNode, Kokkos::Compat::KokkosOpenMPWrapperNode, Kokkos::Compat::KokkosThreadsWrapperNode, or KokkosClassic::SerialNode are enabled, so we do not know what default Node type to pick.")
  ENDIF ()
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
