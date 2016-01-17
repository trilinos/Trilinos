

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
    Kokkos_NodeAPIConfigDefs.hpp
  )
APPEND_SET(SOURCES
    Kokkos_DefaultNode.cpp
  )

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

# mfh 04 Oct 2015: Deprecate KokkosClassic_DefaultNode in favor of
# Tpetra_DefaultNode (or better yet, Tpetra_DefaultExecutionSpace).
IF(NOT DEFINED Tpetra_DefaultNode)
  IF(DEFINED KokkosClassic_DefaultNode AND NOT DEFINED Tpetra_DefaultNode)
    SET(Tpetra_DefaultNode ${KokkosClassic_DefaultNode})
  ENDIF()
ENDIF()

IF(DEFINED Tpetra_DefaultNode)
  MESSAGE(STATUS "User-set default Node: ${Tpetra_DefaultNode}")
ENDIF()

IF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
  SET(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosOpenMPWrapperNode")
  SET (HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosThreadsWrapperNode")
  SET (HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosSerialWrapperNode")
  SET (HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
ELSE()
  IF (Tpetra_INST_CUDA)
    IF(NOT Kokkos_ENABLE_Cuda)
      MESSAGE(ERROR "Tpetra_INST_CUDA is ON, but Kokkos_ENABLE_Cuda is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_OPENMP)
    IF(NOT Kokkos_ENABLE_OpenMP)
      MESSAGE(ERROR "Tpetra_INST_OPENMP is ON, but Kokkos_ENABLE_OpenMP is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_PTHREAD)
    IF(NOT Kokkos_ENABLE_Pthread)
      MESSAGE(ERROR "Tpetra_INST_PTHREAD is ON, but Kokkos_ENABLE_Pthread is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Kokkos_ENABLE_Cuda)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_OpenMP)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Serial)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosSerialWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosSerialWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
  ELSEIF (Kokkos_ENABLE_Pthread)
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
  ELSE ()
    MESSAGE(ERROR "Unrecognized or disabled Node type \"${Tpetra_DefaultNode}\".  None of the possibilities Kokkos::Compat::KokkosSerialWrapperNode, Kokkos::Compat::KokkosOpenMPWrapperNode, Kokkos::Compat::KokkosThreadsWrapperNode, or Kokkos::Compat::KokkosCudaWrapperNode are enabled, so we do not know what default Node type to pick.")
  ENDIF ()
ENDIF()
SET(Tpetra_DefaultNode "${Tpetra_DefaultNode}" CACHE STRING "" FORCE)
# mfh 04 Oct 2015: Set this for backwards compatibility ONLY.
SET(KokkosClassic_DefaultNode "${Tpetra_DefaultNode}" CACHE STRING "" FORCE)

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
