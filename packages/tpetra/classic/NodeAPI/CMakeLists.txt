

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_OpenMPNode)
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::OpenMPNode class was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set ${PROJECT_NAME}_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest.")
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
IF(DEFINED Tpetra_DefaultNode)
  MESSAGE(STATUS "User-set default Node: ${Tpetra_DefaultNode}")
ELSE()
  IF(DEFINED KokkosClassic_DefaultNode AND NOT DEFINED Tpetra_DefaultNode)
    MESSAGE(WARNING "You set Tpetra's default Node type to ${KokkosClassic_DefaultNode}, but you used KokkosClassic_DefaultNode, which has been deprecated.  Please use Tpetra_DefaultNode instead.")
    SET(Tpetra_DefaultNode ${KokkosClassic_DefaultNode})
  ELSE()
    SET(Tpetra_DefaultNode "")
  ENDIF()
ENDIF()

# If the user set a default Node type, Tpetra_DefaultNode will 
# be a nonempty string.  In that case, check that the requested
# default Node type is actually enabled before using it.  If the
# user did not set a default Node type, set one ourselves.  Note
# precedence rule: First Cuda, then OpenMP, then Serial, then 
# Threads (least preferred because it will likely interfere with
# performance of users' own threads -- OpenMP doesn't have
# that problem).

IF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
  IF (NOT Kokkos_ENABLE_Cuda)
    MESSAGE (FATAL_ERROR "Kokkos_ENABLE_Cuda is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosCudaWrapperNode.")
  ENDIF ()
  IF (NOT Tpetra_INST_CUDA)
    MESSAGE (FATAL_ERROR "Tpetra_INST_CUDA is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosCudaWrapperNode.")
  ENDIF ()
  SET (HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosOpenMPWrapperNode")
  IF (NOT Kokkos_ENABLE_OpenMP)
    MESSAGE (FATAL_ERROR "Kokkos_ENABLE_OpenMP is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosOpenMPWrapperNode.")
  ENDIF ()
  IF (NOT Tpetra_INST_OPENMP)
    MESSAGE (FATAL_ERROR "Tpetra_INST_OPENMP is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosOpenMPWrapperNode.")
  ENDIF ()
  SET (HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosThreadsWrapperNode")
  IF (NOT Kokkos_ENABLE_Pthread)
    MESSAGE (FATAL_ERROR "Kokkos_ENABLE_Pthread is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosThreadsWrapperNode.")
  ENDIF ()
  IF (NOT Tpetra_INST_PTHREAD)
    MESSAGE (FATAL_ERROR "Tpetra_INST_PTHREAD is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosThreadsWrapperNode.")
  ENDIF ()
  SET (HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosSerialWrapperNode")
  IF (NOT Kokkos_ENABLE_Serial)
    MESSAGE (FATAL_ERROR "Kokkos_ENABLE_Serial is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosSerialWrapperNode.")
  ENDIF ()
  IF (NOT Tpetra_INST_SERIAL)
    MESSAGE (FATAL_ERROR "Tpetra_INST_SERIAL is OFF, but you set Tpetra's default Node type to Kokkos::Compat::KokkosSerialWrapperNode.")
  ENDIF ()
  SET (HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")
ELSEIF (Tpetra_DefaultNode STREQUAL "")
  #
  # The user didn't set Tpetra_DefaultNode, so we need to pick it.
  #
  IF (Tpetra_INST_CUDA)
    IF(NOT Kokkos_ENABLE_Cuda)
      MESSAGE(FATAL_ERROR "Tpetra_INST_CUDA is ON, but Kokkos_ENABLE_Cuda is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosCudaWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosCudaWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_OPENMP)
    IF(NOT Kokkos_ENABLE_OpenMP)
      MESSAGE(FATAL_ERROR "Tpetra_INST_OPENMP is ON, but Kokkos_ENABLE_OpenMP is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosOpenMPWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosOpenMPWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_SERIAL)
    IF(NOT Kokkos_ENABLE_Serial)
      MESSAGE(FATAL_ERROR "Tpetra_INST_SERIAL is ON, but Kokkos_ENABLE_Serial is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosSerialWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosSerialWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE ON CACHE BOOL "")  
  ELSEIF (Tpetra_INST_PTHREAD)
    IF(NOT Kokkos_ENABLE_Pthread)
      MESSAGE(FATAL_ERROR "Tpetra_INST_PTHREAD is ON, but Kokkos_ENABLE_Pthread is OFF.")
    ENDIF()
    MESSAGE(STATUS "Setting default Node to Kokkos::Compat::KokkosThreadsWrapperNode.")
    SET(Tpetra_DefaultNode "Kokkos::Compat::KokkosThreadsWrapperNode")
    SET(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE ON CACHE BOOL "")  
  ELSE ()
    # Extra check for the fix for Issue #190.
    MESSAGE(FATAL_ERROR "Tpetra cannot pick a default Node type if all Node types are disabled.")
  ENDIF ()
ELSE ()
  #
  # The user set Tpetra_DefaultNode to something other than 
  # the empty string, but we don't recognize their string.
  #
  MESSAGE(FATAL_ERROR "Unrecognized or disabled Node type \"${Tpetra_DefaultNode}\".  None of the possibilities Kokkos::Compat::KokkosSerialWrapperNode, Kokkos::Compat::KokkosOpenMPWrapperNode, Kokkos::Compat::KokkosThreadsWrapperNode, or Kokkos::Compat::KokkosCudaWrapperNode are enabled, so we do not know what default Node type to pick.")
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
