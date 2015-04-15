#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Classic)

#
# Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

# If the "Kokkos refactor" version of Tpetra is enabled,
# disable ALL the "classic" Node types by default.
ASSERT_DEFINED(Tpetra_ENABLE_Kokkos_Refactor)

#
# Decide whether to enable KokkosClassic::OpenMPNode (which is
# DEPRECATED and must be enabled explicitly).  This requires that
# OpenMP be enabled in Trilinos.
#
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_OpenMP)
GLOBAL_SET(${PACKAGE_NAME}_ENABLE_OpenMPNode_DEFAULT OFF)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMPNode
  HAVE_TPETRACLASSIC_OPENMP
  "Enable the KokkosClassic::OpenMPNode Node type in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  This Node type is DEPRECATED and OFF by default."
  ${${PACKAGE_NAME}_ENABLE_OpenMPNode_DEFAULT}
  )
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_OpenMPNode)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_OpenMPNode is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode AND NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  MESSAGE(ERROR "${PROJECT_NAME}_ENABLE_OpenMP is OFF, so it is not legal to set ${PACKAGE_NAME}_ENABLE_OpenMPNode to ON.")
ENDIF ()

#
# Decide whether to enable KokkosClassic::DoNotUseSerialNode.  It is
# DEPRECATED and will go away ANY TIME.  It will also be disabled by
# default if the Kokkos refactor version of Tpetra is enabled.
#
IF(NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_SerialNode_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_SerialNode_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_SerialNode
  HAVE_TPETRACLASSIC_SERIAL
  "Add KokkosClassic::DoNotUse::SerialNode to the list of enabled Node types.  (DO NOT USE THIS!)"
  ${${PACKAGE_NAME}_ENABLE_SerialNode_DEFAULT}
  )
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_SerialNode)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_SerialNode is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()

#
# Users may specify TpetraClassic_ENABLE_TPINode to enable
# KokkosClassic::TPINode explicitly.  It will be DEPRECATED in the
# next Trilinos release.  This Node will also be disabled by default
# if the Kokkos refactor version of Tpetra is enabled.
#
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_ThreadPool)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThreadPool)
IF(${PROJECT_NAME}_ENABLE_ThreadPool AND ${PACKAGE_NAME}_ENABLE_ThreadPool AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TPINode_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TPINode_DEFAULT OFF)
ENDIF()
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TPINode_DEFAULT)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_TPINode_DEFAULT is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_TPINode
  HAVE_TPETRACLASSIC_THREADPOOL
  "Enable KokkosClassic::TPINode in ${PACKAGE_NAME}.  This Node type and this CMake option will be DEPRECATED in the next Trilinos release.  Please use Kokkos::Compat::KokkosThreadsWrapperNode instead."  
  ${${PACKAGE_NAME}_ENABLE_TPINode_DEFAULT}
  )
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TPINode)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_TPINode is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()


#
# Users may specify TpetraClassic_ENABLE_TBBNode to enable
# KokkosClassic::TBBNode explicitly.  It will be DEPRECATED in the
# next Trilinos release.  This Node will also be disabled by default
# if the Kokkos refactor version of Tpetra is enabled.
#
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TBB)
IF(${PACKAGE_NAME}_ENABLE_TBB AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TBBNode_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TBBNode_DEFAULT OFF)
ENDIF()
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TBBNode_DEFAULT)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_TBBNode_DEFAULT is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_TBBNode
  HAVE_TPETRACLASSIC_TBB
  "Enable KokkosClassic::TBBNode in ${PACKAGE_NAME}.  This Node type and this CMake option will be DEPRECATED in the next Trilinos release."
  ${${PACKAGE_NAME}_ENABLE_TBBNode_DEFAULT}
  )
IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TBBNode)
  MESSAGE(ERROR "${PACKAGE_NAME}_ENABLE_TBBNode is not DEFINED.  This should never happen.  Please report this bug to the Tpetra developers.")
ENDIF ()


TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_TREAT_SERIALNODE_AS_DEVICE
  HAVE_KOKKOSCLASSIC_TREAT_SERIALNODE_AS_DEVICE
  "Treat KokkosClassic::DoNotUse::SerialNode as a device node instead of a host node, for kernels where there is a difference."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
  )

#
# Add libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(NodeAPI)
ADD_SUBDIRECTORY(LinAlg)

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
