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
# Decide whether to enable KokkosClassic::DoNotUse::OpenMPNode (which
# is DEPRECATED and must be enabled explicitly).  This requires that
# OpenMP be enabled in Trilinos.
#
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_OpenMP)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMPNode
  HAVE_TPETRACLASSIC_OPENMP
  "Enable the KokkosClassic::DoNotUse::OpenMPNode Node type in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  DO NOT SET THIS!  This Node type was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set Trilinos_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest."
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::OpenMPNode class was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set Trilinos_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest.")
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

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_TPINode
  HAVE_TPETRACLASSIC_THREADPOOL
  "Enable KokkosClassic::DoNotUse::TPINode in ${PACKAGE_NAME}.  DO NOT USE THIS!  This Node type was deprecated long aga and has since been removed.  Please use Kokkos::Compat::KokkosThreadsWrapperNode instead."  
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_TPINode)
  MESSAGE(FATAL_ERROR "KokkosClassic::DoNotUse::TPINode was deprecated long ago and has been removed.")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TBB)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_TBBNode
  HAVE_TPETRACLASSIC_TBB
  "Enable KokkosClassic::DoNotUse::TBBNode in ${PACKAGE_NAME}.  This Node type and this CMake option has been DEPRECATED."
  OFF
  )
#IF (${PACKAGE_NAME}_ENABLE_TBBNode)
#  MESSAGE(FATAL_ERROR "KokkosClassic::DoNotUse::TBBNode was deprecated long ago and has been removed.")
#ENDIF()

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
