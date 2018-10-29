#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Classic)

#
# Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

ASSERT_DEFINED(Tpetra_ENABLE_Kokkos_Refactor)
IF (NOT Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(FATAL_ERROR "Tpetra_ENABLE_Kokkos_Refactor MUST be set to ON, if you set it at all.  You don't need to set this variable any more!")
ENDIF()

#
# Decide whether to enable KokkosClassic::DoNotUse::OpenMPNode (which
# is DEPRECATED and must be enabled explicitly).  This requires that
# OpenMP be enabled in Trilinos.
#
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_OpenMP)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMPNode
  HAVE_TPETRACLASSIC_OPENMP
  "Enable the KokkosClassic::DoNotUse::OpenMPNode Node type in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  DO NOT SET THIS!  This Node type was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set ${PROJECT_NAME}_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest."
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_OpenMPNode)
  MESSAGE(FATAL_ERROR "The KokkosClassic::DoNotUse::OpenMPNode class was deprecated long ago and has since been removed.  If you want to use OpenMP with Tpetra, just set ${PROJECT_NAME}_ENABLE_OpenMP:BOOL=ON and Tpetra will do the rest.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_SerialNode
  HAVE_TPETRACLASSIC_SERIAL
  "Add KokkosClassic::DoNotUse::SerialNode to the list of enabled Node types.  (DO NOT USE THIS!)  Please use Kokkos::Compat::KokkosSerialWrapperNode instead."
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_SerialNode)
  MESSAGE(FATAL_ERROR "${PACKAGE_NAME}_ENABLE_SerialNode is no longer available.  Please use Kokkos::Compat::KokkosSerialWrapperNode instead.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_TBBNode
  HAVE_TPETRACLASSIC_TBB
  "Enable KokkosClassic::DoNotUse::TBBNode in ${PACKAGE_NAME}.  This Node type and this CMake option was deprecated long ago and has since been REMOVED."
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_TBBNode)
  MESSAGE(FATAL_ERROR "KokkosClassic::DoNotUse::TBBNode was deprecated long ago and has since been removed.")
ENDIF()

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
