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
# Decide whether to enable KokkosClassic::OpenMPNode
# (which is DEPRECATED).
#
# Users may specify KokkosClassic_ENABLE_OpenMP to enable
# the DEPRECATED KokkosClassic::OpenMPNode explicitly.
# Even still, it won't be enabled unless OpenMP is enabled in
# Trilinos.
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_OpenMP)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_OpenMP
  HAVE_TPETRACLASSIC_OPENMP
  "Enable KokkosClassic::OpenMPNode in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  This option is DEPRECATED and is OFF by default.  It will be removed in future Trilinos releases.  If the Kokkos refactor version of Tpetra is enabled (it will be enabled by default in the next Trilinos release), you may use Kokkos::Compat::KokkosOpenMPWrapperNode instead of KokkosClassic::OpenMPNode."
  OFF
  )
# Set whether to enable KokkosClassic::OpenMPNode by default.
# The default is OFF, since KokkosClassic_ENABLE_OpenMP is OFF
# by default.  We also disable this Node type by default if the
# Kokkos refactor version of Tpetra is enabled.
IF(KokkosClassic_ENABLE_OpenMP AND ${PROJECT_NAME}_ENABLE_OpenMP AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT OFF)
ENDIF()

# This is the final CMake option that controls whether
# KokkosClassic::OpenMPNode is enabled.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMP
  HAVE_TPETRACLASSIC_OPENMP
  "Enable the KokkosClassic::OpenMPNode Node type in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  This Node type is DEPRECATED and OFF by default."
  ${${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT}
  )
IF(${PACKAGE_NAME}_ENABLE_OpenMP AND NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  MESSAGE(ERROR "${PROJECT_NAME}_ENABLE_OpenMP is OFF, so it is not legal to set ${PACKAGE_NAME}_ENABLE_OpenMP to ON.")
ENDIF()

#
# Decide whether to enable KokkosClassic::SerialNode.
# It will be DEPRECATED in the next Trilinos release.
#
# Users may specify KokkosClassic_ENABLE_Serial to enable
# the soon-to-be deprecated KokkosClassic::SerialNode explicitly.
# It will still be disabled by default if the Kokkos refactor
# version of Tpetra is enabled.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_Serial
  HAVE_TPETRACLASSIC_SERIAL
  "Enable KokkosClassic::SerialNode in ${PACKAGE_NAME}."
  ON
  )
# Disable KokkosClassic::SerialNode by default, if the Kokkos
# refactor version of Tpetra is enabled.
IF(KokkosClassic_ENABLE_Serial AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_Serial_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_Serial_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_Serial
  HAVE_TPETRACLASSIC_SERIAL
  "Add KokkosClassic::SerialNode to the list of enabled Node types."
  ${${PACKAGE_NAME}_ENABLE_Serial_DEFAULT}
  )

#
# Decide whether to enable KokkosClassic::TPINode.
# It will be DEPRECATED in the next Trilinos release.
#
# Users may specify KokkosClassic_ENABLE_ThreadPool to enable
# the soon-to-be deprecated KokkosClassic::TPINode explicitly.
# It will still be disabled by default if the Kokkos refactor
# version of Tpetra is enabled.
IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThreadPool)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_ThreadPool
  HAVE_TPETRACLASSIC_THREADPOOL
  "Enable KokkosClassic::TPINode in ${PACKAGE_NAME}.  This Node type and this CMake option will be DEPRECATED in the next Trilinos release.  Please use Kokkos::Compat::KokkosThreadsWrapperNode instead."
  ${${PROJECT_NAME}_ENABLE_ThreadPool}
  )
# Disable KokkosClassic::TPINode by default, if the Kokkos
# refactor version of Tpetra is enabled.
IF(KokkosClassic_ENABLE_ThreadPool AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_ThreadPool_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_ThreadPool_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_ThreadPool
  HAVE_TPETRACLASSIC_THREADPOOL
  "Add KokkosClassic::TPINode to the list of enabled Node types."
  ${${PACKAGE_NAME}_ENABLE_ThreadPool_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_CUDA_FLOAT
  HAVE_TPETRACLASSIC_CUDA_FLOAT
  "Enable float type for CUDA nodes."
  ON
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_CUDA_DOUBLE
  HAVE_TPETRACLASSIC_CUDA_DOUBLE
  "Enable double type for CUDA nodes."
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_CUDA_NODE_MEMORY_PROFILING
  HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
  "Enable profiling of memory transfers for CUDA nodes."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_CUDA_NODE_MEMORY_TRACE
  HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_TRACE
  "Enable memory movement trace for CUDA nodes."
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_TREAT_SERIALNODE_AS_DEVICE
  HAVE_KOKKOSCLASSIC_TREAT_SERIALNODE_AS_DEVICE
  "Treat KokkosClassic::SerialNode as a device node instead of a host node, for kernels where there is a difference."
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
