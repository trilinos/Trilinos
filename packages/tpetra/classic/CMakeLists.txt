#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Classic)

#
# Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

# The "KokkosClassic_ENABLE_OpenMP" option is for for backwards
# compatibility ONLY, so that existing scripts that set this option
# won't suddenly change behavior.  TpetraClassic_ENABLE_OpenMP
# overrides this option.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_OpenMP
  HAVE_TPETRACLASSIC_OPENMP
  "Enable OpenMP in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  This option is DEPRECATED.  Please use ${PACKAGE_NAME}_ENABLE_OpenMP instead."
  ${${PROJECT_NAME}_ENABLE_OpenMP}
  )
IF(KokkosClassic_ENABLE_OpenMP AND ${PROJECT_NAME}_ENABLE_OpenMP)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT OFF)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_OpenMP)
TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMP
  HAVE_TPETRACLASSIC_OPENMP
  "Enable the KokkosClassic::OpenMPNode Node type in ${PACKAGE_NAME} (only works if ${PROJECT_NAME}_ENABLE_OpenMP is enabled).  This Node type is broken and soon to be deprecated anyway."
  ${${PACKAGE_NAME}_ENABLE_OpenMP_DEFAULT}
  )
IF(${PACKAGE_NAME}_ENABLE_OpenMP AND NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  MESSAGE(ERROR "${PROJECT_NAME}_ENABLE_OpenMP is OFF, so it is not legal to set ${PACKAGE_NAME}_ENABLE_OpenMP to ON.")
ENDIF()

# The "KokkosClassic_ENABLE_Serial" option is for for backwards
# compatibility ONLY, so that existing scripts that set this option
# won't suddenly change behavior.  TpetraClassic_ENABLE_Serial
# overrides this option.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_Serial
  HAVE_TPETRACLASSIC_SERIAL
  "Enable KokkosClassic::SerialNode in ${PACKAGE_NAME}.  This option is DEPRECATED.  Please use ${PACKAGE_NAME}_ENABLE_Serial instead."
  ON
  )
IF(KokkosClassic_ENABLE_Serial)
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

# For backwards compatibility only.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_NodeAPI
  HAVE_KOKKOSCLASSIC_NODEAPI
  "Enable classic Kokkos node programming API.  This option is for backwards compatibility only and should never be specified explicitly."
  ON
  )

# For backwards compatibility only.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_LinAlg
  HAVE_KOKKOSCLASSIC_LINALG
  "Enable classic Kokkos local linear algebra.  This option is for backwards compatibility only and should never be specified explicitly."
  ON
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

# TRIBITS_ADD_OPTION_AND_DEFINE(
#   ${PACKAGE_NAME}_ENABLE_CUDA_COMPLEX_FLOAT
#   HAVE_TPETRACLASSIC_CUDA_COMPLEX_FLOAT
#   "Enable complex<float> type for CUDA nodes (requires Teuchos_ENABLE_COMPLEX)"
#   OFF
# )
# TRIBITS_ADD_OPTION_AND_DEFINE(
#   ${PACKAGE_NAME}_ENABLE_CUDA_COMPLEX_DOUBLE
#   HAVE_TPETRACLASSIC_CUDA_COMPLEX_DOUBLE
#   "Enable complex<double> type for CUDA nodes (requires Teuchos_ENABLE_COMPLEX)."
#   OFF
# )
# IF((${PACKAGE_NAME}_ENABLE_CUDA_COMPLEX_FLOAT OR ${PACKAGE_NAME}_ENABLE_CUDA_COMPLEX_DOUBLE) AND NOT (Teuchos_ENABLE_COMPLEX AND ${PACKAGE_NAME}_ENABLE_Cusp))
#   MESSAGE(FATAL_ERROR 
#     "\nKokkos support for Complex requires enabling Teuchos support for Complex via Teuchos_ENABLE_COMPLEX and the Cusp TPL")
# ENDIF()

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

IF (KokkosClassic_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(NodeAPI)
ENDIF()

if (KokkosClassic_ENABLE_LinAlg AND KokkosClassic_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(LinAlg)
ENDIF()

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
