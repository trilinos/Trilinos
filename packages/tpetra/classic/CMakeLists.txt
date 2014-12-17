#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Classic)

#
# Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_OpenMP
  HAVE_KOKKOSCLASSIC_OPENMP
  "Enable OpenMP in ${PACKAGE_NAME}."
  ${${PROJECT_NAME}_ENABLE_OpenMP}
  )
IF(${PACKAGE_NAME}_ENABLE_OpenMP AND NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  MESSAGE(WARNING "Disabling OpenMP in Kokkos because ${PROJECT_NAME}_ENABLE_OpenMP is ${${PROJECT_NAME}_ENABLE_OpenMP}")
ENDIF()

# For backwards compatibility ONLY.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_OpenMP
  HAVE_KOKKOSCLASSIC_OPENMP
  "Enable OpenMP in ${PACKAGE_NAME}."
  ${${PROJECT_NAME}_ENABLE_OpenMP}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_KokkosCompat
  HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT
  "Enable use of the KokkosCompat subpackage in the ${PACKAGE_NAME} subpackage.  This option is only for use by expert developers.  You normally do not need to set this option explicitly, even if you are using the Kokkos refactor version of Tpetra."
  "${${PROJECT_NAME}_ENABLE_KokkosCompat}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_NodeAPI
  HAVE_KOKKOSCLASSIC_NODEAPI
  "Enable Kokkos node programming API."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE( 
  ${PACKAGE_NAME}_ENABLE_Serial
  HAVE_KOKKOSCLASSIC_SERIAL
  "Add SerialNode to the list of enabled Node types"
  ON
  )

# For backwards compatibility ONLY.
TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_Serial
  HAVE_KOKKOSCLASSIC_SERIAL
  "Add SerialNode to the list of enabled Node types"
  ${${PACKAGE_NAME}_ENABLE_Serial}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_LinAlg
  HAVE_KOKKOSCLASSIC_LINALG
  "Enable Kokkos local linear algebra."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_CUDA_FLOAT
  HAVE_KOKKOSCLASSIC_CUDA_FLOAT
  "Enable float type for CUDA nodes."
  ON
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_CUDA_DOUBLE
  HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
  "Enable double type for CUDA nodes."
  OFF
)

# TRIBITS_ADD_OPTION_AND_DEFINE(
#   KokkosClassic_ENABLE_CUDA_COMPLEX_FLOAT
#   HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT
#   "Enable complex<float> type for CUDA nodes (requires Teuchos_ENABLE_COMPLEX)"
#   OFF
# )
# TRIBITS_ADD_OPTION_AND_DEFINE(
#   KokkosClassic_ENABLE_CUDA_COMPLEX_DOUBLE
#   HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE
#   "Enable complex<double> type for CUDA nodes (requires Teuchos_ENABLE_COMPLEX)."
#   OFF
# )
# IF((KokkosClassic_ENABLE_CUDA_COMPLEX_FLOAT OR KokkosClassic_ENABLE_CUDA_COMPLEX_DOUBLE) AND NOT (Teuchos_ENABLE_COMPLEX AND ${PACKAGE_NAME}_ENABLE_Cusp))
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
