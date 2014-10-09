
#
# A) Define the package
#

TRIBITS_SUBPACKAGE(Classic)


SET(KokkosClassic_VERSION "1.3d")

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_OpenMP
  HAVE_KOKKOSCLASSIC_OPENMP
  "Enable OpenMP in Kokkos."
  ${${PROJECT_NAME}_ENABLE_OpenMP}
  )
IF(KokkosClassic_ENABLE_OpenMP AND NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  MESSAGE(WARNING "Disabling OpenMP in Kokkos because ${PROJECT_NAME}_ENABLE_OpenMP is ${${PROJECT_NAME}_ENABLE_OpenMP}")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_KokkosCompat
  HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT
  "Enable use of the KokkosCompat subpackage in the KokkosClassic subpackage.  This option is only for use by expert developers.  You normally do not need to set this option explicitly, even if you are using the Kokkos refactor version of Tpetra."
  "${${PROJECT_NAME}_ENABLE_KokkosCompat}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_NodeAPI
  HAVE_KOKKOSCLASSIC_NODEAPI
  "Enable Kokkos node programming API."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_Serial
  HAVE_KOKKOSCLASSIC_SERIAL
  "Add SerialNode to the list of enabled Node types"
  ON
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
# IF((KokkosClassic_ENABLE_CUDA_COMPLEX_FLOAT OR KokkosClassic_ENABLE_CUDA_COMPLEX_DOUBLE) AND NOT (Teuchos_ENABLE_COMPLEX AND KokkosClassic_ENABLE_Cusp))
#   MESSAGE(FATAL_ERROR 
#     "\nKokkos support for Complex requires enabling Teuchos support for Complex via Teuchos_ENABLE_COMPLEX and the Cusp TPL")
# ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_CUDA_NODE_MEMORY_PROFILING
  HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
  "Enable profiling of memory transfers for CUDA nodes."
  ${KokkosClassic_ENABLE_DEBUG}
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
  ${KokkosClassic_ENABLE_DEBUG}
  )

# don't need this; TriBITS does it for us courtesy of the dependency on MKL
# TRIBITS_ADD_OPTION_AND_DEFINE( 
#   KokkosClassic_ENABLE_MKL
#   HAVE_KOKKOSCLASSIC_MKL
#   "Enable Intel Math Kernel Library (MKL) support for sparse matrix kernels in Kokkos."
#   ${${PROJECT_NAME}_ENABLE_MKL}
#   )
# IF(KokkosClassic_ENABLE_MKL AND NOT ${PROJECT_NAME}_ENABLE_MKL)
#   MESSAGE(WARNING "Disabling Intel Math Kernel Library (MKL) support in Kokkos because ${PROJECT_NAME}_ENABLE_MKL is ${${PROJECT_NAME}_ENABLE_MKL}")
# ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (KokkosClassic_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(NodeAPI)
ENDIF()

if (KokkosClassic_ENABLE_LinAlg AND KokkosClassic_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(LinAlg)
ENDIF()

#
# D) Do standard postprocessing
#

TRIBITS_SUBPACKAGE_POSTPROCESS()
