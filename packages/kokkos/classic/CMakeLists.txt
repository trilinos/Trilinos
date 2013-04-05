INCLUDE(TribitsSubPackageMacros)

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
  KokkosClassic_ENABLE_NodeAPI
  HAVE_KOKKOSCLASSIC_NODEAPI
  "Enable Kokkos node programming API."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosClassic_ENABLE_Serial
  HAVE_KOKKOSCLASSIC_SERIAL
  "Add SerialNode to the list of enabled nodes (no by default)"
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
  "Treat Kokkos::SerialNode as a device node instead of a host node, for kernels where there is a difference."
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

######################################################################
## BEGIN TSQR options
######################################################################

# Whether to enable TSQR (Tall Skinny QR factorization).
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_TSQR
  HAVE_KOKKOSCLASSIC_TSQR
  "Enable TSQR (Tall Skinny QR factorization) in Kokkos"
  ON
  )

# Whether to enable TSQR's complex arithmetic (std::complex<T>)
# support.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Teuchos is built with complex arithmetic support.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_TSQR_Complex
  HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
  "Enable complex arithmetic (std::complex<T>) support for TSQR"
  "${Teuchos_ENABLE_COMPLEX}" 
  )

# Whether to enable TSQR's Fortran backend (TSQR::CombineFortran).
#
# Defaults to off, even if Trilinos is built with Fortran support.
# This is because TSQR's Fortran backend must support the Fortran 2003
# language standard.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_TSQR_Fortran
  HAVE_KOKKOSCLASSIC_TSQR_FORTRAN
  "For TSQR: Enable the Fortran back end option for local kernels.  This requires Fortran 90 support, so it is disabled by default, even if you have enabled Fortran in your build.  TSQR will work fine without the Fortran back end, though it may perform better on some platforms.  This option is mainly only of interest to performance tuners."
  OFF
  )

# Whether to build TbbTsqr and related classes.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Trilinos is built with the TBB (Intel's Threading Building
# Blocks) TPL (third-party library) enabled.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosClassic_ENABLE_TSQR_Intel_TBB
  HAVE_KOKKOSCLASSIC_TSQR_INTEL_TBB
  "For TSQR: enable Intel Threading Building Blocks (TBB) intranode parallelization.  This option is enabled by default if you are building Trilinos with TBB enabled as a 'third-party library' (TPL), so you should not have to enable this option manually.  TSQR will work without this, but enabling it will make TSQR go faster on multicore CPU nodes."
  "${TPL_ENABLE_TBB}" 
  )

######################################################################
## END TSQR options
######################################################################

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

######################################################################
## BEGIN TSQR setup
######################################################################

if (KokkosClassic_ENABLE_TSQR)
  # Intranode (shared-memory only) part of TSQR.
  ADD_SUBDIRECTORY(NodeTSQR)
  # Internode (distributed-memory) part of TSQR.  Depends only on
  # Teuchos for distributed-memory communication.  We put this in
  # Kokkos for lack of a better place: it's independent of Tpetra /
  # Epetra and can work for both, so it shouldn't go in either
  # package.
  ADD_SUBDIRECTORY(DistTSQR)  
ENDIF()

######################################################################
## END TSQR setup
######################################################################

#
# D) Do standard postprocessing
#

TRIBITS_SUBPACKAGE_POSTPROCESS()
