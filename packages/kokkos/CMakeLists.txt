INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(Kokkos)

SET(${PACKAGE_NAME}_VERSION "1.3d")

#
# B) Set up package-specific options
#

PACKAGE_ADD_DEBUG_OPTION()

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_NodeAPI
  HAVE_KOKKOS_NODEAPI
  "Enable Kokkos node programming API."
  ON
  )

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_LinAlg
  HAVE_KOKKOS_LINALG
  "Enable Kokkos local linear algebra."
  ON
  )
ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Array
  HAVE_KOKKOS_ARRAY
  "Enable Kokkos multidimensional arrays."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_CUDA_FLOAT
  HAVE_KOKKOS_CUDA_FLOAT
  "Enable float type for CUDA nodes."
  ON
)

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_CUDA_DOUBLE
  HAVE_KOKKOS_CUDA_DOUBLE
  "Enable double type for CUDA nodes."
  OFF
)

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_CUDA_COMPLEX_FLOAT
  HAVE_KOKKOS_CUDA_COMPLEX_FLOAT
  "Enable complex<float> type for CUDA nodes."
  OFF
)

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_CUDA_COMPLEX_DOUBLE
  HAVE_KOKKOS_CUDA_COMPLEX_DOUBLE
  "Enable complex<double> type for CUDA nodes."
  OFF
)

ADD_OPTION_AND_DEFINE(
  Kokkos_CUDA_NODE_MEMORY_PROFILING
  HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
  "Enable profiling of memory transfers for CUDA nodes."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
)

ADD_OPTION_AND_DEFINE(
  Kokkos_CUDA_NODE_MEMORY_TRACE
  HAVE_KOKKOS_CUDA_NODE_MEMORY_TRACE
  "Enable memory movement trace for CUDA nodes."
  OFF
)

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_TREAT_SERIALNODE_AS_DEVICE
  HAVE_KOKKOS_TREAT_SERIALNODE_AS_DEVICE
  "Treat Kokkos::SerialNode as a device node instead of a host node, for kernels where there is a difference."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
  )

######################################################################
## BEGIN TSQR options
######################################################################

# Whether to enable TSQR (Tall Skinny QR factorization).
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_KOKKOS_TSQR
  "Enable TSQR (Tall Skinny QR factorization) in Kokkos"
  OFF
  )

# Whether to enable TSQR's complex arithmetic (std::complex<T>)
# support.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Teuchos is built with complex arithmetic support.
ADD_OPTION_AND_DEFINE(
  TSQR_ENABLE_Complex
  HAVE_TSQR_COMPLEX
  "Enable complex arithmetic (std::complex<T>) support for TSQR"
  "${Teuchos_ENABLE_COMPLEX}" 
  )

# Whether to enable TSQR's Fortran backend (TSQR::CombineFortran).
#
# Defaults to off, even if Trilinos is built with Fortran support.
# This is because TSQR's Fortran backend must support the Fortran 2003
# language standard.
ADD_OPTION_AND_DEFINE(
  TSQR_ENABLE_Fortran
  HAVE_TSQR_FORTRAN
  "For TSQR: Enable the CombineFortran backend option for Combine.  This requires Fortran 2003 support, so it is disabled by default, even if you have enabled Fortran in your build.  TSQR will work fine without CombineFortran; this option is only of interest to performance tuners."
  OFF
  )

# Whether to build TbbTsqr and related classes.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Trilinos is built with the TBB (Intel's Threading Building
# Blocks) TPL (third-party library) enabled.
ADD_OPTION_AND_DEFINE(
  TSQR_ENABLE_Intel_TBB
  HAVE_TSQR_INTEL_TBB
  "For TSQR: enable Intel Threading Building Blocks intranode parallelization.  This option is enabled by default if you are building Trilinos with TBB enabled as a 'third-party library' (TPL), so you should not have to enable this option manually.  TSQR will work without this, but enabling it will make TSQR go faster on multicore CPU nodes."
  "${TPL_ENABLE_TBB}" 
  )

######################################################################
## END TSQR options
######################################################################

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (Kokkos_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(NodeAPI)
ENDIF()

if (Kokkos_ENABLE_LinAlg AND Kokkos_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(LinAlg)
ENDIF()

IF (Kokkos_ENABLE_Array)
  ADD_SUBDIRECTORY(array)
ENDIF()

######################################################################
## BEGIN TSQR setup
######################################################################

if (Kokkos_ENABLE_TSQR)
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
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc
  DistTSQR/doc
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
