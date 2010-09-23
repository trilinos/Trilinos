INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

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

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
