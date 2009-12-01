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
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
