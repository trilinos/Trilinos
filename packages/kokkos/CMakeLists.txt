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

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (Kokkos_ENABLE_NodeAPI)
  ADD_SUBDIRECTORY(NodeAPI)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
