

#
# A) Define the package
#

TRIBITS_PACKAGE(Ifpack  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
  HAVE_IFPACK_METIS
  "Enable support for Metis in Ifpack."
  OFF )

# Add option to enable or disable support for MPI parallel subdomain solvers
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PARALLEL_SUBDOMAIN_SOLVERS
  HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
  "Enable support for MPI parallel subdomain solvers."
  OFF )

# Add option to enable or disable the Ifpack_DynamicFactory class
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DYNAMIC_FACTORY
  HAVE_IFPACK_DYNAMIC_FACTORY
  "Enable usage of the Ifpack_DynamicFactory class."
  OFF )

#Adding a define for when the TPL Metis is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_HYPRE "${${PACKAGE_NAME}_ENABLE_HYPRE}")

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Euclid
  HAVE_EUCLID
  "Enable Ifpack/Euclid adapters"
  OFF )

#Adding a define for when the TPL HIPS is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_IFPACK_HIPS "${${PACKAGE_NAME}_ENABLE_HIPS}")

#Adding a define for when the TPL SuperLU is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_IFPACK_SUPERLU "${${PACKAGE_NAME}_ENABLE_SuperLU}")

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  IFPACK_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_FLOPCOUNTERS
  IFPACK_FLOPCOUNTERS
  "Enable Internal Flop Counters in Ifpack."
  OFF )

IF(${PACKAGE_NAME}_ENABLE_SPARSKIT)
  SET(HAVE_IFPACK_SPARSKIT ON )
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Boost)
  SET(HAVE_IFPACK_SUPPORTGRAPH ON )
ENDIF()

# "Experimental" code lives in Ifpack, but is not quite ready for use
# by a general audience.  You have to enable the CMake option below in
# order to install experimental headers, include experimental source
# files in the libraries, and build experimental tests.
#
# It's up to Ifpack developers to define what headers and source files
# are experimental, and what tests depend on experimental code.  They
# must use ${PACKAGE_NAME}_ENABLE_Experimental to conditionally
# disable tests or examples that have a required dependency on
# experimental code.  (${PACKAGE_NAME} evaluates to Ifpack, but we
# write it out so you know what to look for in the file.)  Code that
# _optionally_ depends on experimental code must be protected by
# #ifdef HAVE_IFPACK_EXPERIMENTAL ... #endif.
#
# Ifpack developers must maintain the lists of experimental header
# resp. source files in src/CMakeLists.txt, as
# ${PACKAGE_NAME}_EXPERIMENTAL_HEADERS
# resp. ${PACKAGE_NAME}_EXPERIMENTAL_SOURCES.
#
# Note that some experimental files may be excluded from the release
# tarball.  Ifpack uses the TRIBITS_EXCLUDE_FILES macro later in
# this file to list the files and directories to exclude from the
# release tarball.  That macro must NOT include conditional
# expressions; it must evaluate to the same list, no matter how
# Trilinos is configured.  We ensure this in Ifpack by listing
# explicitly the files and directories to exclude from the release
# tarball, rather than using variables or conditional CMake
# expressions.
#
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_IFPACK_EXPERIMENTAL
  "Enable building and installation of experimental Ifpack features."
  NO )

# See above discussion of "experimental" Ifpack headers and source
# files.  src/CMakeLists.txt may append to these lists.
SET(${PACKAGE_NAME}_EXPERIMENTAL_HEADERS "")
SET(${PACKAGE_NAME}_EXPERIMENTAL_SOURCES "")


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  doc/UsersGuide
  example/Ifpack_ex_ScalarLaplacian_FEM.cpp
  example/Ifpack_ex_VectorLaplacian_FEM.cpp
  example/ifpack_hb
  example/ifpack_threaded_hb
  src/Ifpack_CrsGraph.h
  src/Ifpack_CrsIlut.cpp
  src/Ifpack_CrsIlut.h
  src/Ifpack_CrsRick.cpp
  src/Ifpack_CrsRick.h
  src/Ifpack_HashTable.cpp
  src/Ifpack_OverlapFactor.*
  src/Ifpack_OverlapSolveObject..*
  src/Ifpack_PerturbedMatrix.h
  src/az_ifpack.*
  src/ifp_Block
  src/ifp_DenseMat.*
  src/ifp_GlobalPrecon.h
  src/ifp_Local.*
  src/ifp_Matrix.h
  src/ifp_Precon..*
  src/ifp_SparseUtil..*
  src/ifp_arch.h
  src/ifp_b.*
  src/ifp_c_wrappers..*
  src/ifp_ifpack.h
  src/ifp_lapackd.h
  src/ifp_sp.*
  src/old.Makefile
  src/stamp-h.in
  src/xxemacs
  test/PointPreconditioner
  test/scripts
  test/scripts/run-tests.sh
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
