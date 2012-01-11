
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

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
