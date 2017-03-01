INCLUDE(TribitsETISupport)

#
# Define the package
#
TRIBITS_PACKAGE(KokkosKernels)

# Define ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION,
# which defaults to ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Add this subpackage to the list of subpackages that do explicit
# template instantiation (ETI).
TRIBITS_ADD_ETI_SUPPORT()

# If building in debug mode, define the HAVE_KOKKOSKERNELS_DEBUG macro.
TRIBITS_ADD_DEBUG_OPTION()

#
# Set up subpackage-specific configuration options
#

#
# "Optimization level" for KokkosKernels computational kernels.  The
# higher the level, the more code variants get generated, and thus the
# longer the compile times.  However, more code variants mean both
# better performance overall, and more uniform performance for corner
# cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
#
TRIBITS_ADD_OPTION_AND_DEFINE( KokkosLinAlg_Opt_Level
  KOKKOSLINALG_OPT_LEVEL
  "Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems."
  "1"
  )

# Enable experimental features of Tpetrakernels if set at configure
# time. Default is no.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_KOKKOSKERNELS_EXPERIMENTAL
  "Enable building and installation of experimental KokkosKernels features."
  NO )


ASSERT_DEFINED (Kokkos_ENABLE_Cuda)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_CUDA
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_CUDA
  "Whether to build kernels for the execution space Kokkos::Cuda.  This option is Tpetra_INST_CUDA=${Tpetra_INST_CUDA} by default.  If expliict template instantiation (ETI) is enabled in Trilinos, disabling this when Tpetra_INST_CUDA is enabled will make Tpetra build times longer."
  ${Kokkos_ENABLE_Cuda}
  )
ASSERT_DEFINED (Kokkos_ENABLE_OpenMP)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_OPENMP
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_OPENMP
  "Whether to build kernels for the execution space Kokkos::OpenMP.  This option is Tpetra_INST_OPENMP=${Tpetra_INST_OPENMP} by default.  If expliict template instantiation (ETI) is enabled in Trilinos, disabling this when Tpetra_INST_OPENMP is enabled will make Tpetra build times longer."
  ${Kokkos_ENABLE_OpenMP}
  )
ASSERT_DEFINED (Kokkos_ENABLE_Pthread)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_PTHREAD
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_PTHREAD
  "Whether to build kernels for the execution space Kokkos::Threads.  This option is Tpetra_INST_PTHREAD=${Tpetra_INST_PTHREAD} by default.  If expliict template instantiation (ETI) is enabled in Trilinos, disabling this when Tpetra_INST_PTHREAD is enabled will make Tpetra build times longer."
  ${Kokkos_ENABLE_Pthread}
  )
ASSERT_DEFINED (Kokkos_ENABLE_Serial)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_SERIAL
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_SERIAL
  "Whether to build kernels for the execution space Kokkos::Serial.  This option is Tpetra_INST_SERIAL=${Tpetra_INST_SERIAL} by default.  If expliict template instantiation (ETI) is enabled in Trilinos, disabling this when Tpetra_INST_SERIAL is enabled will make Tpetra build times longer."
  ${Kokkos_ENABLE_Serial}
  )

# ==================================================================
# Explicit template instantiation (ETI) and test instantiation logic
# ==================================================================
SET(${PACKAGE_NAME}_ETI_SCALARS "")
SET(${PACKAGE_NAME}_ETI_LORDS "")
SET(${PACKAGE_NAME}_ETI_DEVICES "")

# ============================================================
# Kokkos execution space types
# ============================================================


#
# Construct the list of enabled Device types.
#
IF(KokkosKernels_ENABLE_EXECUTION_SPACE_SERIAL) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>")
ENDIF()
IF(KokkosKernels_ENABLE_EXECUTION_SPACE_PTHREAD) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>")
ENDIF()
IF(KokkosKernels_ENABLE_EXECUTION_SPACE_OPENMP) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>")
ENDIF()
IF(KokkosKernels_ENABLE_EXECUTION_SPACE_CUDA) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>")
ENDIF()

# ============================================================
# LocalOrdinal types
# ============================================================

# FIXME (mfh 16 Oct 2015) Tpetra currently only supports LocalOrdinal
# = int.  At some point, the list of allowed LocalOrdinal types will
# expand to include 64-bit signed integer types.

LIST (APPEND ${PACKAGE_NAME}_ETI_LORDS "int")
LIST (APPEND ${PACKAGE_NAME}_ETI_SCALARS "double")

# "Export" the names for use in the ETI system.
# If we don't do this, ETI won't see these variables.

GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS ${${PACKAGE_NAME}_ETI_SCALARS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_LORDS   ${${PACKAGE_NAME}_ETI_LORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_DEVICES ${${PACKAGE_NAME}_ETI_DEVICES})

# ==================================================================
# Process subdirectories
# ==================================================================

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PACKAGE_POSTPROCESS()

