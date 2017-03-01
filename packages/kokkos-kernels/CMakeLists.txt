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



# Enable experimental features of KokkosKernels if set at configure
# time. Default is no.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_KOKKOSKERNELS_EXPERIMENTAL
  "Enable building and installation of experimental KokkosKernels features."
  NO )

# Define what execution spaces KokkosKernels enables.
# KokkosKernels may enable fewer execution spaces than
# Kokkos enables.  This can reduce build and test times.

ASSERT_DEFINED (Kokkos_ENABLE_Cuda)
SET(${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_CUDA_DEFAULT ${Kokkos_ENABLE_Cuda})
ASSERT_DEFINED (Kokkos_ENABLE_OpenMP)
SET(${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_OPENMP_DEFAULT ${Kokkos_ENABLE_OpenMP})
ASSERT_DEFINED (Kokkos_ENABLE_Pthread)
SET(${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_PTHREAD_DEFAULT ${Kokkos_ENABLE_Pthread})
ASSERT_DEFINED (Kokkos_ENABLE_Serial)
SET(${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_SERIAL_DEFAULT ${Kokkos_ENABLE_Serial})

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_CUDA
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_CUDA
  "Whether to build kernels for the execution space Kokkos::Cuda.  This option is Kokkos_ENABLE_Cuda=${Kokkos_ENABLE_Cuda} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_Cuda is enabled may increase build times."
  ${${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_CUDA_DEFAULT}}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_OPENMP
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_OPENMP
  "Whether to build kernels for the execution space Kokkos::OpenMP.  This option is Kokkos_ENABLE_OpenMP=${Kokkos_ENABLE_OpenMP} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_OpenMP is enabled may increase build times."
  ${${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_OPENMP_DEFAULT}}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_PTHREAD
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_PTHREAD
  "Whether to build kernels for the execution space Kokkos::Threads.  This option is Kokkos_ENABLE_Pthread=${Kokkos_ENABLE_Pthread} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_Pthread is enabled may increase build times."
  ${${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_PTHREAD_DEFAULT}}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_SERIAL
  KOKKOSKERNELS_BUILD_EXECUTION_SPACE_SERIAL
  "Whether to build kernels for the execution space Kokkos::Serial.  This option is Kokkos_ENABLE_Serial=${Kokkos_ENABLE_Serial} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_Serial is enabled may increase build times."
  ${${PACKAGE_NAME}_ENABLE_EXECUTION_SPACE_SERIAL_DEFAULT}}
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

