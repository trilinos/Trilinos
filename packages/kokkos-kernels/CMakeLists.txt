INCLUDE(TribitsETISupport)

#
# Define the package
#
TRIBITS_PACKAGE(KokkosKernels)

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
SET(${PACKAGE_NAME}_INST_EXECSPACE_CUDA_DEFAULT ${Kokkos_ENABLE_Cuda})
ASSERT_DEFINED (Kokkos_ENABLE_OpenMP)
SET(${PACKAGE_NAME}_INST_EXECSPACE_OPENMP_DEFAULT ${Kokkos_ENABLE_OpenMP})
ASSERT_DEFINED (Kokkos_ENABLE_Pthread)
SET(${PACKAGE_NAME}_INST_EXECSPACE_THREADS_DEFAULT ${Kokkos_ENABLE_Pthread})
ASSERT_DEFINED (Kokkos_ENABLE_Serial)
SET(${PACKAGE_NAME}_INST_EXECSPACE_SERIAL_DEFAULT ${Kokkos_ENABLE_Serial})

IF(${Kokkos_ENABLE_Cuda})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_EXECSPACE_CUDA
  KOKKOSKERNELS_INST_EXECSPACE_CUDA
  "Whether to pre instantiate kernels for the execution space Kokkos::Cuda.  This option is KokkosKernels_INST_EXECSPACE_CUDA=${Kokkos_ENABLE_Cuda} by default.  Disabling this when Kokkos_ENABLE_Cuda is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_CUDA_DEFAULT}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_MEMSPACE_CUDAUVMSPACE
  KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE
  "Whether to pre instantiate kernels for the memory space Kokkos::CudaUVMSpace.  This option is KokkosKernels_INST_EXECSPACE_CUDAUVMSPACE=${Kokkos_ENABLE_Cuda} by default.  Disabling this when Kokkos_ENABLE_Cuda is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_CUDA_DEFAULT}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_MEMSPACE_CUDASPACE
  KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE
  "Whether to pre instantiate kernels for the memory space Kokkos::CudaSpace.  This option is KokkosKernels_INST_EXECSPACE_CUDASPACE=${Kokkos_ENABLE_Cuda} by default.  Disabling this when Kokkos_ENABLE_Cuda is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_CUDA_DEFAULT}
  )

  IF(${${PACKAGE_NAME}_INST_EXECSPACE_CUDA} AND ${${PACKAGE_NAME}_INST_MEMSPACE_CUDASPACE})
    LIST(APPEND DEVICE_LIST "<Cuda,CudaSpace>")
  ENDIF()
  IF(${${PACKAGE_NAME}_INST_EXECSPACE_CUDA} AND ${${PACKAGE_NAME}_INST_MEMSPACE_CUDAUVMSPACE})
    LIST(APPEND DEVICE_LIST "<Cuda,CudaUVMSpace>")
  ENDIF()

  IF( Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)) )
    MESSAGE( WARNING "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with CUDA and complex double enabled.")
  ENDIF()

ENDIF()

IF(${Kokkos_ENABLE_Serial} OR ${Kokkos_ENABLE_Pthread} OR ${Kokkos_ENABLE_OpenMP})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_MEMSPACE_HOSTSPACE
  KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE
  "Whether to pre instantiate kernels for the memory space Kokkos::HostSpace.  This option is KokkosKernels_INST_EXECSPACE_HOSTSPACE=(${Kokkos_ENABLE_OpenMP} OR ${Kokkos_ENABLE_Pthread} OR ${Kokkos_ENABLE_Serial}) by default.  Disabling this when one of the Host execution spaces is enabled may increase build times."
  ON
  )
ENDIF()

IF(${Kokkos_ENABLE_OpenMP})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_EXECSPACE_OPENMP
  KOKKOSKERNELS_INST_EXECSPACE_OPENMP
  "Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  This option is KokkosKernels_INST_EXECSPACE_OPENMP=${Kokkos_ENABLE_OpenMP} by default.  Disabling this when Kokkos_ENABLE_OpenMP is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_OPENMP_DEFAULT}
  )
  IF(${${PACKAGE_NAME}_INST_EXECSPACE_OPENMP} AND ${${PACKAGE_NAME}_INST_MEMSPACE_HOSTSPACE})
    LIST(APPEND DEVICE_LIST "<OpenMP,HostSpace>")
  ENDIF()
ENDIF()


IF(${Kokkos_ENABLE_Pthread})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_EXECSPACE_THREADS
  KOKKOSKERNELS_INST_EXECSPACE_THREADS
  "Whether to build kernels for the execution space Kokkos::Threads.  This option is Kokkos_ENABLE_Pthread=${Kokkos_ENABLE_Pthread} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_Pthread is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_THREADS_DEFAULT}
  )
  IF(${${PACKAGE_NAME}_INST_EXECSPACE_THREADS} AND ${${PACKAGE_NAME}_INST_MEMSPACE_HOSTSPACE})
    LIST(APPEND DEVICE_LIST "<Threads,HostSpace>")
  ENDIF()
ENDIF()


IF(${Kokkos_ENABLE_Serial})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_EXECSPACE_SERIAL
  KOKKOSKERNELS_INST_EXECSPACE_SERIAL
  "Whether to build kernels for the execution space Kokkos::Serial.  This option is Kokkos_ENABLE_Serial=${Kokkos_ENABLE_Serial} by default.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_Serial is enabled may increase build times."
  ${${PACKAGE_NAME}_INST_EXECSPACE_SERIAL_DEFAULT}
  )
  IF(${${PACKAGE_NAME}_INST_EXECSPACE_SERIAL} AND ${${PACKAGE_NAME}_INST_MEMSPACE_HOSTSPACE})
    LIST(APPEND DEVICE_LIST "<Serial,HostSpace>")
  ENDIF()
ENDIF()

# ==================================================================
# Enable Scalar Types for ETI
# ==================================================================

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ETI_Only
  KOKKOSKERNELS_ETI_ONLY
  "Whether to restrict availability of kernels to ETI types only. This is off by default, i.e. any type combination can be instantiated. Turning this on guarantees that kernels are never build inside of object files which simply call KokkosKernels functions."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Test_ETI_Only
  KOKKOSKERNELS_TEST_ETI_ONLY
  "Whether to restrict testing to ETI types."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_DOUBLE
  KOKKOSKERNELS_INST_DOUBLE
  "Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times."
  ON
  )

IF (Trilinos_ENABLE_FLOAT)
  GLOBAL_SET(KokkosKernels_INST_FLOAT_DEFAULT  ON)
ELSE()
  GLOBAL_SET(KokkosKernels_INST_FLOAT_DEFAULT  OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_FLOAT
  KOKKOSKERNELS_INST_FLOAT
  "Whether to pre instantiate kernels for the scalar type float.  This option is KokkosKernels_INST_FLOAT=Trilinos_ENABLE_FLOAT by default.  Disabling this may increase build times."
  ${KokkosKernels_INST_FLOAT_DEFAULT}
  )

IF (KokkosKernels_INST_DOUBLE AND Trilinos_ENABLE_COMPLEX_DOUBLE)
  GLOBAL_SET(KokkosKernels_INST_COMPLEX_DOUBLE_DEFAULT ON)
ELSE()
  GLOBAL_SET(KokkosKernels_INST_COMPLEX_DOUBLE_DEFAULT OFF)
ENDIF()
IF (KokkosKernels_INST_FLOAT AND Trilinos_ENABLE_COMPLEX_FLOAT)
  GLOBAL_SET(KokkosKernels_INST_COMPLEX_FLOAT_DEFAULT ON)
ELSE()
  GLOBAL_SET(KokkosKernels_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_DOUBLE
  KOKKOSKERNELS_INST_COMPLEX_DOUBLE
  "Whether to pre instantiate kernels for the scalar type complex<double>.  This option is KokkosKernels_INST_COMPLEX_DOUBLE=${Trilinos_ENABLE_COMPLEX_DOUBLE} by default.  Disabling this may increase build times."
  ${KokkosKernels_INST_COMPLEX_DOUBLE_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_FLOAT
  KOKKOSKERNELS_INST_COMPLEX_FLOAT
  "Whether to pre instantiate kernels for the scalar type complex<float>.  This option is KokkosKernels_INST_COMPLEX_FLOAT=${Trilinos_ENABLE_COMPLEX_FLOAT} by default.  Disabling this may increase build times."
  ${KokkosKernels_INST_COMPLEX_FLOAT_DEFAULT}
  )

IF (KokkosKernels_INST_DOUBLE)
  LIST(APPEND SCALAR_LIST "double")
ENDIF()

IF (KokkosKernels_INST_FLOAT)
  LIST(APPEND SCALAR_LIST "float")
ENDIF()

IF (KokkosKernels_INST_COMPLEX_DOUBLE)
  LIST(APPEND SCALAR_LIST "complex<double>")
ENDIF()

IF (KokkosKernels_INST_FLOAT)
  LIST(APPEND SCALAR_LIST "complex<float>")
ENDIF()


TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_LAYOUTLEFT
  KOKKOSKERNELS_INST_LAYOUTLEFT
  "Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is KokkosKernels_INST_LAYOUTLEFT=ON by default.  Disabling this may increase build times."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_LAYOUTRIGHT
  KOKKOSKERNELS_INST_LAYOUTRIGHT
  "Whether to pre instantiate kernels for the view layout LayoutRight.  This option is KokkosKernels_INST_LAYOUTRIGHT=OFF by default.  Disabling this may increase build times."
  OFF
  )

# ==================================================================
# Enable Ordinal Types for ETI
# ==================================================================


GLOBAL_SET(KokkosKernels_INST_ORDINAL_INT_DEFAULT  ON)
GLOBAL_SET(KokkosKernels_INST_ORDINAL_INT64_T_DEFAULT  OFF)


TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_ORDINAL_INT
  KOKKOSKERNELS_INST_ORDINAL_INT
  "Whether to pre instantiate kernels for the ordinal type int.  This option is KokkosKernels_INST_ORDINAL_INT=ON by default."
  ${KokkosKernels_INST_ORDINAL_INT_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_ORDINAL_INT64_T
  KOKKOSKERNELS_INST_ORDINAL_INT64_T
  "Whether to pre instantiate kernels for the ordinal type int64_t.  This option is KokkosKernels_INST_ORDINAL_INT64_T=OFF by default."
  ${KokkosKernels_INST_ORDINAL_INT64_T_DEFAULT}
  )

IF (KokkosKernels_INST_ORDINAL_INT)
  LIST(APPEND ORDINAL_LIST "int")
ENDIF()

IF (KokkosKernels_INST_ORDINAL_INT64_T)
  LIST(APPEND ORDINAL_LIST "int64_t")
ENDIF()


# ==================================================================
# Enable Offset Types for ETI
# ==================================================================


GLOBAL_SET(KokkosKernels_INST_OFFSET_SIZE_T_DEFAULT  ON)
GLOBAL_SET(KokkosKernels_INST_OFFSET_INT_DEFAULT  ON)
#GLOBAL_SET(KokkosKernels_INST_OFFSET_INT64_T_DEFAULT  OFF)


TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_OFFSET_INT
  KOKKOSKERNELS_INST_OFFSET_INT
  "Whether to pre instantiate kernels for the offset type int.  This option is KokkosKernels_INST_OFFSET_INT=ON by default."
  ${KokkosKernels_INST_OFFSET_INT_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_OFFSET_SIZE_T
  KOKKOSKERNELS_INST_OFFSET_SIZE_T
  "Whether to pre instantiate kernels for the offset type size_t.  This option is KokkosKernels_INST_OFFSET_SIZE_T=OFF by default."
  ${KokkosKernels_INST_OFFSET_SIZE_T_DEFAULT}
  )

IF (KokkosKernels_INST_OFFSET_INT)
  LIST(APPEND OFFSET_LIST "int")
ENDIF()

IF (KokkosKernels_INST_OFFSET_SIZE_T)
  LIST(APPEND OFFSET_LIST "size_t")
ENDIF()

# ==================================================================
# Enable Third Party Libraries
# ==================================================================

IF (TPL_ENABLE_BLAS)
  SET(KOKKOSKERNELS_ENABLE_TPL_BLAS ${TPL_ENABLE_BLAS})
ENDIF()
IF (TPL_ENABLE_MKL)
  SET(KOKKOSKERNELS_ENABLE_TPL_MKL ${TPL_ENABLE_MKL})
ENDIF()

IF(${Kokkos_ENABLE_Cuda})
  IF (NOT KOKKOSKERNELS_ENABLE_TPL_BLAS)
    SET(KOKKOSKERNELS_ENABLE_TPL_BLAS ON)
    LIST(APPEND TPL_LIST "BLAS")
  ENDIF()
  # CUBLAS is ON by default when CUDA is enabled
  SET(KOKKOSKERNELS_ENABLE_TPL_CUBLAS ON)
  # Tribit provides TPL mechanism for CUSPARSE; thus, use it
  IF (TPL_ENABLE_CUSPARSE)
    SET(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE ${TPL_ENABLE_CUSPARSE})
  ENDIF()
ENDIF()

IF (KOKKOSKERNELS_ENABLE_TPL_BLAS)
  LIST(APPEND TPL_LIST "BLAS")
ENDIF()
IF (KOKKOSKERNELS_ENABLE_TPL_MKL)
  LIST(APPEND TPL_LIST "MKL")
ENDIF()
IF (KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  LIST(APPEND TPL_LIST "CUSPARSE")
ENDIF()
IF (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
  LIST(APPEND TPL_LIST "CUBLAS")
ENDIF()

# ==================================================================
# CMake Summary
# ==================================================================


MESSAGE("")
MESSAGE("=======================")
MESSAGE("KokkosKernels ETI Types")
MESSAGE("   Devices:  ${DEVICE_LIST}")
MESSAGE("   Scalars:  ${SCALAR_LIST}")
MESSAGE("   Ordinals: ${ORDINAL_LIST}")
MESSAGE("   Offsets:  ${OFFSET_LIST}")
MESSAGE("")
MESSAGE("KokkosKernels TPLs")
MESSAGE("   ${TPL_LIST}")
MESSAGE("=======================")
MESSAGE("")

# ==================================================================
# Process subdirectories
# ==================================================================

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test_common)
TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PACKAGE_POSTPROCESS()

