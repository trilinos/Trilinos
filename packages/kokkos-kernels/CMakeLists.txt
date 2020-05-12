

IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOSKERNELS_HAS_TRILINOS ON CACHE BOOL "")
ELSE()
  SET(KOKKOSKERNELS_HAS_TRILINOS OFF CACHE BOOL "")
ENDIF()

GET_DIRECTORY_PROPERTY(KOKKOSKERNELS_HAS_PARENT PARENT_DIRECTORY)

SET(KOKKOSKERNELS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(KOKKOSKERNELS_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})


IF(NOT KOKKOSKERNELS_HAS_TRILINOS)
  cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
  IF (Spack_WORKAROUND) 
    #if we are explicitly using Spack for development,
    #nuke the Spack compiler
    SET(SPACK_CXX $ENV{SPACK_CXX})
    IF(SPACK_CXX)
      SET(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
      SET(ENV{CXX} ${SPACK_CXX})
    ENDIF()
  ENDIF()
  IF(NOT DEFINED ${PROJECT_NAME})
    PROJECT(KokkosKernels CXX)
  ENDIF()
  SET(KokkosKernels_VERSION_MAJOR 3)
  SET(KokkosKernels_VERSION_MINOR 1)
  SET(KokkosKernels_VERSION_PATCH 1)
ENDIF()

IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0") 
  MESSAGE(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()

INCLUDE(GNUInstallDirs)
IF (KOKKOSKERNELS_HAS_TRILINOS)
 SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR ${TRILINOS_INCDIR})
 SET(KOKKOS_ENABLE_CUDA_UVM ${Kokkos_ENABLE_CUDA_UVM})
ELSEIF(KOKKOSKERNELS_HAS_PARENT)
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/kokkos-kernels")
 SET(KOKKOS_ENABLE_CUDA_UVM ${Kokkos_ENABLE_CUDA_UVM})
ELSE()
 FIND_PACKAGE(Kokkos REQUIRED)
 MESSAGE(STATUS "Found Kokkos at ${Kokkos_DIR}")
 SET(KOKKOSKERNELS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
 KOKKOS_CHECK(OPTIONS CUDA_UVM RETURN_VALUE KOKKOS_ENABLE_CUDA_UVM)
ENDIF()

INCLUDE(cmake/fake_tribits.cmake)
INCLUDE(cmake/kokkoskernels_tribits.cmake)

#This will collect imported TPLs we need to export
#in our config file
GLOBAL_SET(KOKKOSKERNELS_TPL_EXPORTS)

#
# Define the package
#
KOKKOSKERNELS_PACKAGE()

# If building in debug mode, define the HAVE_KOKKOSKERNELS_DEBUG macro.
KOKKOSKERNELS_ADD_DEBUG_OPTION()

#
# "Optimization level" for KokkosKernels computational kernels.  The
# higher the level, the more code variants get generated, and thus the
# longer the compile times.  However, more code variants mean both
# better performance overall, and more uniform performance for corner
# cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
#
KOKKOSKERNELS_ADD_OPTION_AND_DEFINE(
  LINALG_OPT_LEVEL
  KOKKOSLINALG_OPT_LEVEL
  "Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems. Default: 1"
  "1")



# Enable experimental features of KokkosKernels if set at configure
# time. Default is no.
KOKKOSKERNELS_ADD_OPTION_AND_DEFINE(
  ENABLE_EXPERIMENTAL
  HAVE_KOKKOSKERNELS_EXPERIMENTAL
  "Enable building and installation of experimental KokkosKernels features. Default: OFF"
  OFF)

# Define what execution spaces KokkosKernels enables.
# KokkosKernels may enable fewer execution spaces than
# Kokkos enables.  This can reduce build and test times.

# Kokkos only defines the variables if the backends are ON
# Define aux variables that exist as on/off
MACRO(CHECK_KOKKOS_BACKEND BE)
  IF (Kokkos_ENABLE_${BE})
    SET(KOKKOS_ENABLE_${BE} ON)
  ELSE()
    SET(KOKKOS_ENABLE_${BE} OFF)
  ENDIF()
  SET(KOKKOSKERNELS_INST_EXECSPACE_${BE}_DEFAULT ${KOKKOS_ENABLE_${BE}})
ENDMACRO(CHECK_KOKKOS_BACKEND)


SET(EXEC_SPACES
  EXECSPACE_CUDA
  EXECSPACE_OPENMP
  EXECSPACE_PTHREAD
  EXECSPACE_SERIAL
)
SET(EXECSPACE_CUDA_CPP_TYPE Kokkos::Cuda)
SET(EXECSPACE_OPENMP_CPP_TYPE Kokkos::OpenMP)
SET(EXECSPACE_PTHREAD_CPP_TYPE Kokkos::Threads)
SET(EXECSPACE_SERIAL_CPP_TYPE Kokkos::Serial)

SET(MEM_SPACES
  MEMSPACE_CUDASPACE
  MEMSPACE_CUDAUVMSPACE
  MEMSPACE_HOSTSPACE
  MEMSPACE_HBWSPACE
)
SET(MEMSPACE_CUDASPACE_CPP_TYPE    Kokkos::CudaSpace)
SET(MEMSPACE_CUDAUVMSPACE_CPP_TYPE Kokkos::CudaUVMSpace)
SET(MEMSPACE_HOSTSPACE_CPP_TYPE    Kokkos::HostSpace)
SET(MEMSPACE_HBWSPACE_CPP_TYPE     Kokkos::HBWSpace)

CHECK_KOKKOS_BACKEND(CUDA)
CHECK_KOKKOS_BACKEND(OPENMP)
CHECK_KOKKOS_BACKEND(PTHREAD)
CHECK_KOKKOS_BACKEND(SERIAL)

KOKKOSKERNELS_ADD_OPTION(
  ADD_DEFAULT_ETI
  ON
  BOOL
  "Whether to include a set of default ETI instantiations (otherwise only those explicitly requested will be included"
)

IF(KOKKOS_ENABLE_CUDA)
 KOKKOSKERNELS_ADD_OPTION(
   INST_EXECSPACE_CUDA
   ${KOKKOSKERNELS_INST_EXECSPACE_CUDA_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the execution space Kokkos::Cuda. Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_CUDAUVMSPACE
   ${KOKKOSKERNELS_INST_EXECSPACE_CUDA_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::CudaUVMSpace.  Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_CUDASPACE
   ${KOKKOSKERNELS_INST_EXECSPACE_CUDA_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::CudaSpace.  Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise."
   )

  IF(KOKKOSKERNELS_INST_EXECSPACE_CUDA AND KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE)
    LIST(APPEND DEVICE_LIST "<Cuda,CudaSpace>")
  ENDIF()
  IF(KOKKOSKERNELS_INST_EXECSPACE_CUDA AND KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
    LIST(APPEND DEVICE_LIST "<Cuda,CudaUVMSpace>")
  ENDIF()

  IF( Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)) )
    MESSAGE( WARNING "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with CUDA and complex double enabled.")
  ENDIF()

ENDIF()

KOKKOSKERNELS_ADD_OPTION(
 INST_MEMSPACE_HOSTSPACE
 ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
 BOOL
 "Whether to pre instantiate kernels for the memory space Kokkos::HostSpace.  Disabling this when one of the Host execution spaces is enabled may increase build times. Default: ON"
)

KOKKOSKERNELS_ADD_OPTION(
 INST_MEMSPACE_HBWSPACE
 OFF
 BOOL
 "Whether to pre instantiate kernels for the memory space Kokkos::HBWSpace."
)

KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_OPENMP
  ${KOKKOSKERNELS_INST_EXECSPACE_OPENMP_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  Disabling this when Kokkos_ENABLE_OpenMP is enabled may increase build times. Default: ON if Kokkos is OpenMP-enabled, OFF otherwise."
)
IF(KOKKOSKERNELS_INST_EXECSPACE_OPENMP AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  LIST(APPEND DEVICE_LIST "<OpenMP,HostSpace>")
  IF(NOT KOKKOS_ENABLE_OPENMP)
    MESSAGE(FATAL_ERROR "Set ETI on for OPENMP, but Kokkos was not configured with the OPENMP backend")
  ENDIF()
ENDIF()


KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_THREADS
  ${KOKKOSKERNELS_INST_EXECSPACE_PTHREAD_DEFAULT}
  BOOL
  "Whether to build kernels for the execution space Kokkos::Threads.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_PTHREAD is enabled may increase build times. Default: ON if Kokkos is Threads-enabled, OFF otherwise."
)
#There continues to be name ambiguity with threads vs pthreads
SET(KOKKOSKERNELS_INST_EXECSPACE_PTHREAD ${KOKKOSKERNELS_INST_EXECSPACE_THREADS})

IF(KOKKOSKERNELS_INST_EXECSPACE_PTHREAD AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  LIST(APPEND DEVICE_LIST "<Threads,HostSpace>")
  IF(NOT KOKKOS_ENABLE_PTHREAD)
    MESSAGE(FATAL_ERROR "Set ETI on for PTHREAD, but Kokkos was not configured with the PTHREAD backend")
  ENDIF()
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_SERIAL
  ${KOKKOSKERNELS_INST_EXECSPACE_SERIAL_DEFAULT}
  BOOL
  "Whether to build kernels for the execution space Kokkos::Serial.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_SERIAL is enabled may increase build times. Default: ON when Kokkos is Serial-enabled, OFF otherwise."
)

SET(EXECSPACE_CUDA_VALID_MEM_SPACES      CUDASPACE CUDAUVMSPACE)
SET(EXECSPACE_SERIAL_VALID_MEM_SPACES    HBWSPACE HOSTSPACE)
SET(EXECSPACE_OPENMP_VALID_MEM_SPACES    HBWSPACE HOSTSPACE)
SET(EXECSPACE_PTHREAD_VALID_MEM_SPACES   HBWSPACE HOSTSPACE)
SET(DEVICES)
FOREACH(EXEC ${EXEC_SPACES})
  IF (KOKKOSKERNELS_INST_${EXEC})
    FOREACH(MEM ${${EXEC}_VALID_MEM_SPACES})
      IF (KOKKOSKERNELS_INST_MEMSPACE_${MEM})
        LIST(APPEND DEVICES ${EXEC}_MEMSPACE_${MEM})
        SET(${EXEC}_MEMSPACE_${MEM}_CPP_TYPE "${${EXEC}_CPP_TYPE},${MEMSPACE_${MEM}_CPP_TYPE}")
        SET(KOKKOSKERNELS_INST_${EXEC}_MEMSPACE_${MEM} ON)
        FOREACH(MEM2 ${${EXEC}_VALID_MEM_SPACES})
          IF (KOKKOSKERNELS_INST_MEMSPACE_${MEM2})
            LIST(APPEND DEVICES_W_SLOW_SPACE ${EXEC}_MEMSPACE_${MEM}_MEMSPACE_${MEM2})
            SET(${EXEC}_MEMSPACE_${MEM}_MEMSPACE_${MEM2}_CPP_TYPE
              "${${EXEC}_CPP_TYPE},${MEMSPACE_${MEM}_CPP_TYPE},${MEMSPACE_${MEM2}_CPP_TYPE}")
            SET(KOKKOSKERNELS_INST_${EXEC}_MEMSPACE_${MEM}_MEMSPACE_${MEM2} ON)
          ENDIF()
        ENDFOREACH()
      ENDIF()
    ENDFOREACH()
  ENDIF()
ENDFOREACH()


IF(KOKKOSKERNELS_INST_EXECSPACE_SERIAL AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  LIST(APPEND DEVICE_LIST "<Serial,HostSpace>")
  IF(NOT KOKKOS_ENABLE_SERIAL)
    MESSAGE(FATAL_ERROR "Set ETI on for SERIAL, but Kokkos was not configured with the SERIAL backend")
  ENDIF()
ENDIF()

# ==================================================================
# Enable Scalar Types for ETI
# ==================================================================

KOKKOSKERNELS_ADD_OPTION(
  ETI_ONLY
  OFF
  BOOL
  "Whether to restrict availability of kernels to ETI types only. Turning this on guarantees that kernels are never built inside of object files which simply call KokkosKernels functions. Default: OFF"
  )

KOKKOSKERNELS_ADD_OPTION(
  TEST_ETI_ONLY
  ON
  BOOL
  "Whether to restrict testing to ETI types. Default: ON"
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_DOUBLE
  ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
  BOOL
  "Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times. Default: ON"
  )

IF (Trilinos_ENABLE_FLOAT)
  SET(KOKKOSKERNELS_INST_FLOAT_DEFAULT  ON)
ELSE()
  SET(KOKKOSKERNELS_INST_FLOAT_DEFAULT  OFF)
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_FLOAT
  ${KOKKOSKERNELS_INST_FLOAT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type float.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_FLOAT."
  )

SET(FLOATS
  FLOAT
  DOUBLE
  COMPLEX_FLOAT
  COMPLEX_DOUBLE)
SET(DOUBLE_CPP_TYPE "double")
SET(FLOAT_CPP_TYPE "float")
SET(COMPLEX_FLOAT_CPP_TYPE "Kokkos::complex<float>")
SET(COMPLEX_DOUBLE_CPP_TYPE "Kokkos::complex<double>")

#Just leave the Trilinos logic here alone
#Outside of Trilinos, let the defaults be always OFF
IF (KOKKOSKERNELS_INST_DOUBLE AND Trilinos_ENABLE_COMPLEX_DOUBLE)
  SET(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT ON)
ELSE()
  SET(KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT OFF)
ENDIF()
IF (KOKKOSKERNELS_INST_FLOAT AND Trilinos_ENABLE_COMPLEX_FLOAT)
  SET(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT ON)
ELSE()
  SET(KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_COMPLEX_DOUBLE
  ${KOKKOSKERNELS_INST_COMPLEX_DOUBLE_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type complex<double>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_DOUBLE."
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_COMPLEX_FLOAT
  ${KOKKOSKERNELS_INST_COMPLEX_FLOAT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the scalar type complex<float>.  Disabling this may increase build times. Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_FLOAT."
  )

IF (KOKKOSKERNELS_INST_DOUBLE)
  LIST(APPEND SCALAR_LIST "double")
ENDIF()

IF (KOKKOSKERNELS_INST_FLOAT)
  LIST(APPEND SCALAR_LIST "float")
ENDIF()

IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
  LIST(APPEND SCALAR_LIST "complex<double>")
ENDIF()

IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
  LIST(APPEND SCALAR_LIST "complex<float>")
ENDIF()

SET(LAYOUTS
  LAYOUTLEFT
  LAYOUTRIGHT)
SET(LAYOUTLEFT_CPP_TYPE Kokkos::LayoutLeft)
SET(LAYOUTRIGHT_CPP_TYPE Kokkos::LayoutRight)

KOKKOSKERNELS_ADD_OPTION(
  INST_LAYOUTLEFT
  ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
  BOOL
  "Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is KokkosKernels_INST_LAYOUTLEFT=ON by default.  Disabling this may increase build times. Default: ON"
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_LAYOUTRIGHT
  OFF
  BOOL
  "Whether to pre instantiate kernels for the view layout LayoutRight.  This option is KokkosKernels_INST_LAYOUTRIGHT=OFF by default.  Disabling this may increase build times. Default: OFF"
  )

IF (KOKKOSKERNELS_INST_LAYOUTLEFT)
  LIST(APPEND LAYOUT_LIST "LayoutLeft")
ENDIF()

IF (KOKKOSKERNELS_INST_LAYOUTRIGHT)
  LIST(APPEND LAYOUT_LIST "LayoutRight")
ENDIF()

# ==================================================================
# Enable Ordinal Types for ETI
# ==================================================================

SET(ORDINALS
  ORDINAL_INT
  ORDINAL_INT64_T
)
SET(ORDINAL_INT_CPP_TYPE int)
SET(ORDINAL_INT64_T_CPP_TYPE int64_t)

SET(KOKKOSKERNELS_INST_ORDINAL_INT_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
SET(KOKKOSKERNELS_INST_ORDINAL_INT64_T_DEFAULT  OFF)


KOKKOSKERNELS_ADD_OPTION(
  INST_ORDINAL_INT
  ${KOKKOSKERNELS_INST_ORDINAL_INT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the ordinal type int.  This option is KokkosKernels_INST_ORDINAL_INT=ON by default. Default: ON"
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_ORDINAL_INT64_T
  ${KOKKOSKERNELS_INST_ORDINAL_INT64_T_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the ordinal type int64_t.  This option is KokkosKernels_INST_ORDINAL_INT64_T=OFF by default. Default: OFF"
  )

IF (KOKKOSKERNELS_INST_ORDINAL_INT)
  LIST(APPEND ORDINAL_LIST "int")
ENDIF()

IF (KOKKOSKERNELS_INST_ORDINAL_INT64_T)
  LIST(APPEND ORDINAL_LIST "int64_t")
ENDIF()


# ==================================================================
# Enable Offset Types for ETI
# ==================================================================


SET(KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
SET(KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT ${KOKKOSKERNELS_ADD_DEFAULT_ETI})
SET(OFFSETS
  OFFSET_INT
  OFFSET_SIZE_T
)
SET(OFFSET_INT_CPP_TYPE int)
SET(OFFSET_SIZE_T_CPP_TYPE size_t)
#GLOBAL_SET(KokkosKernels_INST_OFFSET_INT64_T_DEFAULT  OFF)

KOKKOSKERNELS_ADD_OPTION(
  INST_OFFSET_INT
  ${KOKKOSKERNELS_INST_OFFSET_INT_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the offset type int.  This option is KokkosKernels_INST_OFFSET_INT=ON by default. Default: ON"
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_OFFSET_SIZE_T
  ${KOKKOSKERNELS_INST_OFFSET_SIZE_T_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the offset type size_t.  This option is KokkosKernels_INST_OFFSET_SIZE_T=OFF by default. Default: ON"
  )

IF (KOKKOSKERNELS_INST_OFFSET_INT)
  LIST(APPEND OFFSET_LIST "int")
ENDIF()

IF (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  LIST(APPEND OFFSET_LIST "size_t")
ENDIF()

# ==================================================================
# Enable Third Party Libraries
# ==================================================================

IF (NOT KOKKOSKERNELS_HAS_TRILINOS)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
  KOKKOSKERNELS_ADD_OPTION(
    "ENABLE_EXAMPLES"
    OFF
    BOOL
    "Whether to build examples. Default: OFF"
  )
  KOKKOSKERNELS_ADD_OPTION(
    "ENABLE_TESTS"
    OFF
    BOOL
    "Whether to build tests. Default: OFF"
  )
ENDIF()
INCLUDE(cmake/kokkoskernels_tpls.cmake)
INCLUDE(cmake/kokkoskernels_features.cmake)

# ==================================================================
# Fortran Complex BLAS
# ==================================================================

IF (KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL)
  INCLUDE(CheckHostBlasReturnComplex.cmake)
  CHECK_HOST_BLAS_RETURN_COMPLEX(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
ENDIF()

# ==================================================================
# CMake Summary
# ==================================================================

INCLUDE(cmake/kokkos_requirements.cmake)

MESSAGE("")
MESSAGE("=======================")
MESSAGE("KokkosKernels ETI Types")
MESSAGE("   Devices:  ${DEVICE_LIST}")
MESSAGE("   Scalars:  ${SCALAR_LIST}")
MESSAGE("   Ordinals: ${ORDINAL_LIST}")
MESSAGE("   Offsets:  ${OFFSET_LIST}")
MESSAGE("   Layouts:  ${LAYOUT_LIST}")
MESSAGE("")
MESSAGE("KokkosKernels TPLs")
FOREACH(TPL ${KOKKOSKERNELS_TPL_LIST})
  PAD_STRING("${TPL}:" TPL_PADDED 12)
  MESSAGE("   ${TPL_PADDED} ${${TPL}_FOUND_INFO}")
ENDFOREACH()
MESSAGE("=======================")
MESSAGE("")

# ==================================================================
# Process subdirectories
# ==================================================================

ADD_SUBDIRECTORY(src)

KOKKOSKERNELS_ADD_TEST_DIRECTORIES(test_common)
KOKKOSKERNELS_ADD_TEST_DIRECTORIES(perf_test)
KOKKOSKERNELS_ADD_TEST_DIRECTORIES(unit_test)
#KOKKOSKERNELS_ADD_EXAMPLE_DIRECTORIES(example)

KOKKOSKERNELS_PACKAGE_POSTPROCESS()


