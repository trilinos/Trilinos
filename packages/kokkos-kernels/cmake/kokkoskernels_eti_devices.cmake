# Define what execution spaces KokkosKernels enables.
# KokkosKernels may enable fewer execution spaces than
# Kokkos enables.  This can reduce build and test times.

SET(EXEC_SPACES
  EXECSPACE_CUDA
  EXECSPACE_HIP
  EXECSPACE_SYCL
  EXECSPACE_OPENMPTARGET
  EXECSPACE_OPENMP
  EXECSPACE_THREADS
  EXECSPACE_SERIAL
)
SET(EXECSPACE_CUDA_CPP_TYPE         Kokkos::Cuda)
SET(EXECSPACE_HIP_CPP_TYPE          Kokkos::HIP)
SET(EXECSPACE_SYCL_CPP_TYPE         Kokkos::Experimental::SYCL)
SET(EXECSPACE_OPENMPTARGET_CPP_TYPE Kokkos::Experimental::OpenMPTarget)
SET(EXECSPACE_OPENMP_CPP_TYPE       Kokkos::OpenMP)
SET(EXECSPACE_THREADS_CPP_TYPE      Kokkos::Threads)
SET(EXECSPACE_SERIAL_CPP_TYPE       Kokkos::Serial)

SET(MEM_SPACES
  MEMSPACE_CUDASPACE
  MEMSPACE_CUDAUVMSPACE
  MEMSPACE_HIPSPACE
  MEMSPACE_HIPMANAGEDSPACE
  MEMSPACE_SYCLSPACE
  MEMSPACE_SYCLSHAREDSPACE
  MEMSPACE_OPENMPTARGET
  MEMSPACE_HOSTSPACE
)
SET(MEMSPACE_CUDASPACE_CPP_TYPE         Kokkos::CudaSpace)
SET(MEMSPACE_CUDAUVMSPACE_CPP_TYPE      Kokkos::CudaUVMSpace)
SET(MEMSPACE_HIPSPACE_CPP_TYPE          Kokkos::HIPSpace)
SET(MEMSPACE_HIPMANAGEDSPACE_CPP_TYPE   Kokkos::HIPManagedSpace)
SET(MEMSPACE_SYCLSPACE_CPP_TYPE         Kokkos::Experimental::SYCLDeviceUSMSpace)
SET(MEMSPACE_SYCLSHAREDSPACE_CPP_TYPE   Kokkos::Experimental::SYCLSharedUSMSpace)
SET(MEMSPACE_OPENMPTARGETSPACE_CPP_TYPE Kokkos::Experimental::OpenMPTargetSpace)
SET(MEMSPACE_HOSTSPACE_CPP_TYPE         Kokkos::HostSpace)

IF(KOKKOS_ENABLE_CUDA)
 KOKKOSKERNELS_ADD_OPTION(
   INST_EXECSPACE_CUDA
   ON
   BOOL
   "Whether to pre instantiate kernels for the execution space Kokkos::Cuda. Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise."
   )

 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_CUDAUVMSPACE
   OFF
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::CudaUVMSpace.  Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: OFF."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_CUDASPACE
   ON
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

IF(KOKKOS_ENABLE_HIP)
 KOKKOSKERNELS_ADD_OPTION(
   INST_EXECSPACE_HIP
   ${KOKKOSKERNELS_INST_EXECSPACE_HIP_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the execution space Kokkos::HIP. Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: ON if Kokkos is HIP-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_HIPSPACE
   ${KOKKOSKERNELS_INST_EXECSPACE_HIP_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::HIPSpace.  Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: ON if Kokkos is HIP-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_HIPMANAGEDSPACE
   OFF
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::HIPManagedSpace.  Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: OFF."
   )

  IF(KOKKOSKERNELS_INST_EXECSPACE_HIP AND KOKKOSKERNELS_INST_MEMSPACE_HIPSPACE)
    LIST(APPEND DEVICE_LIST "<HIP,HIPSpace>")
  ENDIF()
  IF(KOKKOSKERNELS_INST_EXECSPACE_HIP AND KOKKOSKERNELS_INST_MEMSPACE_HIPMANAGEDSPACE)
    LIST(APPEND DEVICE_LIST "<HIP,HIPManagedSpace>")
  ENDIF()

  IF( Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)) )
    MESSAGE( WARNING "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with HIP and complex double enabled.")
  ENDIF()

ENDIF()

IF(KOKKOS_ENABLE_SYCL)
 KOKKOSKERNELS_ADD_OPTION(
   INST_EXECSPACE_SYCL
   ${KOKKOSKERNELS_INST_EXECSPACE_SYCL_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the execution space Kokkos::Experimental::SYCL. Disabling this when Kokkos_ENABLE_SYCL is enabled may increase build times. Default: ON if Kokkos is SYCL-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_SYCLSPACE
   ${KOKKOSKERNELS_INST_EXECSPACE_SYCL_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::Experimental::SYCLSpace.  Disabling this when Kokkos_ENABLE_SYCL is enabled may increase build times. Default: ON if Kokkos is SYCL-enabled, OFF otherwise."
   )

  IF(KOKKOSKERNELS_INST_EXECSPACE_SYCL AND KOKKOSKERNELS_INST_MEMSPACE_SYCLSPACE)
    LIST(APPEND DEVICE_LIST "<SYCL,SYCLDeviceUSMSpace>")
  ENDIF()
  IF(KOKKOSKERNELS_INST_EXECSPACE_SYCL AND KOKKOSKERNELS_INST_MEMSPACE_SYCLSHAREDSPACE)
    LIST(APPEND DEVICE_LIST "<SYCL,SYCLSharedUSMSpace>")
  ENDIF()

  IF( Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)) )
    MESSAGE( WARNING "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with SYCL and complex double	 enabled.")
  ENDIF()
ENDIF()

IF(KOKKOS_ENABLE_OPENMPTARGET)
 KOKKOSKERNELS_ADD_OPTION(
   INST_EXECSPACE_OPENMPTARGET
   ${KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the execution space Kokkos::Experimental::OpenMPTarget. Disabling this when Kokkos_ENABLE_OPENMPTARGET is enabled may increase build times. Default: ON if Kokkos is OpenMPTarget-enabled, OFF otherwise."
   )
 KOKKOSKERNELS_ADD_OPTION(
   INST_MEMSPACE_OPENMPTARGETSPACE
   ${KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET_DEFAULT}
   BOOL
   "Whether to pre instantiate kernels for the memory space Kokkos::Experimental::OpenMPTargetSpace.  Disabling this when Kokkos_ENABLE_OPENMPTARGET is enabled may increase build times. Default: ON if Kokkos is OpenMPTarget-enabled, OFF otherwise."
   )

  IF(KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET AND KOKKOSKERNELS_INST_MEMSPACE_OPENMPTARGETSPACE)
    LIST(APPEND DEVICE_LIST "<OpenMPTarget,OpenMPTargetSpace>")
  ENDIF()

  IF( Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)) )
    MESSAGE( WARNING "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with OpenMPTarget and complex double enabled.")
  ENDIF()
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
 INST_MEMSPACE_HOSTSPACE
 ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
 BOOL
 "Whether to pre instantiate kernels for the memory space Kokkos::HostSpace.  Disabling this when one of the Host execution spaces is enabled may increase build times. Default: ON"
)

KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_OPENMP
  ${KOKKOSKERNELS_INST_EXECSPACE_OPENMP_DEFAULT}
  BOOL
  "Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  Disabling this when Kokkos_ENABLE_OPENMP is enabled may increase build times. Default: ON if Kokkos is OpenMP-enabled, OFF otherwise."
)
IF(KOKKOSKERNELS_INST_EXECSPACE_OPENMP AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  LIST(APPEND DEVICE_LIST "<OpenMP,HostSpace>")
  IF(NOT KOKKOS_ENABLE_OPENMP)
    MESSAGE(FATAL_ERROR "Set ETI on for OPENMP, but Kokkos was not configured with the OPENMP backend")
  ENDIF()
ENDIF()


KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_THREADS
  ${KOKKOSKERNELS_INST_EXECSPACE_THREADS_DEFAULT}
  BOOL
  "Whether to build kernels for the execution space Kokkos::Threads.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_THREADS is enabled may increase build times. Default: ON if Kokkos is Threads-enabled, OFF otherwise."
)
#There continues to be name ambiguity with threads vs pthreads
SET(KOKKOSKERNELS_INST_EXECSPACE_THREADS ${KOKKOSKERNELS_INST_EXECSPACE_THREADS})

IF(KOKKOSKERNELS_INST_EXECSPACE_THREADS AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  LIST(APPEND DEVICE_LIST "<Threads,HostSpace>")
  IF(NOT KOKKOS_ENABLE_THREADS)
    MESSAGE(FATAL_ERROR "Set ETI on for THREADS, but Kokkos was not configured with the THREADS backend")
  ENDIF()
ENDIF()

KOKKOSKERNELS_ADD_OPTION(
  INST_EXECSPACE_SERIAL
  ${KOKKOSKERNELS_INST_EXECSPACE_SERIAL_DEFAULT}
  BOOL
  "Whether to build kernels for the execution space Kokkos::Serial.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_SERIAL is enabled may increase build times. Default: ON when Kokkos is Serial-enabled, OFF otherwise."
)

SET(EXECSPACE_CUDA_VALID_MEM_SPACES               CUDASPACE CUDAUVMSPACE)
SET(EXECSPACE_HIP_VALID_MEM_SPACES                HIPSPACE HIPMANAGEDSPACE)
SET(EXECSPACE_SYCL_VALID_MEM_SPACES               SYCLSPACE SYCLSHAREDSPACE)
SET(EXECSPACE_OPENMPTARGET_VALID_MEM_SPACES       OPENMPTARGETSPACE)
SET(EXECSPACE_SERIAL_VALID_MEM_SPACES             HOSTSPACE)
SET(EXECSPACE_OPENMP_VALID_MEM_SPACES             HOSTSPACE)
SET(EXECSPACE_THREADS_VALID_MEM_SPACES            HOSTSPACE)
SET(DEVICES)
FOREACH(EXEC ${EXEC_SPACES})
  IF (KOKKOSKERNELS_INST_${EXEC})
    FOREACH(MEM ${${EXEC}_VALID_MEM_SPACES})
      IF (KOKKOSKERNELS_INST_MEMSPACE_${MEM})
        LIST(APPEND DEVICES ${EXEC}_MEMSPACE_${MEM})
        SET(${EXEC}_MEMSPACE_${MEM}_CPP_TYPE "${${EXEC}_CPP_TYPE},${MEMSPACE_${MEM}_CPP_TYPE}")
        SET(KOKKOSKERNELS_INST_${EXEC}_MEMSPACE_${MEM} ON)
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
