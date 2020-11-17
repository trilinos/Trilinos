# Define what execution spaces KokkosKernels enables.
# KokkosKernels may enable fewer execution spaces than
# Kokkos enables.  This can reduce build and test times.

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
