# Define what execution spaces KokkosKernels enables.
# KokkosKernels may enable fewer execution spaces than
# Kokkos enables.  This can reduce build and test times.

set(EXEC_SPACES
    EXECSPACE_CUDA
    EXECSPACE_HIP
    EXECSPACE_SYCL
    EXECSPACE_OPENMPTARGET
    EXECSPACE_OPENMP
    EXECSPACE_THREADS
    EXECSPACE_SERIAL)
set(EXECSPACE_CUDA_CPP_TYPE          Kokkos::Cuda)
set(EXECSPACE_HIP_CPP_TYPE           Kokkos::HIP)
set(EXECSPACE_SYCL_CPP_TYPE          Kokkos::Experimental::SYCL)
set(EXECSPACE_OPENMPTARGET_CPP_TYPE  Kokkos::Experimental::OpenMPTarget)
set(EXECSPACE_OPENMP_CPP_TYPE        Kokkos::OpenMP)
set(EXECSPACE_THREADS_CPP_TYPE       Kokkos::Threads)
set(EXECSPACE_SERIAL_CPP_TYPE        Kokkos::Serial)

set(MEM_SPACES
    MEMSPACE_CUDASPACE
    MEMSPACE_CUDAUVMSPACE
    MEMSPACE_HIPSPACE
    MEMSPACE_HIPMANAGEDSPACE
    MEMSPACE_SYCLSPACE
    MEMSPACE_SYCLSHAREDSPACE
    MEMSPACE_OPENMPTARGET
    MEMSPACE_HOSTSPACE)
set(MEMSPACE_CUDASPACE_CPP_TYPE          Kokkos::CudaSpace)
set(MEMSPACE_CUDAUVMSPACE_CPP_TYPE       Kokkos::CudaUVMSpace)
set(MEMSPACE_HIPSPACE_CPP_TYPE           Kokkos::HIPSpace)
set(MEMSPACE_HIPMANAGEDSPACE_CPP_TYPE    Kokkos::HIPManagedSpace)
set(MEMSPACE_SYCLSPACE_CPP_TYPE          Kokkos::Experimental::SYCLDeviceUSMSpace)
set(MEMSPACE_SYCLSHAREDSPACE_CPP_TYPE    Kokkos::Experimental::SYCLSharedUSMSpace)
set(MEMSPACE_OPENMPTARGETSPACE_CPP_TYPE  Kokkos::Experimental::OpenMPTargetSpace)
set(MEMSPACE_HOSTSPACE_CPP_TYPE          Kokkos::HostSpace)

if(KOKKOS_ENABLE_CUDA)
  kokkoskernels_add_option("INST_EXECSPACE_CUDA" ON BOOL
    "Whether to pre instantiate kernels for the execution space Kokkos::Cuda. Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise.")

  kokkoskernels_add_option("INST_MEMSPACE_CUDAUVMSPACE" OFF BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::CudaUVMSpace.  Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: OFF.")

  kokkoskernels_add_option("INST_MEMSPACE_CUDASPACE" ON BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::CudaSpace.  Disabling this when Kokkos_ENABLE_CUDA is enabled may increase build times. Default: ON if Kokkos is CUDA-enabled, OFF otherwise.")

  if(KOKKOSKERNELS_INST_EXECSPACE_CUDA AND KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE)
    list(APPEND DEVICE_LIST "<Cuda,CudaSpace>")
  endif()
  if(KOKKOSKERNELS_INST_EXECSPACE_CUDA AND KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
    list(APPEND DEVICE_LIST "<Cuda,CudaUVMSpace>")
  endif()

  if(Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)))
    message(WARNING
      "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with CUDA and complex double enabled."
    )
  endif()

endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_option("INST_EXECSPACE_HIP" ${KOKKOSKERNELS_INST_EXECSPACE_HIP_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the execution space Kokkos::HIP. Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: ON if Kokkos is HIP-enabled, OFF otherwise.")

  kokkoskernels_add_option("INST_MEMSPACE_HIPSPACE" ${KOKKOSKERNELS_INST_EXECSPACE_HIP_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::HIPSpace.  Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: ON if Kokkos is HIP-enabled, OFF otherwise.")
  kokkoskernels_add_option("INST_MEMSPACE_HIPMANAGEDSPACE" OFF BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::HIPManagedSpace.  Disabling this when Kokkos_ENABLE_HIP is enabled may increase build times. Default: OFF.")

  if(KOKKOSKERNELS_INST_EXECSPACE_HIP AND KOKKOSKERNELS_INST_MEMSPACE_HIPSPACE)
    list(APPEND DEVICE_LIST "<HIP,HIPSpace>")
  endif()
  if(KOKKOSKERNELS_INST_EXECSPACE_HIP AND KOKKOSKERNELS_INST_MEMSPACE_HIPMANAGEDSPACE)
    list(APPEND DEVICE_LIST "<HIP,HIPManagedSpace>")
  endif()

  if(Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)))
    message(WARNING
      "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with HIP and complex double enabled.")
  endif()

endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_option("INST_EXECSPACE_SYCL" ${KOKKOSKERNELS_INST_EXECSPACE_SYCL_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the execution space Kokkos::Experimental::SYCL. Disabling this when Kokkos_ENABLE_SYCL is enabled may increase build times. Default: ON if Kokkos is SYCL-enabled, OFF otherwise.")

  kokkoskernels_add_option("INST_MEMSPACE_SYCLSPACE" ${KOKKOSKERNELS_INST_EXECSPACE_SYCL_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::Experimental::SYCLSpace.  Disabling this when Kokkos_ENABLE_SYCL is enabled may increase build times. Default: ON if Kokkos is SYCL-enabled, OFF otherwise.")

  if(KOKKOSKERNELS_INST_EXECSPACE_SYCL AND KOKKOSKERNELS_INST_MEMSPACE_SYCLSPACE)
    list(APPEND DEVICE_LIST "<SYCL,SYCLDeviceUSMSpace>")
  endif()
  if(KOKKOSKERNELS_INST_EXECSPACE_SYCL AND KOKKOSKERNELS_INST_MEMSPACE_SYCLSHAREDSPACE)
    list(APPEND DEVICE_LIST "<SYCL,SYCLSharedUSMSpace>")
  endif()

  if(Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)))
    message(WARNING
      "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with SYCL and complex double enabled.")
  endif()
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  kokkoskernels_add_option("INST_EXECSPACE_OPENMPTARGET" ${KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the execution space Kokkos::Experimental::OpenMPTarget. Disabling this when Kokkos_ENABLE_OPENMPTARGET is enabled may increase build times. Default: ON if Kokkos is OpenMPTarget-enabled, OFF otherwise.")

  kokkoskernels_add_option("INST_MEMSPACE_OPENMPTARGETSPACE" ${KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET_DEFAULT} BOOL
    "Whether to pre instantiate kernels for the memory space Kokkos::Experimental::OpenMPTargetSpace.  Disabling this when Kokkos_ENABLE_OPENMPTARGET is enabled may increase build times. Default: ON if Kokkos is OpenMPTarget-enabled, OFF otherwise.")

  if(KOKKOSKERNELS_INST_EXECSPACE_OPENMPTARGET AND KOKKOSKERNELS_INST_MEMSPACE_OPENMPTARGETSPACE)
    list(APPEND DEVICE_LIST "<OpenMPTarget,OpenMPTargetSpace>")
  endif()

  if(Trilinos_ENABLE_COMPLEX_DOUBLE AND ((NOT DEFINED CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS) OR (NOT CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS)))
    message(WARNING
      "The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either undefined or OFF.  Please set CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with OpenMPTarget and complex double enabled.")
  endif()
endif()

kokkoskernels_add_option("INST_MEMSPACE_HOSTSPACE" ${KOKKOSKERNELS_ADD_DEFAULT_ETI} BOOL
  "Whether to pre instantiate kernels for the memory space Kokkos::HostSpace.  Disabling this when one of the Host execution spaces is enabled may increase build times. Default: ON")

kokkoskernels_add_option("INST_EXECSPACE_OPENMP" ${KOKKOSKERNELS_INST_EXECSPACE_OPENMP_DEFAULT} BOOL
  "Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  Disabling this when Kokkos_ENABLE_OPENMP is enabled may increase build times. Default: ON if Kokkos is OpenMP-enabled, OFF otherwise.")

if(KOKKOSKERNELS_INST_EXECSPACE_OPENMP AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  list(APPEND DEVICE_LIST "<OpenMP,HostSpace>")
  if(NOT KOKKOS_ENABLE_OPENMP)
    message(FATAL_ERROR "Set ETI on for OPENMP, but Kokkos was not configured with the OPENMP backend")
  endif()
endif()

kokkoskernels_add_option("INST_EXECSPACE_THREADS" ${KOKKOSKERNELS_INST_EXECSPACE_THREADS_DEFAULT} BOOL
  "Whether to build kernels for the execution space Kokkos::Threads.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_THREADS is enabled may increase build times. Default: ON if Kokkos is Threads-enabled, OFF otherwise.")

#There continues to be name ambiguity with threads vs pthreads
set(KOKKOSKERNELS_INST_EXECSPACE_THREADS ${KOKKOSKERNELS_INST_EXECSPACE_THREADS})

if(KOKKOSKERNELS_INST_EXECSPACE_THREADS AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  list(APPEND DEVICE_LIST "<Threads,HostSpace>")
  if(NOT KOKKOS_ENABLE_THREADS)
    message(FATAL_ERROR "Set ETI on for THREADS, but Kokkos was not configured with the THREADS backend")
  endif()
endif()

kokkoskernels_add_option("INST_EXECSPACE_SERIAL" ${KOKKOSKERNELS_INST_EXECSPACE_SERIAL_DEFAULT} BOOL
  "Whether to build kernels for the execution space Kokkos::Serial.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_SERIAL is enabled may increase build times. Default: ON when Kokkos is Serial-enabled, OFF otherwise.")

set(EXECSPACE_CUDA_VALID_MEM_SPACES          CUDASPACE CUDAUVMSPACE)
set(EXECSPACE_HIP_VALID_MEM_SPACES           HIPSPACE HIPMANAGEDSPACE)
set(EXECSPACE_SYCL_VALID_MEM_SPACES          SYCLSPACE SYCLSHAREDSPACE)
set(EXECSPACE_OPENMPTARGET_VALID_MEM_SPACES  OPENMPTARGETSPACE)
set(EXECSPACE_SERIAL_VALID_MEM_SPACES        HOSTSPACE)
set(EXECSPACE_OPENMP_VALID_MEM_SPACES        HOSTSPACE)
set(EXECSPACE_THREADS_VALID_MEM_SPACES       HOSTSPACE)
set(DEVICES)
foreach(EXEC ${EXEC_SPACES})
  if(KOKKOSKERNELS_INST_${EXEC})
    foreach(MEM ${${EXEC}_VALID_MEM_SPACES})
      if(KOKKOSKERNELS_INST_MEMSPACE_${MEM})
        list(APPEND DEVICES ${EXEC}_MEMSPACE_${MEM})
        set(${EXEC}_MEMSPACE_${MEM}_CPP_TYPE "${${EXEC}_CPP_TYPE},${MEMSPACE_${MEM}_CPP_TYPE}")
        set(KOKKOSKERNELS_INST_${EXEC}_MEMSPACE_${MEM} ON)
      endif()
    endforeach()
  endif()
endforeach()

if(KOKKOSKERNELS_INST_EXECSPACE_SERIAL AND KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE)
  list(APPEND DEVICE_LIST "<Serial,HostSpace>")
  if(NOT KOKKOS_ENABLE_SERIAL)
    message(FATAL_ERROR "Set ETI on for SERIAL, but Kokkos was not configured with the SERIAL backend")
  endif()
endif()
