
INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})

SET(HEADERS ${HEADERS}
 ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
    Kokkos_DeviceCuda_DeepCopy.hpp
    Kokkos_DeviceCuda_MDArrayView.hpp
    Kokkos_DeviceCuda_MultiVectorView.hpp
    Kokkos_DeviceCuda_ParallelDriver.hpp
    Kokkos_DeviceCuda_ParallelFor.hpp
    Kokkos_DeviceCuda_ParallelReduce.hpp
    Kokkos_DeviceCuda_ValueView.hpp
  )

APPEND_SET(SOURCES
    Kokkos_DeviceCuda.cu
  )

IF(Kokkos_ENABLE_CUDA)
  #set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  #CUDA_COMPILE(CUDA_FILE Kokkos_DeviceCuda.cu)
  PACKAGE_ADD_CUDA_LIBRARY(
    kokkosdevicecuda_lib
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
  )
  SET(CUDA_LIB kokkosdevicecuda_lib)
ELSE()
  SET(CUDA_LIB "")
ENDIF()


PACKAGE_ADD_LIBRARY(
  kokkosdevicecuda
  HEADERS ${HEADERS}
  DEPLIBS ${CUDA_LIB}
  )


# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()
