# This is not a Trilinos-standard CMakeLists.txt file. It should be updated.
PROJECT(MDARRAY_UNIT_TESTS)

# Trilinos sets this. Other packages should not.
# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../ThreadPool/src
  )

FIND_PACKAGE(CUDA)

INCLUDE(FindCUDA)

SET(HEADERS "")
SET(SOURCES "")

LIST(APPEND HEADERS
  Kokkos_test_hex_grad.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_MDArrayViewRawData.hpp
  )

LIST(APPEND SOURCES
  )


SET(CUDA_HEADERS "")
SET(CUDA_SOURCES "")

LIST(APPEND CUDA_HEADERS
  ${HEADERS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_CudaMDArrayView.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_CudaMap.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_CudaDevice.hpp
  )


LIST(APPEND CUDA_SOURCES
  ${SOURCES}
  Kokkos_test_mdarray_cuda.cu
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_CudaDevice.cu
  )

CUDA_ADD_EXECUTABLE(
  bin/cuda_test.exe
  ${CUDA_SOURCES}
  )


SET(HOST_HEADERS "")
SET(HOST_SOURCES "")

LIST(APPEND HOST_HEADERS
  ${HEADERS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_HostMDArrayView.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_HostMap.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_HostDevice.hpp
  )


LIST(APPEND HOST_SOURCES
  ${SOURCES}
  Kokkos_test_mdarray_host.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/Kokkos_HostDevice.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../ThreadPool/src/TPI.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../ThreadPool/src/TPI_Walltime.c
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../ThreadPool/src/ThreadPool_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/ThreadPool_config.h
  )

ADD_EXECUTABLE(
  bin/host_test.exe
  ${HOST_SOURCES}
  )

