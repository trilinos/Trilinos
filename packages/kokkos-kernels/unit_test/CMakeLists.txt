
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../test_common)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/blas)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/blas)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/sparse)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/sparse)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/graph)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/graph)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/dense)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/impl)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/sparse)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/sparse/impl)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/impl)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/common)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/common)

# Kyungjoo: when kokkoskernels and kokkos are symbolic link and their src directories are overrided
#           the following relative path does not work or users should put kokkoskernels and kokkos
#           at the same place
SET(GTEST_SOURCE_DIR ${PACKAGE_SOURCE_DIR}/tpls/gtest)

KOKKOSKERNELS_ADD_TEST_LIBRARY(
  kokkoskernels_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  )
# Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.
TARGET_COMPILE_DEFINITIONS(kokkoskernels_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
TARGET_INCLUDE_DIRECTORIES(kokkoskernels_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)

#Gtest minimally requires C++11
TARGET_COMPILE_FEATURES(kokkoskernels_gtest PUBLIC cxx_std_11)






IF (KOKKOS_ENABLE_CUDA)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/cuda)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/cuda)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Blas.cpp
    COMPONENTS blas
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Graph.cpp
    COMPONENTS graph
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_cuda
    SOURCES
      Test_Main.cpp
      cuda/Test_Cuda_Common.cpp
    COMPONENTS common
  )
ENDIF ()

IF (KOKKOS_ENABLE_SYCL)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/sycl)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/sycl)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Blas.cpp
    COMPONENTS blas
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Graph.cpp
    COMPONENTS graph
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_sycl
    SOURCES
      Test_Main.cpp
      sycl/Test_SYCL_Common.cpp
    COMPONENTS common
  )
ENDIF ()

IF (KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/hip)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/hip)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Blas.cpp
    COMPONENTS blas
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Graph.cpp
    COMPONENTS graph
  )
  
  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_hip
    SOURCES
      Test_Main.cpp
      hip/Test_HIP_Common.cpp
    COMPONENTS common
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/openmp_target)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/openmptarget)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_openmptarget
    SOURCES
      Test_Main.cpp
      openmptarget/Test_OpenMPTarget_Blas.cpp
    COMPONENTS blas
  )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_dla_openmptarget
  #   SOURCES
  #     Test_Main.cpp
  #     openmptarget/Test_OpenMPTarget_Batched_Dense.cpp
  #   COMPONENTS batched
  # )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   sparse_openmptarget
  #   SOURCES
  #     Test_Main.cpp
  #     openmptarget/Test_OpenMPTarget_Sparse.cpp
  #   COMPONENTS sparse
  # )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_openmptarget
    SOURCES
      Test_Main.cpp
      openmptarget/Test_OpenMPTarget_Graph.cpp
    COMPONENTS graph
  )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   common_openmptarget
  #   SOURCES
  #     Test_Main.cpp
  #     openmptarget/Test_OpenMPTarget_Common.cpp
  #   COMPONENTS common
  # )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMP)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/openmp)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/openmp)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Blas.cpp
    COMPONENTS blas
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Graph.cpp
    COMPONENTS graph
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_openmp
    SOURCES
      Test_Main.cpp
      openmp/Test_OpenMP_Common.cpp
    COMPONENTS common
  )
ENDIF ()

IF (KOKKOS_ENABLE_SERIAL)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/serial)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/serial)

  # IF (KOKKOS_ENABLE_DEBUG)
  #   SET(DISABLE_SLOW_DGEMM_DOUBLE_TEST "--gtest_filter=-serial.gemm_double")
  # ENDIF()

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Blas.cpp
    COMPONENTS blas
  )
    
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Graph.cpp
    COMPONENTS graph
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_serial
    SOURCES
      Test_Main.cpp
      serial/Test_Serial_Common.cpp
    COMPONENTS common
  )
ENDIF ()

IF (KOKKOS_ENABLE_THREADS)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/threads)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/threads)

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Blas.cpp
    COMPONENTS blas
  )
    
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Batched_Dense.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Batched_Sparse.cpp
    COMPONENTS batched
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    graph_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Graph.cpp
    COMPONENTS graph
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    common_threads
    SOURCES
      Test_Main.cpp
      threads/Test_Threads_Common.cpp
    COMPONENTS common
  )
ENDIF ()


 
