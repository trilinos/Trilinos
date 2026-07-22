#####################
#                   #
# Add include files #
#                   #
#####################

kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/test_common)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/src)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/batched)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/src)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/impl)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/blas/src)
  kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/blas/src)
  kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/blas/impl)
  kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/blas/impl)
endif()

#####################
#                   #
# Define unit-tests #
#                   #
#####################

#####################
#                   #
# Add GPU backends  #
#                   #
#####################
set(TEST_MAIN ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp)
if(KOKKOS_ENABLE_CUDA)
  kokkoskernels_add_unit_test(
    batched_dla_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(
    batched_dla_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(
    batched_dla_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_dla_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Batched_Dense.cpp
  #   COMPONENTS batched_dla
  # )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_gemm_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Batched_GEMM.cpp
  #   COMPONENTS batched_dla
  # )
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(
    batched_dla_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(
    batched_dla_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(
    batched_dla_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Batched_Dense.cpp
    COMPONENTS batched_dla)

  kokkoskernels_add_unit_test(
    batched_gemm_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Batched_GEMM.cpp
    COMPONENTS batched_dla)
endif()
