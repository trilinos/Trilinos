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
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/src)
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/impl)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/impl)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  message("batched enabled and blas not enabled, we need to include some include directories manually!")
  kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/src/blas)
  kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/blas)
  kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/src/blas/impl)
  kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/blas/impl)
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
    batched_sla_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(
    batched_sla_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(
    batched_sla_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_sla_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Batched_Sparse.cpp
  #   COMPONENTS batched_sla)
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(
    batched_sla_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(
    batched_sla_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(
    batched_sla_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Batched_Sparse.cpp
    COMPONENTS batched_sla)
endif()
