kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/test_common)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

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
  kokkoskernels_add_unit_test(sparse_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_BlockSparse.cpp
    COMPONENTS sparse)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(sparse_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_BlockSparse.cpp
    COMPONENTS sparse)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(sparse_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_BlockSparse.cpp
    COMPONENTS sparse)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # kokkoskernels_add_unit_test(sparse_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Sparse.cpp
  #   COMPONENTS sparse)

  # kokkoskernels_add_unit_test(blocksparse_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_BlockSparse.cpp
  #   COMPONENTS sparse)
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(sparse_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_BlockSparse.cpp
    COMPONENTS sparse)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(sparse_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_BlockSparse.cpp
    COMPONENTS sparse)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(sparse_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Sparse.cpp
    COMPONENTS sparse)

  kokkoskernels_add_unit_test(blocksparse_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_BlockSparse.cpp
    COMPONENTS sparse)
endif()
