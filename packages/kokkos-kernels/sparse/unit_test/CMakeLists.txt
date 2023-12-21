KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/test_common)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})


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
IF (KOKKOS_ENABLE_CUDA)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

IF (KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

IF (KOKKOS_ENABLE_SYCL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   sparse_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_Sparse.cpp
  #   COMPONENTS sparse
  # )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   blocksparse_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_BlockSparse.cpp
  #   COMPONENTS sparse
  # )
ENDIF ()



#####################
#                   #
# Add CPU backends  #
#                   #
#####################
IF (KOKKOS_ENABLE_SERIAL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

IF (KOKKOS_ENABLE_THREADS)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    sparse_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_Sparse.cpp
    COMPONENTS sparse
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    blocksparse_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_BlockSparse.cpp
    COMPONENTS sparse
  )
ENDIF ()

