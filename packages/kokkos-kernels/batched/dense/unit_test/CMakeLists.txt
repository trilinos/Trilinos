#####################
#                   #
# Add include files #
#                   #
#####################

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/test_common)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/src)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/impl)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF (NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/blas/src)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/blas/src)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/blas/impl)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/blas/impl)
ENDIF()

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
    batched_dla_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

IF (KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

IF (KOKKOS_ENABLE_SYCL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_dla_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_Batched_Dense.cpp
  #   COMPONENTS batched_dla
  # )

  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_gemm_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_Batched_GEMM.cpp
  #   COMPONENTS batched_dla
  # )
ENDIF ()



#####################
#                   #
# Add CPU backends  #
#                   #
#####################
IF (KOKKOS_ENABLE_SERIAL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

IF (KOKKOS_ENABLE_THREADS)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_dla_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_Batched_Dense.cpp
    COMPONENTS batched_dla
  )

  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_gemm_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_Batched_GEMM.cpp
    COMPONENTS batched_dla
  )
ENDIF ()

