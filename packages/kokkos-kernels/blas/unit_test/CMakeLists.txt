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
    blas_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

IF (KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

IF (KOKKOS_ENABLE_SYCL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   blas_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_Blas.cpp
  #   COMPONENTS blas
  # )
ENDIF ()



#####################
#                   #
# Add CPU backends  #
#                   #
#####################
IF (KOKKOS_ENABLE_SERIAL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

IF (KOKKOS_ENABLE_THREADS)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    blas_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_Blas.cpp
    COMPONENTS blas
  )
ENDIF ()

