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
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/src)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/batched/dense/impl)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/impl)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF (NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  MESSAGE("batched enabled and blas not enabled, we need to include some include directories manually!")
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/blas)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/blas)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/blas/impl)
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/blas/impl)
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
    batched_sla_cuda
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Cuda_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

IF (KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_hip
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_HIP_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

IF (KOKKOS_ENABLE_SYCL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_sycl
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_SYCL_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(
  #   batched_sla_openmptarget
  #   SOURCES
  #     ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
  #     backends/Test_OpenMPTarget_Batched_Sparse.cpp
  #   COMPONENTS batched_sla
  # )
ENDIF ()



#####################
#                   #
# Add CPU backends  #
#                   #
#####################
IF (KOKKOS_ENABLE_SERIAL)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_serial
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Serial_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

IF (KOKKOS_ENABLE_OPENMP)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_openmp
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_OpenMP_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

IF (KOKKOS_ENABLE_THREADS)
  KOKKOSKERNELS_ADD_UNIT_TEST(
    batched_sla_threads
    SOURCES
      ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp
      backends/Test_Threads_Batched_Sparse.cpp
    COMPONENTS batched_sla
  )
ENDIF ()

