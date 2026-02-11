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
  kokkoskernels_add_unit_test(lapack_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Lapack.cpp
    COMPONENTS lapack)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(lapack_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Lapack.cpp
    COMPONENTS lapack)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(lapack_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Lapack.cpp
    COMPONENTS lapack)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # kokkoskernels_add_unit_test(lapack_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Lapack.cpp
  #   COMPONENTS lapack)
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(lapack_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Lapack.cpp
    COMPONENTS lapack)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(lapack_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Lapack.cpp
    COMPONENTS lapack)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(lapack_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Lapack.cpp
    COMPONENTS lapack)
endif()
