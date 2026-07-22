#####################
#                   #
# Add include files #
#                   #
#####################

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
  kokkoskernels_add_unit_test(common_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_Common.cpp
    COMPONENTS common)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(common_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_Common.cpp
    COMPONENTS common)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(common_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_Common.cpp
    COMPONENTS common)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # KOKKOSKERNELS_ADD_UNIT_TEST(common_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Common.cpp
  #   COMPONENTS common
  # )
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(common_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_Common.cpp
    COMPONENTS common)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(common_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_Common.cpp
    COMPONENTS common)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(common_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_Common.cpp
    COMPONENTS common)
endif()

# Add eager_initialize test, which is not backend-specific
kokkoskernels_add_unit_test(common_eager_initialize
  SOURCES Test_Common_EagerInitialize.cpp
  COMPONENTS common)

# Testing singletons with atexit(Kokkos::finalize)
# which is also not backend specific
kokkoskernels_add_unit_test(common_atexit_finalize
  SOURCES Test_Common_atexit_finalize.cpp
  COMPONENTS common)
