#####################
#                   #
# Add include files #
#                   #
#####################

kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/test_common)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/src)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/ode)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/ode)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/ode/src)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/ode/src)
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/ode/impl)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/ode/impl)

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
  kokkoskernels_add_unit_test(ode_cuda
    SOURCES ${TEST_MAIN} backends/Test_Cuda_ODE.cpp
    COMPONENTS ode)
endif()

if(KOKKOS_ENABLE_HIP)
  kokkoskernels_add_unit_test(ode_hip
    SOURCES ${TEST_MAIN} backends/Test_HIP_ODE.cpp
    COMPONENTS ode)
endif()

if(KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_unit_test(ode_sycl
    SOURCES ${TEST_MAIN} backends/Test_SYCL_ODE.cpp
    COMPONENTS ode)
endif()

if(KOKKOS_ENABLE_OPENMPTARGET)
  # kokkoskernels_add_unit_test(ode_openmptarget
  #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_ODE.cpp
  #   COMPONENTS ode)
endif()

#####################
#                   #
# Add CPU backends  #
#                   #
#####################
if(KOKKOS_ENABLE_SERIAL)
  kokkoskernels_add_unit_test(ode_serial
    SOURCES ${TEST_MAIN} backends/Test_Serial_ODE.cpp
    COMPONENTS ode)
endif()

if(KOKKOS_ENABLE_OPENMP)
  kokkoskernels_add_unit_test(ode_openmp
    SOURCES ${TEST_MAIN} backends/Test_OpenMP_ODE.cpp
    COMPONENTS ode)
endif()

if(KOKKOS_ENABLE_THREADS)
  kokkoskernels_add_unit_test(ode_threads
    SOURCES ${TEST_MAIN} backends/Test_Threads_ODE.cpp
    COMPONENTS ode)
endif()
