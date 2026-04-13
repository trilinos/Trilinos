kokkoskernels_include_directories(${CMAKE_BINARY_DIR}/test_common)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test_common)

kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

#####################
#                   #
# Define unit-tests #
#                   #
#####################

set(KK_ENABLE_GRAPH_TESTS ON)

if(KokkosKernels_TEST_ETI_ONLY)
  if(NOT KokkosKernels_INST_DOUBLE AND NOT KokkosKernels_INST_FLOAT)
    message(WARNING
      "Because only ETI'd type combinations are enabled for testing, the Kokkos Kernels graph tests require that double or float is enabled in ETI.")
    set(KK_ENABLE_GRAPH_TESTS OFF)
  endif()
endif()

set(TEST_MAIN ${PACKAGE_SOURCE_DIR}/test_common/Test_Main.cpp)
if(KK_ENABLE_GRAPH_TESTS)

  #####################
  #                   #
  # Add GPU backends  #
  #                   #
  #####################
  if(KOKKOS_ENABLE_CUDA)
    kokkoskernels_add_unit_test(graph_cuda
      SOURCES ${TEST_MAIN} backends/Test_Cuda_Graph.cpp
      COMPONENTS graph)
  endif()

  if(KOKKOS_ENABLE_HIP)
    kokkoskernels_add_unit_test(graph_hip
      SOURCES ${TEST_MAIN} backends/Test_HIP_Graph.cpp
      COMPONENTS graph)
  endif()

  if(KOKKOS_ENABLE_SYCL)
    kokkoskernels_add_unit_test(graph_sycl
      SOURCES ${TEST_MAIN} backends/Test_SYCL_Graph.cpp
      COMPONENTS graph)
  endif()

  if(KOKKOS_ENABLE_OPENMPTARGET)
    # KOKKOSKERNELS_ADD_UNIT_TEST(graph_openmptarget
    #   SOURCES ${TEST_MAIN} backends/Test_OpenMPTarget_Graph.cpp
    #   COMPONENTS graph)
  endif()

  #####################
  #                   #
  # Add CPU backends  #
  #                   #
  #####################
  if(KOKKOS_ENABLE_SERIAL)
    kokkoskernels_add_unit_test(graph_serial
      SOURCES ${TEST_MAIN} backends/Test_Serial_Graph.cpp
      COMPONENTS graph)
  endif()

  if(KOKKOS_ENABLE_OPENMP)
    kokkoskernels_add_unit_test(graph_openmp
      SOURCES ${TEST_MAIN} backends/Test_OpenMP_Graph.cpp
      COMPONENTS graph)
  endif()

  if(KOKKOS_ENABLE_THREADS)
    kokkoskernels_add_unit_test(graph_threads
      SOURCES ${TEST_MAIN} backends/Test_Threads_Graph.cpp
      COMPONENTS graph)
  endif()
endif()
