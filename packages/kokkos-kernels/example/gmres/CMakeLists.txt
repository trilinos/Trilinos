kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_add_executable(gmres_ex_real_A SOURCES ex_real_A.cpp)

# FIXME_SYCL CUDA_ERROR_INVALID_ADDRESS_SPACE
if(NOT KOKKOS_ENABLE_SYCL)
  kokkoskernels_add_executable_and_test(gmres_test_prec SOURCES test_prec.cpp)
endif()
