KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_BENCHMARK(
  sparse_par_ilut
  SOURCES KokkosSparse_par_ilut.cpp
)

# Provide -DGinkgo_DIR to cmake to enable the ginkgo test in sparse_par_ilut. Ginkgo_DIR should
# point to the dir in the ginkgo install area that contains the GinkgoConfig.cmake file.
# For me, this was $gingko_install_dir/lib64/cmake/Ginkgo
if (Ginkgo_DIR)
  find_package(Ginkgo REQUIRED)

  target_compile_definitions(KokkosKernels_sparse_par_ilut PRIVATE "USE_GINKGO")
  target_link_libraries(KokkosKernels_sparse_par_ilut PRIVATE Ginkgo::ginkgo)
endif()

KOKKOSKERNELS_ADD_BENCHMARK(
  sparse_spmv_benchmark SOURCES KokkosSparse_spmv_benchmark.cpp
)

KOKKOSKERNELS_ADD_BENCHMARK(
  sparse_spmv_bsr_benchmark SOURCES KokkosSparse_spmv_bsr_benchmark.cpp
)

# hipcc 5.2 has an underlying clang that has the std::filesystem
# in an experimental namespace and a different library
if (Kokkos_CXX_COMPILER_ID STREQUAL HIPCC AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 6.3)
  target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lstdc++fs)
# IntelLLVM < 2023.1.0 (and possible higher versions too) have an underlying clang that has the std::filesystem
# in an experimental namespace and a different library
elseif (Kokkos_CXX_COMPILER_ID STREQUAL IntelLLVM AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 2023.1.0)
  target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lstdc++fs)
# GNU compilers < 9.1 require linking against -lstdc++fs
elseif (Kokkos_CXX_COMPILER_ID STREQUAL GNU AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 9.1)
  target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lstdc++fs)
# Clang versions < 9.0 require linking against -lc++fs
elseif (Kokkos_CXX_COMPILER_ID STREQUAL Clang AND Kokkos_CXX_COMPILER_VERSION VERSION_LESS 9.0)
  target_link_libraries(KokkosKernels_sparse_spmv_bsr_benchmark PRIVATE -lc++fs)
endif()

