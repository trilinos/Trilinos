kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

if( (KOKKOS_ENABLE_CUDA AND KOKKOSKERNELS_ENABLE_TPL_CUSOLVER) OR
    (KOKKOS_ENABLE_HIP AND KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER) OR
    ((KOKKOS_ENABLE_OPENMP OR KOKKOS_ENABLE_SERIAL OR KOKKOS_ENABLE_THREADS) AND
    (KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)) )
  kokkoskernels_add_executable_and_test(docs_geqrf SOURCES KokkosLapack_docs_geqrf.cpp)
  kokkoskernels_add_executable_and_test(docs_gemqr SOURCES KokkosLapack_docs_gemqr.cpp)
  kokkoskernels_add_executable_and_test(docs_gegqr SOURCES KokkosLapack_docs_gegqr.cpp)
  kokkoskernels_add_executable_and_test(docs_potr  SOURCES KokkosLapack_docs_potr.cpp)
endif()
