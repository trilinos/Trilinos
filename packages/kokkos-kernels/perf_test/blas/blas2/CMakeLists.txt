KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_EXECUTABLE(
    KokkosBlas2_gemv_perf_test
    SOURCES KokkosBlas2_gemv_perf_test.cpp
)

IF(KokkosKernels_ENABLE_BENCHMARK)
    KOKKOSKERNELS_ADD_BENCHMARK(
        Blas2_gemv_Benchmark
        SOURCES KokkosBlas2_gemv_perf_test_benchmark.cpp
    )
ENDIF()

IF(KokkosKernels_ENABLE_BENCHMARK)
    KOKKOSKERNELS_ADD_BENCHMARK(
        Blas2_ger_Benchmark
        SOURCES KokkosBlas2_ger_perf_test_benchmark.cpp
    )
ENDIF()
