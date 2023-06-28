KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_EXECUTABLE(
    KokkosBlas3_perf_test
    SOURCES KokkosBlas3_perf_test.cpp
    TESTONLYLIBS kokkoskernelsperf_gtest
)

KOKKOSKERNELS_ADD_EXECUTABLE(
    KokkosBlas3_gemm_perf_test
    SOURCES KokkosBlas3_gemm_standalone_perf_test.cpp
)

IF(KokkosKernels_ENABLE_BENCHMARK)
    KOKKOSKERNELS_ADD_BENCHMARK(
        Blas3_gemm_benchmark
        SOURCES KokkosBlas3_gemm_standalone_perf_test_benchmark.cpp
    )
ENDIF()
