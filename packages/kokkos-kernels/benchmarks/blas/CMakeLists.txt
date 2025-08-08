KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_BENCHMARK(
  PerformanceTest_Benchmark
  SOURCES
    KokkosBlas1_dot_benchmark.cpp
    KokkosBlas1_dot_mv_benchmark.cpp
    KokkosBlas1_team_dot_benchmark.cpp
    BenchmarkMain.cpp
)

KOKKOSKERNELS_ADD_BENCHMARK(
    Blas2_gemv_Benchmark
    SOURCES KokkosBlas2_gemv_benchmark.cpp
)

KOKKOSKERNELS_ADD_BENCHMARK(
    Blas2_ger_Benchmark
    SOURCES KokkosBlas2_ger_benchmark.cpp
)

KOKKOSKERNELS_ADD_BENCHMARK(
    Blas3_gemm_benchmark
    SOURCES KokkosBlas3_gemm_benchmark.cpp
)
