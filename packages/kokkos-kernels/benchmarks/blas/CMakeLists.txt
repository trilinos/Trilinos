kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_add_benchmark(
  PerformanceTest_Benchmark
  SOURCES
    KokkosBlas1_dot_benchmark.cpp
    KokkosBlas1_dot_mv_benchmark.cpp
    KokkosBlas1_team_dot_benchmark.cpp
    BenchmarkMain.cpp)

kokkoskernels_add_benchmark(Blas2_gemv_Benchmark SOURCES KokkosBlas2_gemv_benchmark.cpp)
kokkoskernels_add_benchmark(Blas2_ger_Benchmark  SOURCES KokkosBlas2_ger_benchmark.cpp)
kokkoskernels_add_benchmark(Blas3_gemm_benchmark SOURCES KokkosBlas3_gemm_benchmark.cpp)
