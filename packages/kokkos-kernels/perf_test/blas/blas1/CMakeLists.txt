kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_add_executable(KokkosBlas_perf_test          SOURCES KokkosBlas_perf_test.cpp)
