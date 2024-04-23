if (KokkosKernels_ENABLE_PERFTESTS)
    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

    #At some point, we may wish to make this into real "performance
    #tests, " in the sense that they can pass or fail.  At that point, use
    #"CATEGORIES PERFORMANCE" to mark them as such.For now, we just
    #build the executables for manual use, but don't run the tests.  They
    #build correctly with or without MPI, but only run them with a single
    #MPI process.

    SET(GTEST_SOURCE_DIR ${PACKAGE_SOURCE_DIR}/tpls/gtest)

    KOKKOSKERNELS_ADD_TEST_LIBRARY(
                                kokkoskernelsperf_gtest
                                HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
                                SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
                                    )
    #Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.

    TARGET_COMPILE_DEFINITIONS(kokkoskernelsperf_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
    TARGET_INCLUDE_DIRECTORIES(kokkoskernelsperf_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)

    #Gtest minimally requires C++ 11
    TARGET_COMPILE_FEATURES(kokkoskernelsperf_gtest PUBLIC cxx_std_11)

    KOKKOSKERNELS_INCLUDE_DIRECTORIES(sparse)

    if(KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
        #Add RPS implementations of KK perf tests here
        KOKKOSKERNELS_ADD_EXECUTABLE(
            tracked_testing
            SOURCES KokkosKernelsTrackedTesting.cpp
            sparse/KokkosSparse_spmv_test.cpp
            blas/blas2/KokkosBlas2_gemv_tracked_perf_test.cpp
            blas/blas1/KokkosBlas_dot_tracked_perf_test.cpp
            blas/blas1/KokkosBlas_team_dot_tracked_perf_test.cpp
            blas/blas3/KokkosBlas3_gemm_tracked_perf_test.cpp
            PerfTestUtilities.cpp
            sparse/spmv/OpenMPSmartStatic_SPMV.cpp
            #sparse / KokkosSparse_spgemm_test.cpp
            )
    endif()

    ADD_COMPONENT_SUBDIRECTORY(batched)
    ADD_COMPONENT_SUBDIRECTORY(graph)
    ADD_COMPONENT_SUBDIRECTORY(sparse)
    ADD_COMPONENT_SUBDIRECTORY(blas)
    ADD_COMPONENT_SUBDIRECTORY(ode)
    ADD_COMPONENT_SUBDIRECTORY(lapack)
    ADD_SUBDIRECTORY(performance)
    #ADD_SUBDIRECTORY(common)

endif()

if(KokkosKernels_ENABLE_BENCHMARK)
    KOKKOSKERNELS_ADD_BENCHMARK(
      PerformanceTest_Benchmark
      SOURCES
        blas/blas1/KokkosBlas_dot_perf_test_benchmark.cpp
        blas/blas1/KokkosBlas_dot_mv_perf_test_benchmark.cpp
        blas/blas1/KokkosBlas_team_dot_perf_test_benchmark.cpp
        BenchmarkMain.cpp
    )
endif()
