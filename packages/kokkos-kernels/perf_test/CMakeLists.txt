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
    ADD_SUBDIRECTORY(performance)
    #ADD_SUBDIRECTORY(common)

endif()

IF(KokkosKernels_ENABLE_BENCHMARK)

    IF (KOKKOSKERNELS_HAS_TRILINOS)
      message(FATAL_ERROR "Benchmarks are not supported when building as part of Trilinos")
    ENDIF()

    find_package(benchmark QUIET)

    IF(benchmark_FOUND)
      MESSAGE(STATUS "Using google benchmark found in ${benchmark_DIR}")
    ELSE()
      message(STATUS "No installed google benchmark found, fetching from GitHub")
      include(FetchContent)
      SET(BENCHMARK_ENABLE_TESTING OFF)

      list(APPEND CMAKE_MESSAGE_INDENT "    ")
      #Note: recent bug (google/benchmark#1441) is preventing us from using
      # the latest benchmark release.
      FetchContent_Declare(
        googlebenchmark
        URL https://github.com/google/benchmark/archive/refs/tags/v1.6.2.tar.gz
        URL_HASH MD5=14d14849e075af116143a161bc3b927b
      )
      FetchContent_MakeAvailable(googlebenchmark)
      list(POP_BACK CMAKE_MESSAGE_INDENT)

      include_directories(${benchmark_SOURCE_DIR}/include)

      # Suppress clang-tidy diagnostics on code that we do not have control over
      IF(CMAKE_CXX_CLANG_TIDY)
        SET_TARGET_PROPERTIES(benchmark PROPERTIES CXX_CLANG_TIDY "")
      ENDIF()

      target_compile_options(benchmark PRIVATE -w)
      target_compile_options(benchmark_main PRIVATE -w)
    ENDIF()

    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    FUNCTION(KOKKOSKERNELS_ADD_BENCHMARK NAME)
      CMAKE_PARSE_ARGUMENTS(
        BENCHMARK
        ""
        ""
        "SOURCES"
        ${ARGN}
      )
      IF(DEFINED BENCHMARK_UNPARSED_ARGUMENTS)
        MESSAGE(
          WARNING
          "Unexpected arguments when adding a benchmark: "
          ${BENCHMARK_UNPARSED_ARGUMENTS}
        )
      ENDIF()

      SET(BENCHMARK_NAME ${PACKAGE_NAME}_${NAME})

      ADD_EXECUTABLE(
        ${BENCHMARK_NAME}
        ${BENCHMARK_SOURCES}
      )
      TARGET_LINK_LIBRARIES(
        ${BENCHMARK_NAME}
        PRIVATE benchmark::benchmark Kokkos::kokkoskernels
      )
      FOREACH(SOURCE_FILE ${BENCHMARK_SOURCES})
        SET_SOURCE_FILES_PROPERTIES(
          ${SOURCE_FILE}
          PROPERTIES LANGUAGE CXX
        )
      ENDFOREACH()

      STRING(TIMESTAMP BENCHMARK_TIME "%Y-%m-%d_T%H-%M-%S" UTC)
      SET(
        BENCHMARK_ARGS
        --benchmark_counters_tabular=true
        --benchmark_out=${BENCHMARK_NAME}_${BENCHMARK_TIME}.json
      )

      ADD_TEST(
        NAME ${BENCHMARK_NAME}
        COMMAND ${BENCHMARK_NAME} ${BENCHMARK_ARGS}
      )
    ENDFUNCTION()

    SET(
      BENCHMARK_SOURCES
      BenchmarkMain.cpp
    )

    KOKKOSKERNELS_ADD_BENCHMARK(
      PerformanceTest_Benchmark
      SOURCES ${BENCHMARK_SOURCES}
    )

endif()
