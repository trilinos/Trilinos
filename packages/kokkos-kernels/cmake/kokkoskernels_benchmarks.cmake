IF(KOKKOSKERNELS_HAS_TRILINOS)
    MESSAGE(
        FATAL_ERROR
        "Benchmarks are not supported when building as part of Trilinos")
ENDIF()

FIND_PACKAGE(benchmark QUIET)

IF(benchmark_FOUND)
    MESSAGE(STATUS "Using google benchmark found in ${benchmark_DIR}")
ELSE()
    MESSAGE(STATUS "No installed google benchmark found, fetching from GitHub")
    INCLUDE(FetchContent)
    SET(BENCHMARK_ENABLE_TESTING OFF)

    LIST(APPEND CMAKE_MESSAGE_INDENT "[benchmark] ")

    # Note: recent bug (google/benchmark#1441) is preventing us from using
    # the latest benchmark release.
    SET(BENCHMARK_VERSION 1.6.2)

    # CMake 3.24 introduced DOWNLOAD_EXTRACT_TIMESTAMP, which controls whether
    # extracting this file archive sets the file times to archive time (TRUE),
    # or to extraction time (FALSE).
    # In CMake 3.24+, the default is FALSE
    # Prior, it did not exist, and was effectively TRUE
    # Here, we okay the new default to silence CMP0135 warning
    IF (${CMAKE_VERSION} VERSION_LESS "3.24.0")
        FetchContent_Declare(
            googlebenchmark
            URL https://github.com/google/benchmark/archive/refs/tags/v${BENCHMARK_VERSION}.tar.gz
            URL_HASH MD5=14d14849e075af116143a161bc3b927b
        )
    ELSE()
        FetchContent_Declare(
            googlebenchmark
            URL https://github.com/google/benchmark/archive/refs/tags/v${BENCHMARK_VERSION}.tar.gz
            URL_HASH MD5=14d14849e075af116143a161bc3b927b
            DOWNLOAD_EXTRACT_TIMESTAMP FALSE
        )
    ENDIF()
    FetchContent_MakeAvailable(googlebenchmark)
    LIST(POP_BACK CMAKE_MESSAGE_INDENT)

    # remove the CXX_CLANG_TIDY property from google benchmark
    # when we're running clang-tidy we don't care if google benchmark passes or not
    SET_PROPERTY(TARGET benchmark PROPERTY CXX_CLANG_TIDY "")
    SET_PROPERTY(TARGET benchmark_main PROPERTY CXX_CLANG_TIDY "")

    # disable warnings for google benchmark
    TARGET_COMPILE_OPTIONS(benchmark PRIVATE -w)
    TARGET_COMPILE_OPTIONS(benchmark_main PRIVATE -w)
ENDIF()

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
    TARGET_INCLUDE_DIRECTORIES(
        ${BENCHMARK_NAME}
        SYSTEM PRIVATE ${benchmark_SOURCE_DIR}/include
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

    SET_PROPERTY(TEST ${BENCHMARK_NAME} PROPERTY LABELS Benchmark)

    IF(NOT KokkosKernels_RUN_BENCHMARKS)
        SET_PROPERTY(TEST ${BENCHMARK_NAME} PROPERTY DISABLED TRUE)
    ENDIF()
ENDFUNCTION()
