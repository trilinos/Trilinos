if(KOKKOSKERNELS_HAS_TRILINOS)
  message(FATAL_ERROR "Benchmarks are not supported when building as part of Trilinos")
endif()

find_package(benchmark QUIET)

if(benchmark_FOUND)
  message(STATUS "Using google benchmark found in ${benchmark_DIR}")
else()
  message(STATUS "No installed google benchmark found, fetching from GitHub")
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF)

  list(APPEND CMAKE_MESSAGE_INDENT "[benchmark] ")

  # Note: recent bug (google/benchmark#1441) is preventing us from using
  # the latest benchmark release.
  set(BENCHMARK_VERSION 1.6.2)

  # CMake 3.24 introduced DOWNLOAD_EXTRACT_TIMESTAMP, which controls whether
  # extracting this file archive sets the file times to archive time (TRUE),
  # or to extraction time (FALSE).
  # In CMake 3.24+, the default is FALSE
  # Prior, it did not exist, and was effectively TRUE
  # Here, we okay the new default to silence CMP0135 warning
  if(${CMAKE_VERSION} VERSION_LESS "3.24.0")
    FetchContent_Declare(googlebenchmark
      URL https://github.com/google/benchmark/archive/refs/tags/v${BENCHMARK_VERSION}.tar.gz
      URL_HASH MD5=14d14849e075af116143a161bc3b927b)
  else()
    FetchContent_Declare(googlebenchmark
      URL https://github.com/google/benchmark/archive/refs/tags/v${BENCHMARK_VERSION}.tar.gz
      URL_HASH MD5=14d14849e075af116143a161bc3b927b
      DOWNLOAD_EXTRACT_TIMESTAMP FALSE)
  endif()
  FetchContent_MakeAvailable(googlebenchmark)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # remove the CXX_CLANG_TIDY property from google benchmark
  # when we're running clang-tidy we don't care if google benchmark passes or not
  set_property(TARGET benchmark PROPERTY CXX_CLANG_TIDY "")
  set_property(TARGET benchmark_main PROPERTY CXX_CLANG_TIDY "")

  # disable warnings for google benchmark
  target_compile_options(benchmark PRIVATE -w)
  target_compile_options(benchmark_main PRIVATE -w)
endif()

function(kokkoskernels_add_benchmark NAME)
  cmake_parse_arguments(BENCHMARK "" "" "SOURCES" ${ARGN})

  if(DEFINED BENCHMARK_UNPARSED_ARGUMENTS)
    message(WARNING "Unexpected arguments when adding a benchmark: " ${BENCHMARK_UNPARSED_ARGUMENTS})
  endif()

  set(BENCHMARK_NAME ${PACKAGE_NAME}_${NAME})

  add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCES})
  target_link_libraries(${BENCHMARK_NAME} PRIVATE benchmark::benchmark Kokkos::kokkoskernels)
  target_include_directories(${BENCHMARK_NAME} SYSTEM PRIVATE ${benchmark_SOURCE_DIR}/include)

  foreach(SOURCE_FILE ${BENCHMARK_SOURCES})
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CXX)
  endforeach()

  string(TIMESTAMP BENCHMARK_TIME "%Y-%m-%d_T%H-%M-%S" UTC)
  set(BENCHMARK_ARGS --benchmark_counters_tabular=true --benchmark_out=${BENCHMARK_NAME}_${BENCHMARK_TIME}.json)

  add_test(NAME ${BENCHMARK_NAME} COMMAND ${BENCHMARK_NAME} ${BENCHMARK_ARGS})

  set_property(TEST ${BENCHMARK_NAME} PROPERTY LABELS Benchmark)

  if(NOT KokkosKernels_RUN_BENCHMARKS)
    set_property(TEST ${BENCHMARK_NAME} PROPERTY DISABLED TRUE)
  endif()
endfunction()
