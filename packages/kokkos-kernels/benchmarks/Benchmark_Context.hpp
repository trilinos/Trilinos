/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
*/

#ifndef KOKKOSKERNELS_PERFTEST_BENCHMARK_CONTEXT_HPP
#define KOKKOSKERNELS_PERFTEST_BENCHMARK_CONTEXT_HPP

#include "KokkosKernels_PrintConfiguration.hpp"

#include <cstdlib>
#include <string>

#include <benchmark/benchmark.h>

#include <Kokkos_Core.hpp>
#include <KokkosKernels_PrintConfiguration.hpp>
#include <KokkosKernels_Version_Info.hpp>

namespace KokkosKernelsBenchmark {

/// \brief Remove unwanted spaces and colon signs from input string. In case of
/// invalid input it will return an empty string.
inline std::string remove_unwanted_characters(std::string str) {
  auto from = str.find_first_not_of(" :");
  auto to   = str.find_last_not_of(" :");

  if (from == std::string::npos || to == std::string::npos) {
    return "";
  }

  // return extracted part of string without unwanted spaces and colon signs
  return str.substr(from, to + 1);
}

/// \brief Extract all key:value pairs from kokkos configuration and add it to
/// the benchmark context
inline void add_kokkos_configuration(bool verbose) {
  std::ostringstream msg;
  Kokkos::print_configuration(msg, verbose);
  KokkosKernels::print_configuration(msg);

  // Iterate over lines returned from kokkos and extract key:value pairs
  std::stringstream ss{msg.str()};
  for (std::string line; std::getline(ss, line, '\n');) {
    auto found = line.find_first_of(':');
    if (found != std::string::npos) {
      auto val = remove_unwanted_characters(line.substr(found + 1));
      // Ignore line without value, for example a category name
      if (!val.empty()) {
        benchmark::AddCustomContext(remove_unwanted_characters(line.substr(0, found)), val);
      }
    }
  }
}

/// \brief Add Kokkos Kernels git info and google benchmark release to
/// benchmark context.
inline void add_version_info() {
  using namespace KokkosKernels::Impl;

  if (!GIT_BRANCH.empty()) {
    benchmark::AddCustomContext("GIT_BRANCH", std::string(GIT_BRANCH));
    benchmark::AddCustomContext("GIT_COMMIT_HASH", std::string(GIT_COMMIT_HASH));
    benchmark::AddCustomContext("GIT_CLEAN_STATUS", std::string(GIT_CLEAN_STATUS));
    benchmark::AddCustomContext("GIT_COMMIT_DESCRIPTION", std::string(GIT_COMMIT_DESCRIPTION));
    benchmark::AddCustomContext("GIT_COMMIT_DATE", std::string(GIT_COMMIT_DATE));
  }
  if (!BENCHMARK_VERSION.empty()) {
    benchmark::AddCustomContext("GOOGLE_BENCHMARK_VERSION", std::string(BENCHMARK_VERSION));
  }
}

inline void add_env_info() {
  auto num_threads = std::getenv("OMP_NUM_THREADS");
  if (num_threads) {
    benchmark::AddCustomContext("OMP_NUM_THREADS", num_threads);
  }
  auto dynamic = std::getenv("OMP_DYNAMIC");
  if (dynamic) {
    benchmark::AddCustomContext("OMP_DYNAMIC", dynamic);
  }
  auto proc_bind = std::getenv("OMP_PROC_BIND");
  if (proc_bind) {
    benchmark::AddCustomContext("OMP_PROC_BIND", proc_bind);
  }
  auto places = std::getenv("OMP_PLACES");
  if (places) {
    benchmark::AddCustomContext("OMP_PLACES", places);
  }
}

/// \brief Gather all context information and add it to benchmark context
inline void add_benchmark_context(bool verbose = false) {
  add_kokkos_configuration(verbose);
  add_version_info();
  add_env_info();
}

template <class FuncType, class... ArgsToCallOp>
inline auto register_benchmark(const char* name, FuncType func, std::vector<std::string> arg_names,
                               std::vector<int64_t> args, int repeat, ArgsToCallOp&&... func_args) {
  if (repeat > 0) {
    return benchmark::RegisterBenchmark(name, func, std::forward<ArgsToCallOp>(func_args)...)
        ->ArgNames(arg_names)
        ->Args(args)
        ->UseManualTime()
        ->Iterations(repeat);
  } else {
    return benchmark::RegisterBenchmark(name, func, std::forward<ArgsToCallOp>(func_args)...)
        ->ArgNames(arg_names)
        ->Args(args)
        ->UseManualTime();
  }
}

template <class FuncType, class... ArgsToCallOp>
inline auto register_benchmark_real_time(const char* name, FuncType func, std::vector<std::string> arg_names,
                                         std::vector<int64_t> args, int repeat, ArgsToCallOp&&... func_args) {
  if (repeat > 0) {
    return benchmark::RegisterBenchmark(name, func, std::forward<ArgsToCallOp>(func_args)...)
        ->ArgNames(arg_names)
        ->Args(args)
        ->UseRealTime()
        ->Iterations(repeat);
  } else {
    return benchmark::RegisterBenchmark(name, func, std::forward<ArgsToCallOp>(func_args)...)
        ->ArgNames(arg_names)
        ->Args(args)
        ->UseRealTime();
  }
}

}  // namespace KokkosKernelsBenchmark

#endif
