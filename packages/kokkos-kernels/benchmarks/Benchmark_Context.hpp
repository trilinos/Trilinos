// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_PERFTEST_BENCHMARK_CONTEXT_HPP
#define KOKKOSKERNELS_PERFTEST_BENCHMARK_CONTEXT_HPP

#include "KokkosKernels_PrintConfiguration.hpp"

#include <cstdlib>
#include <string>

#include <benchmark/benchmark.h>

#include <Kokkos_Core.hpp>
#include <KokkosKernels_PrintConfiguration.hpp>
#include <KokkosKernels_Version_Info.hpp>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime.h>
#endif

#if defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

#if defined(KOKKOS_ENABLE_SYCL)
#include <sycl/sycl.hpp>
#endif

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

/// \brief Add GPU architecture information to benchmark context.
/// Queries device properties via CUDA/HIP runtime APIs or SYCL device info.
inline void add_gpu_context() {
#if defined(KOKKOS_ENABLE_CUDA)
  int num_devices = 0;
  if (cudaGetDeviceCount(&num_devices) != cudaSuccess) return;
  for (int dev = 0; dev < num_devices; ++dev) {
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, dev) != cudaSuccess) continue;

    std::string prefix = "GPU" + std::to_string(dev) + "_";

    benchmark::AddCustomContext(prefix + "name", prop.name);
    benchmark::AddCustomContext(prefix + "compute_capability",
                                std::to_string(prop.major) + "." + std::to_string(prop.minor));
    benchmark::AddCustomContext(prefix + "num_multiprocessors", std::to_string(prop.multiProcessorCount));
    benchmark::AddCustomContext(prefix + "total_global_mem_GiB", std::to_string(prop.totalGlobalMem / (1 << 30)));
    benchmark::AddCustomContext(prefix + "l2_cache_size_MiB", std::to_string(prop.l2CacheSize / (1 << 20)));
    benchmark::AddCustomContext(prefix + "clock_rate_MHz", std::to_string(prop.clockRate / 1000));
    benchmark::AddCustomContext(prefix + "memory_clock_rate_MHz", std::to_string(prop.memoryClockRate / 1000));
    benchmark::AddCustomContext(prefix + "memory_bus_width_bits", std::to_string(prop.memoryBusWidth));
    // Peak memory bandwidth in GB/s: 2 * memClk(Hz) * busWidth(bits) / 8
    double peak_bw = 2.0 * (prop.memoryClockRate * 1000.0) * (prop.memoryBusWidth / 8.0) / 1e9;
    benchmark::AddCustomContext(prefix + "peak_memory_bandwidth_GBs", std::to_string(static_cast<int>(peak_bw + 0.5)));
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  int num_devices = 0;
  if (hipGetDeviceCount(&num_devices) != hipSuccess) return;
  for (int dev = 0; dev < num_devices; ++dev) {
    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, dev) != hipSuccess) continue;

    std::string prefix = "GPU" + std::to_string(dev) + "_";

    benchmark::AddCustomContext(prefix + "name", prop.name);
    benchmark::AddCustomContext(prefix + "compute_capability",
                                std::to_string(prop.major) + "." + std::to_string(prop.minor));
    benchmark::AddCustomContext(prefix + "num_multiprocessors", std::to_string(prop.multiProcessorCount));
    benchmark::AddCustomContext(prefix + "total_global_mem_GiB", std::to_string(prop.totalGlobalMem / (1 << 30)));
    benchmark::AddCustomContext(prefix + "l2_cache_size_MiB", std::to_string(prop.l2CacheSize / (1 << 20)));
    benchmark::AddCustomContext(prefix + "clock_rate_MHz", std::to_string(prop.clockRate / 1000));
    benchmark::AddCustomContext(prefix + "memory_clock_rate_MHz", std::to_string(prop.memoryClockRate / 1000));
    benchmark::AddCustomContext(prefix + "memory_bus_width_bits", std::to_string(prop.memoryBusWidth));
    double peak_bw = 2.0 * (prop.memoryClockRate * 1000.0) * (prop.memoryBusWidth / 8.0) / 1e9;
    benchmark::AddCustomContext(prefix + "peak_memory_bandwidth_GBs", std::to_string(static_cast<int>(peak_bw + 0.5)));
  }
#endif

#if defined(KOKKOS_ENABLE_SYCL)
  int dev = 0;
  for (const auto& platform : sycl::platform::get_platforms()) {
    for (const auto& device : platform.get_devices()) {
      if (device.is_gpu()) {
        std::string prefix = "GPU" + std::to_string(dev) + "_";
        benchmark::AddCustomContext(prefix + "name", device.get_info<sycl::info::device::name>());
        benchmark::AddCustomContext(prefix + "vendor", device.get_info<sycl::info::device::vendor>());
        benchmark::AddCustomContext(prefix + "num_compute_units",
                                    std::to_string(device.get_info<sycl::info::device::max_compute_units>()));
        benchmark::AddCustomContext(
            prefix + "total_global_mem_GiB",
            std::to_string(device.get_info<sycl::info::device::global_mem_size>() / (1ULL << 30)));
        benchmark::AddCustomContext(prefix + "max_clock_frequency_MHz",
                                    std::to_string(device.get_info<sycl::info::device::max_clock_frequency>()));
        benchmark::AddCustomContext(prefix + "local_mem_size_KiB",
                                    std::to_string(device.get_info<sycl::info::device::local_mem_size>() / 1024));
        ++dev;
      }
    }
  }
#endif
}

/// \brief Gather all context information and add it to benchmark context
inline void add_benchmark_context(bool verbose = false) {
  add_kokkos_configuration(verbose);
  add_version_info();
  add_env_info();
  add_gpu_context();
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
