//
// Created by Poliakoff, David Zoeller on 4/26/21.
//
#ifndef KOKKOSKERNELS_BLAS_TRACKED_TESTING_HPP
#define KOKKOSKERNELS_BLAS_TRACKED_TESTING_HPP

#include <common/RAJAPerfSuite.hpp>
#include <common/Executor.hpp>

#include "KokkosBlas_dot_perf_test.hpp"
#include "KokkosBlas_team_dot_perf_test.hpp"

namespace test {
namespace blas {

// Register kernels per test

void build_dot_executor(rajaperf::Executor& exec, int, char*[],
                        const rajaperf::RunParams& params) {
  for (auto* kernel : construct_dot_kernel_base(params)) {
    exec.registerKernel("BLAS", kernel);
  }
}
// Team Dot build_executor

void build_team_dot_executor(rajaperf::Executor& exec, int, char*[],
                             const rajaperf::RunParams& params) {
  for (auto* kernel : construct_team_dot_kernel_base(params)) {
    exec.registerKernel("BLAS", kernel);
  }
}

void build_blas_executor(rajaperf::Executor& exec, int argc, char* argv[],
                         const rajaperf::RunParams& params) {
  exec.registerGroup("BLAS");
  build_dot_executor(exec, argc, argv, params);
  build_team_dot_executor(exec, argc, argv, params);
}

}  // namespace blas
}  // namespace test

#endif  // KOKKOSKERNELS_TRACKED_TESTING_HPP
