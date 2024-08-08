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

void build_dot_executor(rajaperf::Executor& exec, int, char*[], const rajaperf::RunParams& params) {
  for (auto* kernel : construct_dot_kernel_base(params)) {
    exec.registerKernel("BLAS", kernel);
  }
}
// Team Dot build_executor

void build_team_dot_executor(rajaperf::Executor& exec, int, char*[], const rajaperf::RunParams& params) {
  for (auto* kernel : construct_team_dot_kernel_base(params)) {
    exec.registerKernel("BLAS", kernel);
  }
}

void build_blas_executor(rajaperf::Executor& exec, int argc, char* argv[], const rajaperf::RunParams& params) {
  exec.registerGroup("BLAS");
  build_dot_executor(exec, argc, argv, params);
  build_team_dot_executor(exec, argc, argv, params);
}

}  // namespace blas
}  // namespace test

#endif  // KOKKOSKERNELS_TRACKED_TESTING_HPP
