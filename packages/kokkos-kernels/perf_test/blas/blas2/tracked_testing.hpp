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
#ifndef KOKKOSKERNELS_BLAS2_TRACKED_TESTING_HPP
#define KOKKOSKERNELS_BLAS2_TRACKED_TESTING_HPP

#include <common/RAJAPerfSuite.hpp>
#include <common/Executor.hpp>

#include "KokkosBlas2_gemv_perf_test.hpp"

namespace test {
namespace blas2 {

// Register kernels for a specific test
void build_gemv_executor(rajaperf::Executor& exec, int argc, char* argv[], const rajaperf::RunParams& params) {
  for (auto* kernel : construct_gemv_kernel_base(params)) {
    exec.registerKernel("BLAS2", kernel);
  }
}

void build_blas2_executor(rajaperf::Executor& exec, int argc, char* argv[], const rajaperf::RunParams& params) {
  exec.registerGroup("BLAS2");
  build_gemv_executor(exec, argc, argv, params);
}

}  // namespace blas2
}  // namespace test

#endif  // KOKKOSKERNELS_TRACKED_TESTING_HPP
