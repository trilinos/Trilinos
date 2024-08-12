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
#ifndef KOKKOSKERNELS_TRACKED_TESTING_HPP
#define KOKKOSKERNELS_TRACKED_TESTING_HPP
#include <common/RAJAPerfSuite.hpp>
#include <common/Executor.hpp>
#include "KokkosSparse_spmv_test.hpp"

namespace test {
namespace sparse {
void build_executor(rajaperf::Executor& exec, int, char*[], const rajaperf::RunParams& params) {
  exec.registerGroup("Sparse");
  for (auto* kernel : make_spmv_kernel_base(params)) {
    exec.registerKernel("Sparse", kernel);
  }
}
}  // namespace sparse
}  // namespace test
#endif  // KOKKOSKERNELS_TRACKED_TESTING_HPP
