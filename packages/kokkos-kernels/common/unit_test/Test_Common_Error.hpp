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

#ifndef TEST_COMMON_ERROR_HPP
#define TEST_COMMON_ERROR_HPP

#include "KokkosKernels_Error.hpp"

void test_kokkoskernels_throw() {
  const std::string my_throw_msg = "Testing Kokkos Kernels' throw_runtime_exception.";
  try {
    KokkosKernels::Impl::throw_runtime_exception(my_throw_msg);
  } catch (const std::runtime_error& e) {
  }
}

TEST_F(TestCategory, common_throw) { test_kokkoskernels_throw(); }

#endif  // TEST_COMMON_ERROR_HPP
