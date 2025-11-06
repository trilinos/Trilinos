// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
