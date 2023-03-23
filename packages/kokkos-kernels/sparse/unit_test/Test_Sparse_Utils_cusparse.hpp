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
// Note: Luc Berger-Vergiat 04/14/21
//       Only include this test if compiling
//       the cuda sparse tests and cuSPARSE
//       is enabled.
#if defined(TEST_CUDA_SPARSE_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosSparse_Utils_cusparse.hpp"

void test_cusparse_safe_call() {
  bool caught_exception = false;

  cusparseStatus_t myStatus = CUSPARSE_STATUS_SUCCESS;
  KOKKOS_CUSPARSE_SAFE_CALL(myStatus);

  try {
    myStatus = CUSPARSE_STATUS_INVALID_VALUE;
    KOKKOS_CUSPARSE_SAFE_CALL(myStatus);
  } catch (std::runtime_error& e) {
    caught_exception = true;
  }

  EXPECT_TRUE(caught_exception == true);
}

TEST_F(TestCategory, sparse_cusparse_safe_call) { test_cusparse_safe_call(); }

#endif  // check for CUDA and cuSPARSE
