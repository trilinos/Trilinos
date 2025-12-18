// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
// Note: Luc Berger-Vergiat 10/25/21
//       Only include this test if compiling
//       the cuda sparse tests and cuSPARSE
//       is enabled.
#if defined(TEST_HIP_SPARSE_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <rocsparse/rocsparse.h>
#include "KokkosSparse_Utils_rocsparse.hpp"

void test_rocsparse_version() {
  // Print version
  rocsparse_handle handle;
  rocsparse_create_handle(&handle);

  int ver;
  char rev[64];

  rocsparse_get_version(handle, &ver);
  rocsparse_get_git_rev(handle, rev);

  std::cout << "rocSPARSE version: " << ver / 100000 << "." << ver / 100 % 1000 << "." << ver % 100 << "-" << rev
            << std::endl;

  rocsparse_destroy_handle(handle);
}

// Check that the wrapper macro
// detects error status correctly
void test_rocsparse_safe_call() {
  bool caught_exception = false;

  rocsparse_status myStatus = rocsparse_status_success;
  KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(myStatus);

  try {
    myStatus = rocsparse_status_internal_error;
    KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(myStatus);
  } catch (std::runtime_error& e) {
    caught_exception = true;
  }

  EXPECT_TRUE(caught_exception == true);
}

// Check that we can create a handle
// using the singleton class if it
// fails it throws an error with the
// KOKKOS_ROCBLAS_SAFE_CALL_IMPL macro
void test_rocsparse_singleton() {
  KokkosKernels::Impl::RocsparseSingleton& s = KokkosKernels::Impl::RocsparseSingleton::singleton();
  (void)s;
}

TEST_F(TestCategory, sparse_rocsparse_version) { test_rocsparse_version(); }
TEST_F(TestCategory, sparse_rocsparse_safe_call) { test_rocsparse_safe_call(); }
TEST_F(TestCategory, sparse_rocsparse_singleton) { test_rocsparse_singleton(); }

#endif  // check for HIP and rocSPARSE
