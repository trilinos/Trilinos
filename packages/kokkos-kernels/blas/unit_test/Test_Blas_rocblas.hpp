// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
// Note: Luc Berger-Vergiat 10/25/21
//       Only include this test if compiling
//       the cuda sparse tests and cuSPARSE
//       is enabled.
#if defined(TEST_HIP_BLAS_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "KokkosBlas_tpl_spec.hpp"
#include <rocblas/rocblas.h>

// Just check if we can build against
// rocblas and get the library version
void test_rocblas_version() {
  // Print version
  size_t size = 128;
  // rocblas_get_version_string_size(&size);
  std::string version(size - 1, '\0');
  rocblas_get_version_string(const_cast<char*>(version.data()), size);
  std::cout << "rocBLAS version: " << version << "\n" << std::endl;
}

// Check that the wrapper macro
// detects error status correctly
void test_rocblas_safe_call() {
  bool caught_exception = false;

  rocblas_status myStatus = rocblas_status_success;
  KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(myStatus);

  try {
    myStatus = rocblas_status_internal_error;
    KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(myStatus);
  } catch (std::runtime_error& e) {
    caught_exception = true;
  }

  EXPECT_TRUE(caught_exception == true);
}

// Check that we can create a handle
// using the singleton class if it
// fails it throws an error with the
// KOKKOS_ROCBLAS_SAFE_CALL_IMPL macro
void test_rocblas_singleton() {
  KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();
  (void)s;
}

TEST_F(TestCategory, blas_rocblas_version) { test_rocblas_version(); }
TEST_F(TestCategory, blas_rocblas_safe_call) { test_rocblas_safe_call(); }
TEST_F(TestCategory, blas_rocblas_singleton) { test_rocblas_singleton(); }

#endif  // check for HIP and rocBLAS
