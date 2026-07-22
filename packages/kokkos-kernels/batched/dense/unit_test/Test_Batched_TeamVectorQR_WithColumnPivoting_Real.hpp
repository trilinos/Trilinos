// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_teamvector_qr_with_columnpivoting_float) {
  typedef Algo::QR::Unblocked algo_tag_type;
  test_batched_qr_with_columnpivoting<TestDevice, float, int, algo_tag_type>();
}
#endif

// FIXME_SYCL timeout
#ifndef KOKKOS_ENABLE_SYCL
#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_teamvector_qr_with_columnpivoting_double) {
  typedef Algo::QR::Unblocked algo_tag_type;
  test_batched_qr_with_columnpivoting<TestDevice, double, int, algo_tag_type>();
}
#endif
#endif
