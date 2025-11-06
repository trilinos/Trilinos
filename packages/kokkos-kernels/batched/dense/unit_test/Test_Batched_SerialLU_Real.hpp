// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_lu_float) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestDevice, float, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_lu_double) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestDevice, double, algo_tag_type>();
}
#endif
