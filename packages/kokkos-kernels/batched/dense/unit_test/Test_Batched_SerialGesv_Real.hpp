// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_gesv_static_pivoting_float) {
  test_batched_gesv<TestDevice, float, KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_serial_gesv_no_pivoting_float) {
  test_batched_gesv<TestDevice, float, KokkosBatched::Gesv::NoPivoting>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_gesv_static_pivoting_double) {
  test_batched_gesv<TestDevice, double, KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_serial_gesv_no_pivoting_double) {
  test_batched_gesv<TestDevice, double, KokkosBatched::Gesv::NoPivoting>();
}
#endif
