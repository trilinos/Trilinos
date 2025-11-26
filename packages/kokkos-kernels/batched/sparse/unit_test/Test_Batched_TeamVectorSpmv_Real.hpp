// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_teamvector_spmv_nt_float_float) {
  typedef ::Test::Spmv::ParamTag<Trans::NoTranspose> param_tag_type;
  test_batched_teamvector_spmv<TestDevice, float, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_teamvector_spmv_nt_double_double) {
  typedef ::Test::Spmv::ParamTag<Trans::NoTranspose> param_tag_type;
  test_batched_teamvector_spmv<TestDevice, double, double, param_tag_type>();
}
#endif
