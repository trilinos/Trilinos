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

#if defined(KOKKOSKERNELS_INST_FLOAT)
// NO TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_nt_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_nt_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_nt_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_nt_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_nt_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_nt_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_t_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_t_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_t_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_t_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_t_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_t_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
// CONJUGATE TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_ct_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_ct_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_ct_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_ct_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_ct_u_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_ct_n_float_float) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_nt_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_nt_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_nt_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_nt_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_nt_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_nt_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_t_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_t_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_t_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_t_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_t_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_t_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
// CONJUGATE TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_ct_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_l_ct_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_ct_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_l_u_ct_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_ct_u_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trmm_r_u_ct_n_double_double) {
  typedef ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;

  test_batched_trmm<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
#endif
