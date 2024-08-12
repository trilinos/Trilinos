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
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif
