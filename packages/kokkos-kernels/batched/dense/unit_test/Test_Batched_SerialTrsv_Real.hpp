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
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_u_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_n_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_u_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_n_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_u_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_n_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_u_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_n_float_float) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_u_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_n_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_u_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_n_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_u_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_n_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_u_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_n_double_double) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
#endif
