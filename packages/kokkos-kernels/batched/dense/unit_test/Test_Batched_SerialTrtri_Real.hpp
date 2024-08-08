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
TEST_F(TestCategory, batched_scalar_serial_trtri_u_n_float_float) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_u_u_float_float) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_n_float_float) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_u_float_float) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, float, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trtri_u_n_double_double) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_u_u_double_double) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_n_double_double) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_u_double_double) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, double, double, param_tag_type, algo_tag_type>();
}
#endif
