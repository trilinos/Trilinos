// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
// NO TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_float) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_float) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_double) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_double) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif
