// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// NO TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trtri_u_n_scomplex_scomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_u_u_scomplex_scomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_n_scomplex_scomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_u_scomplex_scomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>(128);
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_scalar_serial_trtri_u_n_dcomplex_dcomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_u_u_dcomplex_dcomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Upper, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_n_dcomplex_dcomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>(128);
}
TEST_F(TestCategory, batched_scalar_serial_trtri_l_u_dcomplex_dcomplex) {
  typedef ::Test::Trtri::ParamTag<Uplo::Lower, Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;

  test_batched_trtri<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>(128);
}
#endif
