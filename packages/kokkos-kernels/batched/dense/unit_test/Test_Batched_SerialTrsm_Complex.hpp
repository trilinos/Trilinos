// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// LeftSide
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}

// RightSide
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// LeftSide
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_l_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_l_u_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}

// RightSide
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_l_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsm_r_u_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trmm::ParamTag<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsm::Unblocked;
  test_batched_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
