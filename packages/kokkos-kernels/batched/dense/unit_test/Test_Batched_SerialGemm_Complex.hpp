// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)

/// fcomplex, fcomplex

TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_nt_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_t_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_c_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_c_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_c_fcomplex_fcomplex) {
  using param_tag_type =
      ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}

/// fcomplex, float
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_nt_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_t_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_c_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_c_fcomplex_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_c_fcomplex_float) {
  using param_tag_type =
      ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<float>, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}

#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)

/// dcomplex, dcomplex

TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_nt_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_t_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_c_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_c_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_c_dcomplex_dcomplex) {
  using param_tag_type =
      ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}

/// dcomplex, double
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_nt_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_t_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_c_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_c_dcomplex_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_c_c_dcomplex_double) {
  using param_tag_type =
      ::Test::Gemm::ParamTag<KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::ConjTranspose>;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}

#endif
