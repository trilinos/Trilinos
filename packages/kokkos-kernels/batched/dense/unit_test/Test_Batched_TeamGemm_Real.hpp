// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#if defined(KOKKOS_BHALF_T_IS_FLOAT)
TEST_F(TestCategory, batched_scalar_team_gemm_nt_nt_bhalf_bhalf) {
  using param_tag_type =
      ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;

  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_nt_bhalf_bhalf) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;

  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_nt_t_bhalf_bhalf) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;

  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_t_bhalf_bhalf) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;

  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
#endif  // KOKKOS_BHALF_T_IS_FLOAT

#if defined(KOKKOS_HALF_T_IS_FLOAT)
TEST_F(TestCategory, batched_scalar_team_gemm_nt_nt_half_half) {
  using param_tag_type =
      ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;

  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_nt_half_half) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;

  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_nt_t_half_half) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;

  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_t_half_half) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;

  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                        Algo::Gemm::Unblocked>();
}
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_gemm_nt_nt_float_float) {
  using param_tag_type =
      ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_nt_float_float) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_nt_t_float_float) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_t_float_float) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_gemm_nt_nt_double_double) {
  using param_tag_type =
      ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_nt_double_double) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_nt_t_double_double) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_gemm_t_t_double_double) {
  using param_tag_type = ::Test::TeamGemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_teamgemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif
