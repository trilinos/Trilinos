#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F( TestCategory, batched_scalar_team_vector_gemm_nt_nt_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Trans::NoTranspose,Trans::NoTranspose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_t_nt_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Trans::Transpose,Trans::NoTranspose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_nt_t_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Trans::NoTranspose,Trans::Transpose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_t_t_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Trans::Transpose,Trans::Transpose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Unblocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_scalar_team_vector_gemm_nt_nt_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::NoTranspose,Trans::NoTranspose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_t_nt_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::Transpose,Trans::NoTranspose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_nt_t_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::NoTranspose,Trans::Transpose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Unblocked>();
}
TEST_F( TestCategory, batched_scalar_team_vector_gemm_t_t_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::Transpose,Trans::Transpose> param_tag_type;

  //test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Unblocked>();
}
#endif
