
#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// NO TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_nt_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_nt_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_nt_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_nt_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_nt_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_nt_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
// TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_t_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_t_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_t_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_t_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_t_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_t_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
// CONJUGATE TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_ct_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_ct_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_ct_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_ct_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_ct_u_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_ct_n_scomplex_scomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,algo_tag_type>(128);
}
#endif


#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// NO TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_nt_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_nt_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_nt_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_nt_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_nt_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_nt_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
// TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_t_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_t_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_t_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_t_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_t_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::Transpose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_t_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::Transpose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
// CONJUGATE TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_ct_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_l_ct_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_ct_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_l_u_ct_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_ct_u_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::ConjTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
TEST_F( TestCategory, batched_scalar_serial_trmm_r_u_ct_n_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trmm::Unblocked algo_tag_type;
  
  test_batched_trmm<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>(128);
}
#endif

