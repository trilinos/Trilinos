#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)

/// dcomplex, dcomplex

TEST_F( TestCategory, batched_scalar_serial_gemv_nt_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_gemv_t_dcomplex_dcomplex ) {
  typedef ::Test::ParamTag<Trans::Transpose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>();
}
// TEST_F( TestCategory, batched_scalar_serial_gemv_ct_dcomplex_dcomplex ) {
//   typedef ::Test::ParamTag<Trans::ConjTranspose> param_tag_type;
//   typedef Algo::Gemv::Blocked algo_tag_type;
//   test_batched_gemv<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>();
// }

/// dcomplex, double

TEST_F( TestCategory, batched_scalar_serial_gemv_nt_dcomplex_double ) {
  typedef ::Test::ParamTag<Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,Kokkos::complex<double>,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_gemv_t_dcomplex_double ) {
  typedef ::Test::ParamTag<Trans::Transpose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,Kokkos::complex<double>,double,param_tag_type,algo_tag_type>();
}
// TEST_F( TestCategory, batched_scalar_serial_gemv_ct_dcomplex_double ) {
//   typedef ::Test::ParamTag<Trans::ConjTranspose> param_tag_type;
//   typedef Algo::Gemv::Blocked algo_tag_type;
//   test_batched_gemv<TestExecSpace,Kokkos::complex<double>,double,param_tag_type,algo_tag_type>();
// }

#endif
