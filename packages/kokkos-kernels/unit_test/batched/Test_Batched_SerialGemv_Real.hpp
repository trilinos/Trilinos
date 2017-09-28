
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_serial_gemv_nt_float_float ) {
  typedef ::Test::ParamTag<Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_gemv_t_float_float ) {
  typedef ::Test::ParamTag<Trans::Transpose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_serial_gemv_nt_double_double ) {
  typedef ::Test::ParamTag<Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_gemv_t_double_double ) {
  typedef ::Test::ParamTag<Trans::Transpose> param_tag_type;
  typedef Algo::Gemv::Blocked algo_tag_type;
  test_batched_gemv<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
#endif
