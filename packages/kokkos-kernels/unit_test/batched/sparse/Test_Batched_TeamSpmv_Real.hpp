
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_spmv_nt_float_float) {
  typedef ::Test::Spmv::ParamTag<Trans::NoTranspose> param_tag_type;
  test_batched_team_spmv<TestExecSpace, float, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_spmv_nt_double_double) {
  typedef ::Test::Spmv::ParamTag<Trans::NoTranspose> param_tag_type;
  test_batched_team_spmv<TestExecSpace, double, double, param_tag_type>();
}
#endif
