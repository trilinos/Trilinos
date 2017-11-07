
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_team_trsm_l_l_nt_u_float_float ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_l_nt_n_float_float ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_u_nt_u_float_float ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_u_nt_n_float_float ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_r_u_nt_u_float_float ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_r_u_nt_n_float_float ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_team_trsm_l_l_nt_u_double_double ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_l_nt_n_double_double ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_u_nt_u_double_double ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_l_u_nt_n_double_double ) {
  typedef ::Test::ParamTag<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_r_u_nt_u_double_double ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_team_trsm_r_u_nt_n_double_double ) {
  typedef ::Test::ParamTag<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_trsm<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
#endif

