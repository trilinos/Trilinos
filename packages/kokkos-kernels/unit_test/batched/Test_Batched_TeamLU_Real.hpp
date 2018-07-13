
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_team_lu_float ) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestExecSpace,float,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_team_lu_double ) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestExecSpace,double,algo_tag_type>();
}
#endif

