#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_gesv_static_pivoting_float) {
  test_batched_team_gesv<TestExecSpace, float,
                         KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_team_gesv_no_pivoting_float) {
  test_batched_team_gesv<TestExecSpace, float,
                         KokkosBatched::Gesv::NoPivoting>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_gesv_static_pivoting_double) {
  test_batched_team_gesv<TestExecSpace, double,
                         KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_team_gesv_no_pivoting_double) {
  test_batched_team_gesv<TestExecSpace, double,
                         KokkosBatched::Gesv::NoPivoting>();
}
#endif
