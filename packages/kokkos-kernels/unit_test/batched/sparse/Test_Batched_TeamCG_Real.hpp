
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_CG_float) {
  test_batched_team_CG<TestExecSpace, float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_CG_double) {
  test_batched_team_CG<TestExecSpace, double>();
}
#endif
