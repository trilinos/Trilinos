
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_GMRES_float) {
  test_batched_team_GMRES<TestExecSpace, float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_GMRES_double) {
  test_batched_team_GMRES<TestExecSpace, double>();
}
#endif
