
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_GMRES_float) {
  test_batched_serial_GMRES<TestExecSpace, float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_GMRES_double) {
  test_batched_serial_GMRES<TestExecSpace, double>();
}
#endif
