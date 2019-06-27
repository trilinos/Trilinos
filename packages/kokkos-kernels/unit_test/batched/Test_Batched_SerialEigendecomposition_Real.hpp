#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_serial_eigendecomposition_float ) {
  test_batched_eigendecomposition<TestExecSpace,float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_serial_eigendecomposition_double ) {
  test_batched_eigendecomposition<TestExecSpace,double>();
}
#endif


