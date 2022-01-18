  
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_serial_set_float_float ) {
  test_batched_matutil<TestExecSpace,float,float,::Test::BatchedSet>();
}
TEST_F( TestCategory, batched_scalar_serial_scale_float_float ) {
  test_batched_matutil<TestExecSpace,float,float,::Test::BatchedScale>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_serial_set_double_double ) {
  test_batched_matutil<TestExecSpace,double,double,::Test::BatchedSet>();
}
TEST_F( TestCategory, batched_scalar_serial_scale_double_double ) {
  test_batched_matutil<TestExecSpace,double,double,::Test::BatchedScale>();
}
#endif

