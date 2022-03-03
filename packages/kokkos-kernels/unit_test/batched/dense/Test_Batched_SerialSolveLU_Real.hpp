
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_serial_solvelu_float ) {
  //printf("Batched serial solveLU - float - algorithm type: Unblocked\n");
  test_batched_solvelu<TestExecSpace,float,Algo::SolveLU::Unblocked>();
  //printf("Batched serial solveLU - float - algorithm type: Blocked\n");
  test_batched_solvelu<TestExecSpace,float,Algo::SolveLU::Blocked>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_serial_solvelu_double ) {
  //printf("Batched serial solveLU - double - algorithm type: Unblocked\n");
  test_batched_solvelu<TestExecSpace,double,Algo::SolveLU::Unblocked>();
  //printf("Batched serial solveLU - double - algorithm type: Blocked\n");
  test_batched_solvelu<TestExecSpace,double,Algo::SolveLU::Blocked>();
}
#endif

