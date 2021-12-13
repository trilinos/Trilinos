
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_scalar_serial_solvelu_dcomplex ) {
  //printf("Batched serial solveLU - double complex - algorithm type: Unblocked\n");
  test_batched_solvelu<TestExecSpace,Kokkos::complex<double>,Algo::SolveLU::Unblocked>();
  //printf("Batched serial solveLU - double complex - algorithm type: Blocked\n");
  test_batched_solvelu<TestExecSpace,Kokkos::complex<double>,Algo::SolveLU::Blocked>();
}
#endif
