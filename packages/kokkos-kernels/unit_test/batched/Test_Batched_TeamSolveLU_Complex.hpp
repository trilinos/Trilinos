
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_scalar_team_solvelu_dcomplex ) {
  //printf("Batched team solveLU - double complex - algorithm type: Unblocked\n");
  test_batched_solvelu<TestExecSpace,Kokkos::complex<double>,Algo::SolveLU::Unblocked>();
  //printf("Batched team solveLU - double complex - algorithm type: Blocked\n");
  test_batched_solvelu<TestExecSpace,Kokkos::complex<double>,Algo::SolveLU::Blocked>();
}
#endif
