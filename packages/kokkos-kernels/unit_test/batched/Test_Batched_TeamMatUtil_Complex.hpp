
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_scalar_team_set_dcomplex_dcomplex ) {
  test_batched_matutil<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,::Test::BatchedSet>();
}
TEST_F( TestCategory, batched_scalar_team_scale_dcomplex_dcomplex ) {
  test_batched_matutil<TestExecSpace,Kokkos::complex<double>,Kokkos::complex<double>,::Test::BatchedScale>();
}
TEST_F( TestCategory, batched_scalar_team_set_dcomplex_double ) {
  test_batched_matutil<TestExecSpace,Kokkos::complex<double>,double,::Test::BatchedSet>();
}
TEST_F( TestCategory, batched_scalar_team_scale_dcomplex_double ) {
  test_batched_matutil<TestExecSpace,Kokkos::complex<double>,double,::Test::BatchedScale>();
}
#endif
