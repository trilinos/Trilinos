
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_axpy_nt_dcomplex_dcomplex) {
  test_batched_team_axpy<TestExecSpace, Kokkos::complex<double>,
                         Kokkos::complex<double>>();
}

TEST_F(TestCategory, batched_scalar_team_axpy_nt_dcomplex_double) {
  test_batched_team_axpy<TestExecSpace, Kokkos::complex<double>, double>();
}
#endif
