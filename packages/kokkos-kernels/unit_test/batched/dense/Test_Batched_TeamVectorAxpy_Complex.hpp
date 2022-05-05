
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_teamvector_axpy_nt_dcomplex_dcomplex) {
  test_batched_teamvector_axpy<TestExecSpace, Kokkos::complex<double>,
                               Kokkos::complex<double>>();
}

TEST_F(TestCategory, batched_scalar_teamvector_axpy_nt_dcomplex_double) {
  test_batched_teamvector_axpy<TestExecSpace, Kokkos::complex<double>,
                               double>();
}
#endif
