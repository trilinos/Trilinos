
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_teamvector_solve_utv_float ) {
  typedef Algo::UTV::Unblocked algo_tag_type;
  test_batched_solve_utv<TestExecSpace,float,int,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_teamvector_solve_utv_double ) {
  typedef Algo::UTV::Unblocked algo_tag_type;
  test_batched_solve_utv<TestExecSpace,double,int,algo_tag_type>();
}
#endif

