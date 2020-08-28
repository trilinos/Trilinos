
#if defined(KOKKOSKERNELS_INST_FLOAT)
// NO TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trtri_u_n_float_float ) {
  typedef ::Test::ParamTag<Uplo::Upper,Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_u_u_float_float ) {
  typedef ::Test::ParamTag<Uplo::Upper,Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_l_n_float_float ) {
  typedef ::Test::ParamTag<Uplo::Lower,Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_l_u_float_float ) {
  typedef ::Test::ParamTag<Uplo::Lower,Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,float,float,param_tag_type,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
// NO TRANSPOSE
TEST_F( TestCategory, batched_scalar_serial_trtri_u_n_double_double ) {
  typedef ::Test::ParamTag<Uplo::Upper,Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_u_u_double_double ) {
  typedef ::Test::ParamTag<Uplo::Upper,Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_l_n_double_double ) {
  typedef ::Test::ParamTag<Uplo::Lower,Diag::NonUnit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
TEST_F( TestCategory, batched_scalar_serial_trtri_l_u_double_double ) {
  typedef ::Test::ParamTag<Uplo::Lower,Diag::Unit> param_tag_type;
  typedef Algo::Trtri::Unblocked algo_tag_type;
  
  test_batched_trtri<TestExecSpace,double,double,param_tag_type,algo_tag_type>();
}
#endif

