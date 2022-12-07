#if defined(KOKKOS_BHALF_T_IS_FLOAT)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::bhalfScalarType,
                    ::Test::bhalfScalarType, param_tag_type>();
}
#endif  // KOKKOS_BHALF_T_IS_FLOAT

#if defined(KOKKOS_HALF_T_IS_FLOAT)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
  test_batched_gemm<TestExecSpace, ::Test::halfScalarType,
                    ::Test::halfScalarType, param_tag_type>();
}
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if defined(KOKKOSKERNELS_INST_FLOAT)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
  test_batched_gemm<TestExecSpace, float, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Left>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose,
                                 BatchLayout::Right>
      param_tag_type;

  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
  test_batched_gemm<TestExecSpace, double, double, param_tag_type>();
}
#endif
