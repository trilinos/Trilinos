#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>

#include<KokkosSparse_trsv.hpp>
#include<KokkosSparse_spmv.hpp>
#include<KokkosKernels_TestUtils.hpp>
#include<KokkosKernels_IOUtils.hpp>

#include<KokkosKernels_Utils.hpp>

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

namespace Test {
//TODO: remove this once MD develop branch is merge.
//The below functionolity exists in SparseUtils.


template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_trsv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type b, y_vector_type expected_x, int numMV, const char uplo[]){
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename scalar_view_t::value_type ScalarA;
  double eps = (std::is_same<ScalarA,float>::value ? 2*1e-2
               : (std::is_same<ScalarA,std::complex<float>>::value || std::is_same<ScalarA,Kokkos::complex<float>>::value )? 2*1e-1 : 1e-7 );

  Kokkos::fence();
  KokkosSparse::trsv(uplo, "N", "N", input_mat, b, x);

  for (int i = 0; i < numMV; ++i){
    auto x_i = Kokkos::subview (x, Kokkos::ALL (), i);

    auto expected_x_i = Kokkos::subview (expected_x, Kokkos::ALL (), i);

    EXPECT_NEAR_KK_1DVIEW(expected_x_i, x_i, eps);
  }
}
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, class Device>
void test_trsv_mv(lno_t numRows,size_type nnz, lno_t bandwidth, lno_t row_size_variance, int numMV){
  lno_t numCols = numRows;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  //typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeY;

  ViewTypeX b_x("A",numRows,numMV);
  ViewTypeY b_y("B",numCols,numMV);
  ViewTypeX b_x_copy("B",numCols,numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);
  Kokkos::fill_random(b_x_copy,rand_pool,scalar_t(10));

  typename ViewTypeY::non_const_value_type alpha = 1;
  typename ViewTypeY::non_const_value_type beta = 0;


  //this function creates a dense lower and upper triangular matrix.
  //TODO: SHOULD CHANGE IT TO SPARSE
  crsMat_t lower_part = KokkosKernels::Impl::kk_generate_triangular_sparse_matrix<crsMat_t>('L', numRows,numCols,nnz,row_size_variance, bandwidth);
  KokkosSparse::spmv("N", alpha, lower_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(lower_part, b_x, b_y, b_x_copy, numMV, "L");
  //typedef typename Kokkos::View<lno_t*, layout, Device> indexview;

  crsMat_t upper_part = KokkosKernels::Impl::kk_generate_triangular_sparse_matrix<crsMat_t>('U', numRows,numCols,nnz,row_size_variance, bandwidth);
  KokkosSparse::spmv("N", alpha, upper_part, b_x_copy, beta, b_y);
  Test::check_trsv_mv(upper_part, b_x, b_y, b_x_copy, numMV, "U");

}




#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## trsv_mv ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## LAYOUT ## _ ## DEVICE ) { \
  test_trsv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (5000, 5000 * 30, 200, 10, 1); \
  test_trsv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (5000, 5000 * 30, 100, 10, 5); \
  test_trsv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (1000, 1000 * 20, 100, 5, 10); \
}


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif








#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, size_t, LayoutRight, TestExecSpace)
#endif




#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, size_t, LayoutRight, TestExecSpace)
#endif


