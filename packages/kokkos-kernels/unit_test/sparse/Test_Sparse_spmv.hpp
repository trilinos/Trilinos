#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>

#include<KokkosSparse_spmv.hpp>
#include<KokkosKernels_TestUtils.hpp>
#include<KokkosKernels_IOUtils.hpp>
#include<KokkosKernels_Utils.hpp>

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test {



template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void sequential_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
    typename y_vector_type::non_const_value_type alpha, typename y_vector_type::non_const_value_type beta){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type size_type_view_t;
  typedef typename graph_t::entries_type lno_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename size_type_view_t::non_const_value_type size_type;
  typedef typename lno_view_t::non_const_value_type lno_t;
  typedef typename scalar_view_t::non_const_value_type scalar_t;


  typename scalar_view_t::HostMirror h_values = Kokkos::create_mirror_view(input_mat.values);
  Kokkos::deep_copy(h_values,input_mat.values);

  typename lno_view_t::HostMirror h_entries = Kokkos::create_mirror_view(input_mat.graph.entries);
  Kokkos::deep_copy(h_entries,input_mat.graph.entries);

  typename size_type_view_t::HostMirror h_rowmap = Kokkos::create_mirror_view(input_mat.graph.row_map);
  Kokkos::deep_copy(h_rowmap,input_mat.graph.row_map);
  Kokkos::fence();



  typename y_vector_type::HostMirror h_y = Kokkos::create_mirror_view(y);
  typename x_vector_type::HostMirror h_x = Kokkos::create_mirror_view(x);

  KokkosKernels::Impl::safe_device_to_host_deep_copy (x.dimension_0(), x, h_x);
  KokkosKernels::Impl::safe_device_to_host_deep_copy (y.dimension_0(), y, h_y);
  Kokkos::fence();


  lno_t nr = input_mat.numRows();

  for (lno_t i = 0; i < nr; ++i){
    scalar_t result = 0;
    for (size_type j = h_rowmap(i); j < h_rowmap(i+1); ++j){
      lno_t col = h_entries(j);
      scalar_t val = h_values(j);
      scalar_t vector_val = h_x(col);
      result += val * vector_val;
    }
    h_y(i) = beta * h_y(i) + alpha * result;
  }
  KokkosKernels::Impl::safe_host_to_device_deep_copy (y.dimension_0(),  h_y, y);
  Kokkos::fence();
}


template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
    typename y_vector_type::non_const_value_type alpha, typename y_vector_type::non_const_value_type beta){
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename scalar_view_t::value_type ScalarA;
  double eps = std::is_same<ScalarA,float>::value?2*1e-3:1e-7;
  size_t nr = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);
  //KokkosKernels::Impl::print_1Dview(expected_y);
  KokkosSparse::spmv("N", alpha, input_mat, x, beta, y);
  //KokkosKernels::Impl::print_1Dview(y);
  EXPECT_NEAR_KK_1DVIEW(y, expected_y, eps);
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type y, y_vector_type expected_y,
    typename y_vector_type::non_const_value_type alpha,
    typename y_vector_type::non_const_value_type beta, int numMV){
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename scalar_view_t::value_type ScalarA;
  double eps = std::is_same<ScalarA,float>::value?2*1e-3:1e-7;

  Kokkos::deep_copy(expected_y, y);

  Kokkos::fence();

  KokkosSparse::spmv("N", alpha, input_mat, x, beta, y);


  for (int i = 0; i < numMV; ++i){
    auto x_i = Kokkos::subview (x, Kokkos::ALL (), i);

    auto y_i = Kokkos::subview (expected_y, Kokkos::ALL (), i);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta);

    auto y_spmv = Kokkos::subview (y, Kokkos::ALL (), i);
    EXPECT_NEAR_KK_1DVIEW(y_spmv, y_i, eps);
  }
}

}
template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv(lno_t numRows,size_type nnz, lno_t bandwidth, lno_t row_size_variance){

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;




  lno_t numCols = numRows;

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);
  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x ("x", nc);
  y_vector_type output_y ("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x,rand_pool,ScalarX(10));
  Kokkos::fill_random(output_y,rand_pool,ScalarY(10));

  Test::check_spmv(input_mat, input_x, output_y, 1.0, 0.0);
  Test::check_spmv(input_mat, input_x, output_y, 0.0, 1.0);
  Test::check_spmv(input_mat, input_x, output_y, 1.0, 1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, class Device>
void test_spmv_mv(lno_t numRows,size_type nnz, lno_t bandwidth, lno_t row_size_variance, int numMV){
  lno_t numCols = numRows;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  //typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeY;

  ViewTypeX b_x("A",numRows,numMV);
  ViewTypeY b_y("B",numCols,numMV);
  ViewTypeY b_y_copy("B",numCols,numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);
  Kokkos::fill_random(b_x,rand_pool,scalar_t(10));
  Kokkos::fill_random(b_y,rand_pool,scalar_t(10));


  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);
  //lno_t nr = input_mat.numRows();
  //lno_t nc = input_mat.numCols();

  Kokkos::deep_copy(b_y_copy, b_y);


  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 0.0, numMV);
  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 0.0, 1.0, numMV);
  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 1.0, numMV);


}



#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_spmv<SCALAR,ORDINAL,OFFSET,DEVICE> (50000, 50000 * 30, 200, 10); \
  test_spmv<SCALAR,ORDINAL,OFFSET,DEVICE> (50000, 50000 * 30, 100, 10); \
  test_spmv<SCALAR,ORDINAL,OFFSET,DEVICE> (10000, 10000 * 20, 100, 5); \
}

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv_mv ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## LAYOUT ## _ ## DEVICE ) { \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (50000, 50000 * 30, 100, 10, 5); \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (50000, 50000 * 30, 200, 10, 1); \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (10000, 10000 * 20, 100, 5, 10); \
}
  //


#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif



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


