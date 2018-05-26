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

template < class VectorType0, class VectorType1, class AT_Type >
struct fSPMV {
  typedef int value_type;
  typedef Kokkos::Details::ArithTraits<typename AT_Type::non_const_value_type> AT;

  VectorType0 expected_y;
  VectorType1 y;
  double eps;

  fSPMV(const VectorType0 & _ex_y, const VectorType1 & _y, const double _eps)
  : expected_y(_ex_y)
  , y(_y)
  , eps(_eps)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, value_type& err ) const {
    if(AT::abs(expected_y(i)-y(i))>eps) err++;
  }
};


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

  KokkosKernels::Impl::safe_device_to_host_deep_copy (x.extent(0), x, h_x);
  KokkosKernels::Impl::safe_device_to_host_deep_copy (y.extent(0), y, h_y);
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
  KokkosKernels::Impl::safe_host_to_device_deep_copy (y.extent(0),  h_y, y);
  Kokkos::fence();
}


template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
    typename y_vector_type::non_const_value_type alpha, typename y_vector_type::non_const_value_type beta){
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::execution_space ExecSpace;
  typedef Kokkos::RangePolicy<ExecSpace> my_exec_space;

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
  typedef Kokkos::Details::ArithTraits<typename y_vector_type::non_const_value_type> AT;
  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv"
                         ,my_exec_space(0, y.extent(0))
                         ,fSPMV<y_vector_type, y_vector_type, y_vector_type>(expected_y,y,eps)
                         ,num_errors);
  if(num_errors>0) printf("KokkosKernels::UnitTests::spmv: %i errors of %i with params: %lf %lf\n",
      num_errors,y.extent_int(0),AT::abs(alpha),AT::abs(beta));
  EXPECT_TRUE(num_errors==0);
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type y, y_vector_type expected_y,
    typename y_vector_type::non_const_value_type alpha,
    typename y_vector_type::non_const_value_type beta, int numMV){
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::execution_space ExecSpace;
  typedef Kokkos::RangePolicy<ExecSpace> my_exec_space;

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
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv_mv"
                           ,my_exec_space(0,y_i.extent(0))
                           ,fSPMV<decltype(y_i), decltype(y_spmv), y_vector_type>(y_i, y_spmv, eps)
                           ,num_errors);
    if(num_errors>0) printf("KokkosKernels::UnitTests::spmv_mv: %i errors of %i for mv %i\n",
        num_errors,y_i.extent_int(0),i);
    EXPECT_TRUE(num_errors==0);
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

//call it if ordinal int and, scalar float and double are instantiated.
template<class DeviceType>
void test_github_issue_101 ()
{
  typedef KokkosSparse::CrsMatrix<float, int, DeviceType> float_matrix_type;
  typedef KokkosSparse::CrsMatrix<double, int, DeviceType> double_matrix_type;
  static_assert (std::is_same<typename float_matrix_type::StaticCrsGraphType,
    typename double_matrix_type::StaticCrsGraphType>::value,
    "Two KokkosSparse::CrsMatrix types that differ only in the type of "
    "matrix values, appear to have two different StaticCrsGraphType "
    "typedefs.  This should never happen.");
  typedef typename float_matrix_type::StaticCrsGraphType graph_type;

  constexpr int numRows = 1;
  constexpr int numCols = 2;
  constexpr double alpha_d = 1.0;
  constexpr double beta_d = 0.0;
  const float EPS_f = std::numeric_limits<float>::epsilon ();

  graph_type G;
  {
    typename graph_type::entries_type colInds ("colInds", numCols);
    auto colInds_h = Kokkos::create_mirror_view (colInds);
    colInds_h[0] = 0;
    colInds_h[1] = 1;
    Kokkos::deep_copy (colInds, colInds_h);

    typedef typename graph_type::row_map_type::non_const_type row_offsets_type;
    row_offsets_type rowOffsets ("rowOffsets", numRows+1);
    auto rowOffsets_h = Kokkos::create_mirror_view (rowOffsets);
    rowOffsets_h[0] = 0; // Entries start at offset 0
    rowOffsets_h[1] = 2; // 2 entries total in the "sparse" matrix
    Kokkos::deep_copy (rowOffsets, rowOffsets_h);

    G = graph_type (colInds, rowOffsets);
  }

  Kokkos::View<double*, DeviceType> x ("x", numCols);
  Kokkos::deep_copy (x, static_cast<double> (1.0));
  Kokkos::View<double*, DeviceType> y ("y", numRows);
  auto y_h = Kokkos::create_mirror_view (y); // we'll want this later

  // Pick some number large enough to exercise all unrolling cases.
  // Sparse mat-vec does or at least used to unroll for 1, 2, ..., 17
  // vectors.  Include a little extra in case the implementers decide
  // to strip-mine that.
  constexpr int numVecs = 22;
  Kokkos::View<double**, Kokkos::LayoutLeft, DeviceType> X ("X", numCols, numVecs);
  Kokkos::deep_copy (X, static_cast<double> (1.0));
  Kokkos::View<double**, Kokkos::LayoutLeft, DeviceType> Y ("Y", numRows, numVecs);
  auto Y_h = Kokkos::create_mirror_view (Y); // we'll want this later

  // Start with the easy test case, where the matrix and the vectors
  // are all double.
  {
    constexpr double ZERO_d = static_cast<double> (0.0);
    constexpr double ONE_d = static_cast<double> (1.0);
    constexpr double TWO_d = static_cast<double> (2.0);

    double_matrix_type A_d ("A_d", G);
    auto A_d_val_h = Kokkos::create_mirror_view (A_d.values);
    A_d_val_h[0] = ONE_d;
    // This cast is deliberate; we want to use float eps here, but as
    // a double-precision number.  This is just a sanity check for
    // accuracy of the sparse mat-vec when not using mixed precision.
    A_d_val_h[1] = static_cast<double> (EPS_f) / TWO_d;
    EXPECT_NE( A_d_val_h[1], ZERO_d ); // just making sure
    Kokkos::deep_copy (A_d.values, A_d_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy (y, ZERO_d);
    KokkosSparse::spmv ("N", alpha_d, A_d, x, beta_d, y);

    Kokkos::deep_copy (y_h, y);
    const double expectedResult_allDouble = static_cast<double> (1.0) +
      static_cast<double> (EPS_f) / static_cast<double> (2.0);
    EXPECT_NE( expectedResult_allDouble, ZERO_d );
    EXPECT_EQ( y_h[0], expectedResult_allDouble );

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng (0, curNumVecs);
      auto X_sub = Kokkos::subview (X, Kokkos::ALL (), vecRng);
      auto Y_sub = Kokkos::subview (Y, Kokkos::ALL (), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy (Y, ZERO_d);
      KokkosSparse::spmv ("N", alpha_d, A_d, X, beta_d, Y);

      Kokkos::deep_copy (Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0,j);
        EXPECT_EQ( actualResult, expectedResult_allDouble );
      }
    }
  }

  // Now exercise the case where the matrix is in float, but the
  // vectors are in double.
  {
    constexpr float ZERO_f = static_cast<float> (0.0);
    constexpr float ONE_f = static_cast<float> (1.0);
    constexpr float TWO_f = static_cast<float> (2.0);
    constexpr double ZERO_d = static_cast<double> (0.0);

    float_matrix_type A_f ("A_f", G);
    auto A_f_val_h = Kokkos::create_mirror_view (A_f.values);
    A_f_val_h[0] = ONE_f;
    A_f_val_h[1] = EPS_f / TWO_f;
    EXPECT_NE( A_f_val_h[1], ZERO_f ); // just making sure
    Kokkos::deep_copy (A_f.values, A_f_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy (y, ZERO_d);
    KokkosSparse::spmv ("N", alpha_d, A_f, x, beta_d, y);

    Kokkos::deep_copy (y_h, y);
    const double expectedResult_mixed = static_cast<double> (1.0) +
      static_cast<double> (EPS_f) / static_cast<double> (2.0);
    EXPECT_NE( expectedResult_mixed, ZERO_d );
    EXPECT_EQ( y_h[0], expectedResult_mixed );

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng (0, curNumVecs);
      auto X_sub = Kokkos::subview (X, Kokkos::ALL (), vecRng);
      auto Y_sub = Kokkos::subview (Y, Kokkos::ALL (), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy (Y, ZERO_d);
      KokkosSparse::spmv ("N", alpha_d, A_f, X, beta_d, Y);

      Kokkos::deep_copy (Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0,j);
        EXPECT_EQ( actualResult, expectedResult_mixed );
      }
    }
  }
}

#define EXECUTE_TEST_ISSUE_101( DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv_issue_101 ## _ ## OFFSET ## _ ## DEVICE ) { \
	test_github_issue_101<DEVICE> (); \
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

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  EXECUTE_TEST_ISSUE_101(TestExecSpace)
#endif



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


