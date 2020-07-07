#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>

#include<KokkosSparse_spmv.hpp>
#include<KokkosKernels_TestUtils.hpp>
#include<KokkosKernels_Test_Structured_Matrix.hpp>
#include<KokkosKernels_IOUtils.hpp>
#include<KokkosKernels_Utils.hpp>

#include "KokkosKernels_Controls.hpp"

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test {

  template < class VectorType0, class VectorType1 >
struct fSPMV {
  using value_type = int;
  using AT         = Kokkos::ArithTraits<typename VectorType1::non_const_value_type>;
  using ATM        = Kokkos::ArithTraits<typename AT::mag_type>;
  using mag_type   = typename AT::mag_type;

  VectorType0 expected_y;
  VectorType1 y;
  mag_type    eps;

  fSPMV(const VectorType0 & _ex_y, const VectorType1 & _y, const mag_type _eps)
  : expected_y(_ex_y)
  , y(_y)
  , eps(_eps)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, value_type& err ) const {

    const mag_type error = AT::abs(expected_y(i)-y(i))
      / (AT::abs(expected_y(i)) > ATM::zero() ? AT::abs(expected_y(i)) : ATM::one());

    if(error > eps) {
      err++;
      printf("expected_y(%d)=%f, y(%d)=%f\n", i, AT::abs(expected_y(i)), i, AT::abs(y(i)));
    }
  }
};


template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void sequential_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
                     typename y_vector_type::non_const_value_type alpha,
                     typename y_vector_type::non_const_value_type beta){

  using graph_t          = typename crsMat_t::StaticCrsGraphType;
  using size_type_view_t = typename graph_t::row_map_type;
  using lno_view_t       = typename graph_t::entries_type;
  using scalar_view_t    = typename crsMat_t::values_type::non_const_type;

  using size_type = typename size_type_view_t::non_const_value_type;
  using lno_t     = typename lno_view_t::non_const_value_type;
  using scalar_t  = typename scalar_view_t::non_const_value_type;


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
                typename y_vector_type::non_const_value_type alpha,
                typename y_vector_type::non_const_value_type beta) {
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  using ExecSpace = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps = std::is_same<y_value_mag_type, float>::value ? 2*1e-3 : 1e-7;
  const size_t nr = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);
  //KokkosKernels::Impl::print_1Dview(expected_y);
  KokkosSparse::spmv("N", alpha, input_mat, x, beta, y);
  //KokkosKernels::Impl::print_1Dview(y);
  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosSparse::Test::spmv",
                          my_exec_space(0, y.extent(0)),
                          fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps),
                          num_errors);
  if(num_errors>0) printf("KokkosSparse::Test::spmv: %i errors of %i with params: %lf %lf\n",
                          num_errors, y.extent_int(0),
                          y_value_trait::abs(alpha), y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors==0);
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type y, y_vector_type expected_y,
                   typename y_vector_type::non_const_value_type alpha,
                   typename y_vector_type::non_const_value_type beta, int numMV) {
  using ExecSpace = typename crsMat_t::execution_space;
  using my_exec_space = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps = std::is_same<y_value_mag_type, float>::value ? 2*1e-3 : 1e-7;

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
    Kokkos::parallel_reduce("KokkosSparse::Test::spmv_mv",
                            my_exec_space(0,y_i.extent(0)),
                            fSPMV<decltype(y_i), decltype(y_spmv)>(y_i, y_spmv, eps),
                            num_errors);
    if(num_errors>0) printf("KokkosSparse::Test::spmv_mv: %i errors of %i for mv %i\n",
                            num_errors, y_i.extent_int(0), i);
    EXPECT_TRUE(num_errors==0);
  }
}

template <typename crsMat_t,
          typename x_vector_type,
          typename y_vector_type>
void check_spmv_struct(const crsMat_t input_mat,
                       const int stencil_type,
                       const Kokkos::View<typename crsMat_t::non_const_ordinal_type*, Kokkos::HostSpace> structure,
                       x_vector_type x,
                       y_vector_type y,
                       typename y_vector_type::non_const_value_type alpha,
                       typename y_vector_type::non_const_value_type beta) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const double eps = std::is_same<y_value_mag_type, float>::value ? 2*1e-3 : 1e-7;
  const size_t nr = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);
  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure,
                                          alpha, input_mat, x, beta, y);

  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv_struct",
                          my_exec_space(0, y.extent(0)),
                          fSPMV<y_vector_type, y_vector_type>(expected_y,y,eps),
                          num_errors);
  if(num_errors>0) printf("KokkosKernels::UnitTests::spmv_struct: %i errors of %i with params: %d %lf %lf\n",
                          num_errors, y.extent_int(0), stencil_type,
                          y_value_trait::abs(alpha), y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors==0);
} // check_spmv_struct

template <typename crsMat_t,
          typename x_vector_type,
          typename y_vector_type>
void check_spmv_mv_struct(const crsMat_t input_mat,
                          const int stencil_type,
                          const Kokkos::View<typename crsMat_t::non_const_ordinal_type*, Kokkos::HostSpace> structure,
                          x_vector_type x,
                          y_vector_type y,
                          y_vector_type expected_y,
                          typename y_vector_type::non_const_value_type alpha,
                          typename y_vector_type::non_const_value_type beta,
                          int numMV) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const double eps = std::is_same<y_value_mag_type, float>::value ? 2*1e-3 : 1e-7;
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure,
                                          alpha, input_mat, x, beta, y);

  for(int vectorIdx = 0; vectorIdx < numMV; ++vectorIdx) {
    auto x_i = Kokkos::subview (x, Kokkos::ALL (), vectorIdx);
    auto y_i = Kokkos::subview (expected_y, Kokkos::ALL (), vectorIdx);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta);

    auto y_spmv = Kokkos::subview (y, Kokkos::ALL (), vectorIdx);
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv_mv_struct",
                            my_exec_space(0, y.extent(0)),
                            fSPMV<decltype(y_i), decltype(y_spmv)>(y_i,y_spmv,eps),
                            num_errors);
    if(num_errors>0) printf("KokkosKernels::UnitTests::spmv_mv_struct: %i errors of %i with params: %d %lf %lf, in vector %i\n",
                            num_errors, y.extent_int(0), stencil_type,
                            y_value_trait::abs(alpha), y_value_trait::abs(beta), vectorIdx);
    EXPECT_TRUE(num_errors==0);
  }
} // check_spmv_mv_struct

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_controls(KokkosKernels::Experimental::Controls controls,
			 crsMat_t input_mat, x_vector_type x, y_vector_type y,
			 typename y_vector_type::non_const_value_type alpha,
			 typename y_vector_type::non_const_value_type beta) {
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  using ExecSpace = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps = std::is_same<y_value_mag_type, float>::value ? 2*1e-3 : 1e-7;
  const size_t nr = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  controls.setParameter("algorithm", "merge");
  printf("requested merge based algorithm\n");
#endif

  KokkosSparse::spmv(controls, "N", alpha, input_mat, x, beta, y);
  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosSparse::Test::spmv",
                          my_exec_space(0, y.extent(0)),
                          fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps),
                          num_errors);
  if(num_errors>0) printf("KokkosSparse::Test::spmv: %i errors of %i with params: %lf %lf\n",
                          num_errors, y.extent_int(0),
                          y_value_trait::abs(alpha), y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors==0);
} // check_spmv_controls

} // namespace Test


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

  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t**, layout, Device> ViewTypeY;

  ViewTypeX b_x("A",numRows,numMV);
  ViewTypeY b_y("B",numCols,numMV);
  ViewTypeY b_y_copy("B",numCols,numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);
  Kokkos::fill_random(b_x,rand_pool,scalar_t(10));
  Kokkos::fill_random(b_y,rand_pool,scalar_t(10));


  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);

  Kokkos::deep_copy(b_y_copy, b_y);


  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 0.0, numMV);
  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 0.0, 1.0, numMV);
  Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 1.0, numMV);


}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_1D(lno_t nx, lno_t leftBC, lno_t rightBC) {

  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;

  Kokkos::View<lno_t*, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t*[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 1);
  mat_structure(0, 0) = nx;
  if(leftBC  == 1) { mat_structure(0, 1) = 1; }
  if(rightBC == 1) { mat_structure(0, 2) = 1; }

  crsMat_t input_mat = Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x  ("x", nc);
  y_vector_type output_y ("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x,  rand_pool, ScalarX(10));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(10));

  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_2D(lno_t nx, lno_t ny, lno_t horizontalBC, lno_t verticalBC) {

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type  scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  Kokkos::View<lno_t*, Kokkos::HostSpace> structure("Spmv Structure", 2);
  structure(0) = nx;
  structure(1) = ny;
  Kokkos::View<lno_t*[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 2);
  mat_structure(0, 0) = nx;
  if(horizontalBC == 1 || horizontalBC == 3) { mat_structure(0, 1) = 1; }
  if(horizontalBC == 2 || horizontalBC == 3) { mat_structure(0, 2) = 1; }
  mat_structure(1, 0) = ny;
  if(verticalBC == 1 || verticalBC == 3) { mat_structure(1, 1) = 1; }
  if(verticalBC == 2 || verticalBC == 3) { mat_structure(1, 2) = 1; }

  crsMat_t input_mat_FD = Test::generate_structured_matrix2D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE = Test::generate_structured_matrix2D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x ("x", nc);
  y_vector_type output_y ("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x,rand_pool,ScalarX(10));
  Kokkos::fill_random(output_y,rand_pool,ScalarY(10));

  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 1.0);

  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_3D(lno_t nx, lno_t ny, lno_t nz, lno_t horizontal1BC, lno_t horizontal2BC, lno_t verticalBC) {

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type  scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  Kokkos::View<lno_t*, Kokkos::HostSpace> structure("Spmv Structure", 3);
  structure(0) = nx;
  structure(1) = ny;
  structure(2) = nz;
  Kokkos::View<lno_t*[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 3);
  mat_structure(0, 0) = nx;
  if(horizontal1BC == 1 || horizontal1BC == 3) { mat_structure(0, 1) = 1; }
  if(horizontal1BC == 2 || horizontal1BC == 3) { mat_structure(0, 2) = 1; }
  mat_structure(1, 0) = ny;
  if(horizontal2BC == 1 || horizontal2BC == 3) { mat_structure(1, 1) = 1; }
  if(horizontal2BC == 2 || horizontal2BC == 3) { mat_structure(1, 2) = 1; }
  mat_structure(2, 0) = nz;
  if(verticalBC == 1 || verticalBC == 3) { mat_structure(2, 1) = 1; }
  if(verticalBC == 2 || verticalBC == 3) { mat_structure(2, 2) = 1; }

  crsMat_t input_mat_FD = Test::generate_structured_matrix3D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE = Test::generate_structured_matrix3D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x  ("x", nc);
  y_vector_type output_y ("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x,rand_pool,ScalarX(10));
  Kokkos::fill_random(output_y,rand_pool,ScalarY(10));

  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 1.0);

  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, class Device>
void test_spmv_mv_struct_1D(lno_t nx, int numMV) {

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;
  typedef Kokkos::View<scalar_t**, layout, Device> x_multivector_type;
  typedef Kokkos::View<scalar_t**, layout, Device> y_multivector_type;

  Kokkos::View<lno_t*, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t*[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 1);
  mat_structure(0, 0) = nx;
  mat_structure(0, 1) = 1;
  mat_structure(0, 2) = 1;

  crsMat_t input_mat = Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_multivector_type input_x  ("x", nc, numMV);
  y_multivector_type output_y ("y", nr, numMV);
  y_multivector_type output_y_copy ("y_copy", nr, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  typedef typename x_multivector_type::value_type ScalarX;
  typedef typename y_multivector_type::value_type ScalarY;

  Kokkos::fill_random(input_x,  rand_pool, ScalarX(10));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(10));

  Kokkos::deep_copy(output_y_copy, output_y);

  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 1.0, 0.0, numMV);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 0.0, 1.0, numMV);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 1.0, 1.0, numMV);
}

// check that the controls are flowing down correctly in the spmv kernel
template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_controls(lno_t numRows,size_type nnz, lno_t bandwidth, lno_t row_size_variance) {

  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using Controls      = KokkosKernels::Experimental::Controls;


  lno_t numCols = numRows;

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);
  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x ("x", nc);
  y_vector_type output_y ("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  using ScalarX = typename x_vector_type::value_type;
  using ScalarY = typename y_vector_type::value_type;

  Kokkos::fill_random(input_x,rand_pool,ScalarX(10));
  Kokkos::fill_random(output_y,rand_pool,ScalarY(10));

  Controls controls;

  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 1.0, 1.0);
} // test_spmv_controls

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
  test_spmv_controls<SCALAR,ORDINAL,OFFSET,DEVICE> (10000, 10000 * 20, 100, 5); \
}

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv_mv ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## LAYOUT ## _ ## DEVICE ) { \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (50000, 50000 * 30, 100, 10, 5); \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (50000, 50000 * 30, 200, 10, 1); \
  test_spmv_mv<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (10000, 10000 * 20, 100, 5, 10); \
}

#define EXECUTE_TEST_STRUCT(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv_struct ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_spmv_struct_1D<SCALAR,ORDINAL,OFFSET,DEVICE> (10, 1, 1);            \
  test_spmv_struct_2D<SCALAR,ORDINAL,OFFSET,DEVICE> (250, 201, 3, 3);      \
  test_spmv_struct_2D<SCALAR,ORDINAL,OFFSET,DEVICE> (200, 250, 3, 3);      \
  test_spmv_struct_2D<SCALAR,ORDINAL,OFFSET,DEVICE> (251, 251, 3, 3);      \
  test_spmv_struct_3D<SCALAR,ORDINAL,OFFSET,DEVICE> (30, 30, 30, 3, 3, 3); \
  test_spmv_struct_3D<SCALAR,ORDINAL,OFFSET,DEVICE> (40, 40, 40, 3, 3, 3); \
  test_spmv_struct_3D<SCALAR,ORDINAL,OFFSET,DEVICE> (25, 40, 50, 3, 3, 3); \
  test_spmv_struct_3D<SCALAR,ORDINAL,OFFSET,DEVICE> (40, 50, 25, 3, 3, 3); \
  test_spmv_struct_3D<SCALAR,ORDINAL,OFFSET,DEVICE> (50, 24, 40, 3, 3, 3); \
}


#define EXECUTE_TEST_MV_STRUCT(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spmv_mv_struct ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## LAYOUT ## _ ## DEVICE ) { \
  test_spmv_mv_struct_1D<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (10, 1); \
  test_spmv_mv_struct_1D<SCALAR,ORDINAL,OFFSET,Kokkos::LAYOUT,DEVICE> (10, 2); \
}

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  EXECUTE_TEST_ISSUE_101(TestExecSpace)
#endif



#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
 EXECUTE_TEST_STRUCT(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif



#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(double, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(double, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(double, int64_t, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(double, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(float, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(float, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(float, int64_t, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(float, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_double, int64_t, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, int, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST_MV(kokkos_complex_float, int64_t, size_t, LayoutLeft, TestExecSpace)
 EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int64_t, size_t, LayoutLeft, TestExecSpace)
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


