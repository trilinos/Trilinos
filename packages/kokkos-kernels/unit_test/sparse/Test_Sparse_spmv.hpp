#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_Test_Structured_Matrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosKernels_Utils.hpp>

#include "KokkosKernels_Controls.hpp"
#include "KokkosKernels_default_types.hpp"

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;
typedef Kokkos::Experimental::half_t kokkos_half;

namespace Test {

template <class VectorType0, class VectorType1>
struct fSPMV {
  using value_type = int;
  using AT  = Kokkos::ArithTraits<typename VectorType1::non_const_value_type>;
  using ATM = Kokkos::ArithTraits<typename AT::mag_type>;
  using mag_type = typename AT::mag_type;

  VectorType0 expected_y;
  VectorType1 y;
  mag_type eps;

  fSPMV(const VectorType0 &_ex_y, const VectorType1 &_y, const mag_type _eps)
      : expected_y(_ex_y), y(_y), eps(_eps) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &err) const {
    const mag_type error =
        AT::abs(expected_y(i) - y(i)) / (AT::abs(expected_y(i)) > ATM::zero()
                                             ? AT::abs(expected_y(i))
                                             : ATM::one());

    if (error > eps) {
      err++;
      // printf("expected_y(%d)=%f, y(%d)=%f err=%f, eps=%f\n", i,
      //        AT::abs(expected_y(i)), i, AT::abs(y(i)), error, eps);
    }
  }
};

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void sequential_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
                     typename y_vector_type::non_const_value_type alpha,
                     typename y_vector_type::non_const_value_type beta,
                     char mode = 'N') {
  using graph_t          = typename crsMat_t::StaticCrsGraphType;
  using size_type_view_t = typename graph_t::row_map_type;
  using lno_view_t       = typename graph_t::entries_type;
  using scalar_view_t    = typename crsMat_t::values_type::non_const_type;

  using size_type = typename size_type_view_t::non_const_value_type;
  using lno_t     = typename lno_view_t::non_const_value_type;
  using scalar_t  = typename scalar_view_t::non_const_value_type;
  using KAT       = Kokkos::ArithTraits<scalar_t>;

  mode = toupper(mode);

  typename scalar_view_t::HostMirror h_values =
      Kokkos::create_mirror_view(input_mat.values);
  Kokkos::deep_copy(h_values, input_mat.values);

  typename lno_view_t::HostMirror h_entries =
      Kokkos::create_mirror_view(input_mat.graph.entries);
  Kokkos::deep_copy(h_entries, input_mat.graph.entries);

  typename size_type_view_t::HostMirror h_rowmap =
      Kokkos::create_mirror_view(input_mat.graph.row_map);
  Kokkos::deep_copy(h_rowmap, input_mat.graph.row_map);
  Kokkos::fence();

  typename y_vector_type::HostMirror h_y = Kokkos::create_mirror_view(y);
  typename x_vector_type::HostMirror h_x = Kokkos::create_mirror_view(x);

  KokkosKernels::Impl::safe_device_to_host_deep_copy(x.extent(0), x, h_x);
  KokkosKernels::Impl::safe_device_to_host_deep_copy(y.extent(0), y, h_y);
  Kokkos::fence();

  lno_t nr = input_mat.numRows();

  // first, scale y by beta
  for (size_t i = 0; i < h_y.extent(0); i++) h_y(i) *= beta;

  // then go through the matrix and accumulate the matrix-vector product
  for (lno_t row = 0; row < nr; ++row) {
    for (size_type j = h_rowmap(row); j < h_rowmap(row + 1); ++j) {
      lno_t col    = h_entries(j);
      scalar_t val = h_values(j);
      if (mode == 'N')
        h_y(row) += alpha * val * h_x(col);
      else if (mode == 'C')
        h_y(row) += alpha * KAT::conj(val) * h_x(col);
      else if (mode == 'T')
        h_y(col) += alpha * val * h_x(row);
      else if (mode == 'H')
        h_y(col) += alpha * KAT::conj(val) * h_x(row);
    }
  }
  KokkosKernels::Impl::safe_host_to_device_deep_copy(y.extent(0), h_y, y);
  Kokkos::fence();
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
                typename y_vector_type::non_const_value_type alpha,
                typename y_vector_type::non_const_value_type beta, char mode) {
  // typedef typename crsMat_t::StaticCrsGraphType graph_t;
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps =
      std::is_same<y_value_mag_type, float>::value ? 2 * 1e-3 : 1e-7;
  bool transposed = (mode == 'T') || (mode == 'H');
  y_vector_type expected_y(
      "expected", transposed ? input_mat.numCols() : input_mat.numRows());
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta, mode);
  bool threw = false;
  std::string msg;
  try {
    KokkosSparse::spmv(&mode, alpha, input_mat, x, beta, y);
    Kokkos::fence();
  } catch (std::exception &e) {
    threw = true;
    msg   = e.what();
  }
  ASSERT_FALSE(threw) << "KokkosSparse::Test::spmv 1D, mode " << mode
                      << ": threw exception:\n"
                      << msg << '\n';
  int num_errors = 0;
  Kokkos::parallel_reduce(
      "KokkosSparse::Test::spmv", my_exec_space(0, y.extent(0)),
      fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps), num_errors);
  if (num_errors > 0)
    printf("KokkosSparse::Test::spmv: %i errors of %i with params: %lf %lf\n",
           num_errors, y.extent_int(0), y_value_trait::abs(alpha),
           y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors == 0);
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
                   y_vector_type expected_y,
                   typename y_vector_type::non_const_value_type alpha,
                   typename y_vector_type::non_const_value_type beta, int numMV,
                   char mode) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps =
      std::is_same<y_value_mag_type, float>::value ? 2 * 1e-3 : 1e-7;

  Kokkos::deep_copy(expected_y, y);

  Kokkos::fence();

  bool threw = false;
  std::string msg;
  try {
    KokkosSparse::spmv(&mode, alpha, input_mat, x, beta, y);
    Kokkos::fence();
  } catch (std::exception &e) {
    threw = true;
    msg   = e.what();
  }
  ASSERT_FALSE(threw) << "KokkosSparse::Test::spmv 2D, mode " << mode
                      << ": threw exception:\n"
                      << msg << '\n';

  for (int i = 0; i < numMV; ++i) {
    auto x_i = Kokkos::subview(x, Kokkos::ALL(), i);

    auto y_i = Kokkos::subview(expected_y, Kokkos::ALL(), i);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta, mode);

    auto y_spmv    = Kokkos::subview(y, Kokkos::ALL(), i);
    int num_errors = 0;
    Kokkos::parallel_reduce(
        "KokkosSparse::Test::spmv_mv", my_exec_space(0, y_i.extent(0)),
        fSPMV<decltype(y_i), decltype(y_spmv)>(y_i, y_spmv, eps), num_errors);
    if (num_errors > 0)
      std::cout << "KokkosSparse::Test::spmv_mv: " << num_errors
                << " errors of " << y_i.extent_int(0) << " for mv " << i
                << " (alpha=" << alpha << ", beta=" << beta
                << ", mode = " << mode << ")\n";
    EXPECT_TRUE(num_errors == 0);
  }
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_struct(
    const crsMat_t input_mat, const int stencil_type,
    const Kokkos::View<typename crsMat_t::non_const_ordinal_type *,
                       Kokkos::HostSpace>
        structure,
    x_vector_type x, y_vector_type y,
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
  const double eps =
      std::is_same<y_value_mag_type, float>::value ? 2 * 1e-3 : 1e-7;
  const size_t nr = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);
  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, alpha,
                                          input_mat, x, beta, y);

  int num_errors = 0;
  Kokkos::parallel_reduce(
      "KokkosKernels::UnitTests::spmv_struct", my_exec_space(0, y.extent(0)),
      fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps), num_errors);
  if (num_errors > 0)
    printf(
        "KokkosKernels::UnitTests::spmv_struct: %i errors of %i with params: "
        "%d %lf %lf\n",
        num_errors, y.extent_int(0), stencil_type, y_value_trait::abs(alpha),
        y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors == 0);
}  // check_spmv_struct

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv_struct(
    const crsMat_t input_mat, const int stencil_type,
    const Kokkos::View<typename crsMat_t::non_const_ordinal_type *,
                       Kokkos::HostSpace>
        structure,
    x_vector_type x, y_vector_type y, y_vector_type expected_y,
    typename y_vector_type::non_const_value_type alpha,
    typename y_vector_type::non_const_value_type beta, int numMV) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const double eps =
      std::is_same<y_value_mag_type, float>::value ? 2 * 1e-3 : 1e-7;
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, alpha,
                                          input_mat, x, beta, y);

  for (int vectorIdx = 0; vectorIdx < numMV; ++vectorIdx) {
    auto x_i = Kokkos::subview(x, Kokkos::ALL(), vectorIdx);
    auto y_i = Kokkos::subview(expected_y, Kokkos::ALL(), vectorIdx);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta);

    auto y_spmv    = Kokkos::subview(y, Kokkos::ALL(), vectorIdx);
    int num_errors = 0;
    Kokkos::parallel_reduce(
        "KokkosKernels::UnitTests::spmv_mv_struct",
        my_exec_space(0, y.extent(0)),
        fSPMV<decltype(y_i), decltype(y_spmv)>(y_i, y_spmv, eps), num_errors);
    if (num_errors > 0)
      printf(
          "KokkosKernels::UnitTests::spmv_mv_struct: %i errors of %i with "
          "params: %d %lf %lf, in vector %i\n",
          num_errors, y.extent_int(0), stencil_type, y_value_trait::abs(alpha),
          y_value_trait::abs(beta), vectorIdx);
    EXPECT_TRUE(num_errors == 0);
  }
}  // check_spmv_mv_struct

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_controls(KokkosKernels::Experimental::Controls controls,
                         crsMat_t input_mat, x_vector_type x, y_vector_type y,
                         typename y_vector_type::non_const_value_type alpha,
                         typename y_vector_type::non_const_value_type beta) {
  // typedef typename crsMat_t::StaticCrsGraphType graph_t;
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps =
      std::is_same<y_value_mag_type, float>::value ? 2 * 1e-3 : 1e-7;
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
  Kokkos::parallel_reduce(
      "KokkosSparse::Test::spmv", my_exec_space(0, y.extent(0)),
      fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps), num_errors);
  if (num_errors > 0)
    printf("KokkosSparse::Test::spmv: %i errors of %i with params: %lf %lf\n",
           num_errors, y.extent_int(0), y_value_trait::abs(alpha),
           y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors == 0);
}  // check_spmv_controls

}  // namespace Test

template <typename scalar_t>
scalar_t randomUpperBound(int mag) {
  return (scalar_t)mag;
}

template <>
Kokkos::complex<double> randomUpperBound<Kokkos::complex<double>>(int mag) {
  return Kokkos::complex<double>(mag, mag);
}

template <>
Kokkos::complex<float> randomUpperBound<Kokkos::complex<float>>(int mag) {
  return Kokkos::complex<float>(mag, mag);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv(lno_t numRows, size_type nnz, lno_t bandwidth,
               lno_t row_size_variance, bool heavy) {
  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  lno_t numCols = numRows;

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);
  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);
  x_vector_type input_xt("x", nr);
  y_vector_type output_yt("y", nc);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x, rand_pool, randomUpperBound<ScalarX>(1));
  Kokkos::fill_random(output_y, rand_pool, randomUpperBound<ScalarY>(1));
  Kokkos::fill_random(input_xt, rand_pool, randomUpperBound<ScalarX>(1));
  Kokkos::fill_random(output_yt, rand_pool, randomUpperBound<ScalarY>(1));

  std::vector<char> nonTransModes   = {'N'};
  std::vector<char> transModes      = {'T'};
  std::vector<double> testAlphaBeta = {0.0, 1.0};
  if (heavy) {
    nonTransModes.push_back('C');
    transModes.push_back('H');
    testAlphaBeta.push_back(-1.0);
    testAlphaBeta.push_back(2.5);
  }
  for (auto mode : nonTransModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        Test::check_spmv(input_mat, input_x, output_y, alpha, beta, mode);
      }
    }
  }
  for (auto mode : transModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        Test::check_spmv(input_mat, input_xt, output_yt, alpha, beta, mode);
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename layout, class Device>
void test_spmv_mv(lno_t numRows, size_type nnz, lno_t bandwidth,
                  lno_t row_size_variance, bool heavy, int numMV) {
  lno_t numCols = numRows;

  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;

  typedef Kokkos::View<scalar_t **, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t **, layout, Device> ViewTypeY;

  ViewTypeX b_x("A", numRows, numMV);
  ViewTypeY b_y("B", numCols, numMV);
  ViewTypeY b_y_copy("B", numCols, numMV);

  ViewTypeX b_xt("A", numCols, numMV);
  ViewTypeY b_yt("B", numRows, numMV);
  ViewTypeY b_yt_copy("B", numRows, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);
  Kokkos::fill_random(b_x, rand_pool, randomUpperBound<scalar_t>(1));
  Kokkos::fill_random(b_y, rand_pool, randomUpperBound<scalar_t>(1));
  Kokkos::fill_random(b_xt, rand_pool, randomUpperBound<scalar_t>(1));
  Kokkos::fill_random(b_yt, rand_pool, randomUpperBound<scalar_t>(1));

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);

  Kokkos::deep_copy(b_y_copy, b_y);
  Kokkos::deep_copy(b_yt_copy, b_yt);

  std::vector<char> nonTransModes   = {'N'};
  std::vector<char> transModes      = {'T'};
  std::vector<double> testAlphaBeta = {0.0, 1.0};
  if (heavy) {
    nonTransModes.push_back('C');
    transModes.push_back('H');
    testAlphaBeta.push_back(-1.0);
    testAlphaBeta.push_back(2.5);
  }
  for (auto mode : nonTransModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, alpha, beta, numMV,
                            mode);
      }
    }
  }
  for (auto mode : transModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        Test::check_spmv_mv(input_mat, b_xt, b_yt, b_yt_copy, alpha, beta,
                            numMV, mode);
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename layout, class Device>
void test_spmv_mv_heavy(lno_t numRows, size_type nnz, lno_t bandwidth,
                        lno_t row_size_variance, int numMV) {
  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;

  typedef Kokkos::View<scalar_t **, layout, Device> ViewTypeX;
  typedef Kokkos::View<scalar_t **, layout, Device> ViewTypeY;

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numRows, nnz, row_size_variance, bandwidth);
  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  for (int nv = 1; nv <= numMV; nv++) {
    ViewTypeX b_x("A", numRows, nv);
    ViewTypeY b_y("B", numRows, nv);
    ViewTypeY b_y_copy("B", numRows, nv);

    Kokkos::fill_random(b_x, rand_pool, scalar_t(10));
    Kokkos::fill_random(b_y, rand_pool, scalar_t(10));

    Kokkos::deep_copy(b_y_copy, b_y);

    Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 0.0, nv, 'N');
    Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 0.0, 1.0, nv, 'N');
    Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 1.0, nv, 'N');
    Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 1.0, 0.0, nv, 'T');
    Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, 0.0, 1.0, nv, 'T');
    // Testing all modes together, since matrix is square
    std::vector<char> modes           = {'N', 'C', 'T', 'H'};
    std::vector<double> testAlphaBeta = {0.0, 1.0, -1.0, 2.5};
    for (auto mode : modes) {
      for (double alpha : testAlphaBeta) {
        for (double beta : testAlphaBeta) {
          Test::check_spmv_mv(input_mat, b_x, b_y, b_y_copy, alpha, beta, nv,
                              mode);
        }
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_1D(lno_t nx, lno_t leftBC, lno_t rightBC) {
  using crsMat_t = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device,
                                                    void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t * [3], Kokkos::HostSpace> mat_structure("Matrix Structure",
                                                             1);
  mat_structure(0, 0) = nx;
  if (leftBC == 1) {
    mat_structure(0, 1) = 1;
  }
  if (rightBC == 1) {
    mat_structure(0, 2) = 1;
  }

  crsMat_t input_mat =
      Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x, rand_pool, ScalarX(1));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(1));

  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_2D(lno_t nx, lno_t ny, lno_t horizontalBC,
                         lno_t verticalBC) {
  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 2);
  structure(0) = nx;
  structure(1) = ny;
  Kokkos::View<lno_t * [3], Kokkos::HostSpace> mat_structure("Matrix Structure",
                                                             2);
  mat_structure(0, 0) = nx;
  if (horizontalBC == 1 || horizontalBC == 3) {
    mat_structure(0, 1) = 1;
  }
  if (horizontalBC == 2 || horizontalBC == 3) {
    mat_structure(0, 2) = 1;
  }
  mat_structure(1, 0) = ny;
  if (verticalBC == 1 || verticalBC == 3) {
    mat_structure(1, 1) = 1;
  }
  if (verticalBC == 2 || verticalBC == 3) {
    mat_structure(1, 2) = 1;
  }

  crsMat_t input_mat_FD =
      Test::generate_structured_matrix2D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE =
      Test::generate_structured_matrix2D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x, rand_pool, ScalarX(1));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(1));

  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0,
                          0.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0,
                          1.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0,
                          1.0);

  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0,
                          0.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0,
                          1.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0,
                          1.0);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_3D(lno_t nx, lno_t ny, lno_t nz, lno_t horizontal1BC,
                         lno_t horizontal2BC, lno_t verticalBC) {
  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 3);
  structure(0) = nx;
  structure(1) = ny;
  structure(2) = nz;
  Kokkos::View<lno_t * [3], Kokkos::HostSpace> mat_structure("Matrix Structure",
                                                             3);
  mat_structure(0, 0) = nx;
  if (horizontal1BC == 1 || horizontal1BC == 3) {
    mat_structure(0, 1) = 1;
  }
  if (horizontal1BC == 2 || horizontal1BC == 3) {
    mat_structure(0, 2) = 1;
  }
  mat_structure(1, 0) = ny;
  if (horizontal2BC == 1 || horizontal2BC == 3) {
    mat_structure(1, 1) = 1;
  }
  if (horizontal2BC == 2 || horizontal2BC == 3) {
    mat_structure(1, 2) = 1;
  }
  mat_structure(2, 0) = nz;
  if (verticalBC == 1 || verticalBC == 3) {
    mat_structure(2, 1) = 1;
  }
  if (verticalBC == 2 || verticalBC == 3) {
    mat_structure(2, 2) = 1;
  }

  crsMat_t input_mat_FD =
      Test::generate_structured_matrix3D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE =
      Test::generate_structured_matrix3D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  typedef typename x_vector_type::value_type ScalarX;
  typedef typename y_vector_type::value_type ScalarY;

  Kokkos::fill_random(input_x, rand_pool, ScalarX(1));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(1));

  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0,
                          0.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0,
                          1.0);
  Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0,
                          1.0);

  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0,
                          0.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0,
                          1.0);
  Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0,
                          1.0);
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename layout, class Device>
void test_spmv_mv_struct_1D(lno_t nx, int numMV) {
  typedef
      typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>
          crsMat_t;
  typedef Kokkos::View<scalar_t **, layout, Device> x_multivector_type;
  typedef Kokkos::View<scalar_t **, layout, Device> y_multivector_type;

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t * [3], Kokkos::HostSpace> mat_structure("Matrix Structure",
                                                             1);
  mat_structure(0, 0) = nx;
  mat_structure(0, 1) = 1;
  mat_structure(0, 2) = 1;

  crsMat_t input_mat =
      Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_multivector_type input_x("x", nc, numMV);
  y_multivector_type output_y("y", nr, numMV);
  y_multivector_type output_y_copy("y_copy", nr, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  typedef typename x_multivector_type::value_type ScalarX;
  typedef typename y_multivector_type::value_type ScalarY;

  Kokkos::fill_random(input_x, rand_pool, ScalarX(10));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(10));

  Kokkos::deep_copy(output_y_copy, output_y);

  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y,
                             output_y_copy, 1.0, 0.0, numMV);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y,
                             output_y_copy, 0.0, 1.0, numMV);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y,
                             output_y_copy, 1.0, 1.0, numMV);
}

// check that the controls are flowing down correctly in the spmv kernel
template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_controls(lno_t numRows, size_type nnz, lno_t bandwidth,
                        lno_t row_size_variance) {
  using crsMat_t = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device,
                                                    void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using Controls      = KokkosKernels::Experimental::Controls;

  lno_t numCols = numRows;

  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);
  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  using ScalarX = typename x_vector_type::value_type;
  using ScalarY = typename y_vector_type::value_type;

  Kokkos::fill_random(input_x, rand_pool, ScalarX(10));
  Kokkos::fill_random(output_y, rand_pool, ScalarY(10));

  Controls controls;

  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 1.0, 0.0);
  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 0.0, 1.0);
  Test::check_spmv_controls(controls, input_mat, input_x, output_y, 1.0, 1.0);
}  // test_spmv_controls

// call it if ordinal int and, scalar float and double are instantiated.
template <class DeviceType>
void test_github_issue_101() {
  typedef KokkosSparse::CrsMatrix<float, int, DeviceType> float_matrix_type;
  typedef KokkosSparse::CrsMatrix<double, int, DeviceType> double_matrix_type;
  static_assert(
      std::is_same<typename float_matrix_type::StaticCrsGraphType,
                   typename double_matrix_type::StaticCrsGraphType>::value,
      "Two KokkosSparse::CrsMatrix types that differ only in the type of "
      "matrix values, appear to have two different StaticCrsGraphType "
      "typedefs.  This should never happen.");
  typedef typename float_matrix_type::StaticCrsGraphType graph_type;

  constexpr int numRows    = 1;
  constexpr int numCols    = 2;
  constexpr double alpha_d = 1.0;
  constexpr double beta_d  = 0.0;
  const float EPS_f        = std::numeric_limits<float>::epsilon();

  graph_type G;
  {
    typename graph_type::entries_type colInds("colInds", numCols);
    auto colInds_h = Kokkos::create_mirror_view(colInds);
    colInds_h[0]   = 0;
    colInds_h[1]   = 1;
    Kokkos::deep_copy(colInds, colInds_h);

    typedef typename graph_type::row_map_type::non_const_type row_offsets_type;
    row_offsets_type rowOffsets("rowOffsets", numRows + 1);
    auto rowOffsets_h = Kokkos::create_mirror_view(rowOffsets);
    rowOffsets_h[0]   = 0;  // Entries start at offset 0
    rowOffsets_h[1]   = 2;  // 2 entries total in the "sparse" matrix
    Kokkos::deep_copy(rowOffsets, rowOffsets_h);

    G = graph_type(colInds, rowOffsets);
  }

  Kokkos::View<double *, DeviceType> x("x", numCols);
  Kokkos::deep_copy(x, static_cast<double>(1.0));
  Kokkos::View<double *, DeviceType> y("y", numRows);
  auto y_h = Kokkos::create_mirror_view(y);  // we'll want this later

  // Pick some number large enough to exercise all unrolling cases.
  // Sparse mat-vec does or at least used to unroll for 1, 2, ..., 17
  // vectors.  Include a little extra in case the implementers decide
  // to strip-mine that.
  constexpr int numVecs = 22;
  Kokkos::View<double **, default_layout, DeviceType> X("X", numCols, numVecs);
  Kokkos::deep_copy(X, static_cast<double>(1.0));
  Kokkos::View<double **, default_layout, DeviceType> Y("Y", numRows, numVecs);
  auto Y_h = Kokkos::create_mirror_view(Y);  // we'll want this later

  // Start with the easy test case, where the matrix and the vectors
  // are all double.
  {
    constexpr double ZERO_d = static_cast<double>(0.0);
    constexpr double ONE_d  = static_cast<double>(1.0);
    constexpr double TWO_d  = static_cast<double>(2.0);

    double_matrix_type A_d("A_d", G);
    auto A_d_val_h = Kokkos::create_mirror_view(A_d.values);
    A_d_val_h[0]   = ONE_d;
    // This cast is deliberate; we want to use float eps here, but as
    // a double-precision number.  This is just a sanity check for
    // accuracy of the sparse mat-vec when not using mixed precision.
    A_d_val_h[1] = static_cast<double>(EPS_f) / TWO_d;
    EXPECT_NE(A_d_val_h[1], ZERO_d);  // just making sure
    Kokkos::deep_copy(A_d.values, A_d_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy(y, ZERO_d);
    KokkosSparse::spmv("N", alpha_d, A_d, x, beta_d, y);

    Kokkos::deep_copy(y_h, y);
    const double expectedResult_allDouble =
        static_cast<double>(1.0) +
        static_cast<double>(EPS_f) / static_cast<double>(2.0);
    EXPECT_NE(expectedResult_allDouble, ZERO_d);
    EXPECT_EQ(y_h[0], expectedResult_allDouble);

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng(0, curNumVecs);
      auto X_sub = Kokkos::subview(X, Kokkos::ALL(), vecRng);
      auto Y_sub = Kokkos::subview(Y, Kokkos::ALL(), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy(Y, ZERO_d);
      KokkosSparse::spmv("N", alpha_d, A_d, X, beta_d, Y);

      Kokkos::deep_copy(Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0, j);
        EXPECT_EQ(actualResult, expectedResult_allDouble);
      }
    }
  }

  // Now exercise the case where the matrix is in float, but the
  // vectors are in double.
  {
    constexpr float ZERO_f  = static_cast<float>(0.0);
    constexpr float ONE_f   = static_cast<float>(1.0);
    constexpr float TWO_f   = static_cast<float>(2.0);
    constexpr double ZERO_d = static_cast<double>(0.0);

    float_matrix_type A_f("A_f", G);
    auto A_f_val_h = Kokkos::create_mirror_view(A_f.values);
    A_f_val_h[0]   = ONE_f;
    A_f_val_h[1]   = EPS_f / TWO_f;
    EXPECT_NE(A_f_val_h[1], ZERO_f);  // just making sure
    Kokkos::deep_copy(A_f.values, A_f_val_h);

    // Just to make sure, we purge the previous contents of y,
    // before doing the sparse mat-vec.
    Kokkos::deep_copy(y, ZERO_d);
    KokkosSparse::spmv("N", alpha_d, A_f, x, beta_d, y);

    Kokkos::deep_copy(y_h, y);
    const double expectedResult_mixed =
        static_cast<double>(1.0) +
        static_cast<double>(EPS_f) / static_cast<double>(2.0);
    EXPECT_NE(expectedResult_mixed, ZERO_d);
    EXPECT_EQ(y_h[0], expectedResult_mixed);

    for (int curNumVecs = 1; curNumVecs <= numVecs; ++curNumVecs) {
      const Kokkos::pair<int, int> vecRng(0, curNumVecs);
      auto X_sub = Kokkos::subview(X, Kokkos::ALL(), vecRng);
      auto Y_sub = Kokkos::subview(Y, Kokkos::ALL(), vecRng);

      // Just to make sure, we purge the previous contents of Y,
      // before doing the sparse mat-vec.
      Kokkos::deep_copy(Y, ZERO_d);
      KokkosSparse::spmv("N", alpha_d, A_f, X, beta_d, Y);

      Kokkos::deep_copy(Y_h, Y);
      for (int j = 0; j < curNumVecs; ++j) {
        const double actualResult = Y_h(0, j);
        EXPECT_EQ(actualResult, expectedResult_mixed);
      }
    }
  }
}

#define EXECUTE_TEST_ISSUE_101(DEVICE)                                    \
  TEST_F(TestCategory, sparse##_##spmv_issue_101##_##OFFSET##_##DEVICE) { \
    test_github_issue_101<DEVICE>();                                      \
  }

template <typename CrsMat>
CrsMat make_block_matrix(typename CrsMat::ordinal_type &numRows,
                         typename CrsMat::ordinal_type &numCols,
                         typename CrsMat::ordinal_type &blockSize) {
#if 0
    typedef typename CrsMat::StaticCrsGraphType::row_map_type::non_const_type ptr_type ;
    typedef typename CrsMat::StaticCrsGraphType::entries_type::non_const_type ind_type ;
    typedef typename CrsMat::values_type::non_const_type val_type ;
    typedef typename CrsMat::size_type size_type;
#endif
  typedef typename CrsMat::ordinal_type lno_t;
  typedef typename CrsMat::value_type scalar_t;

  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;

  Kokkos::Random_XorShift64<Kokkos::HostSpace> rand(13718);

  // fill outputs with random values
  // Kokkos::Random_XorShift64_Pool<Kokkos::HostSpace> rand_pool(13718);
  // Kokkos::fill_random(hi_x, rand_pool, randomUpperBound<typename
  // hi_scalar_view_t::value_type>(10));

  std::vector<scalar_t> values;
  std::vector<lno_t> rowmap;
  std::vector<lno_t> entries;

  // each row of blocks
  for (lno_t bi = 0; bi < numRows; bi += blockSize) {
    // target number of blocks in the row
    lno_t rowBlockCount = 3;
    {
      // cap the number of blocks in the row
      lno_t maxBlocksInRow = numCols / blockSize;
      rowBlockCount        = std::min(maxBlocksInRow, rowBlockCount);
    }

    // where the blocks in this row of blocks start
    // add that many blocks at random positions in the row
    std::vector<lno_t> bjs;
    for (int _ = 0; _ < rowBlockCount; ++_) {
      bjs.push_back(rand.rand(numCols / blockSize) * blockSize);
    }

    // remove duplicates
    {
      std::sort(bjs.begin(), bjs.end());
      auto it = std::unique(bjs.begin(), bjs.end());
      bjs.resize(it - bjs.begin());
    }

    for (lno_t i = bi; i < bi + blockSize; ++i) {
      rowmap.push_back(entries.size());  // where this row starts

      // for each block
      for (size_t block = 0; block < bjs.size(); ++block) {
        lno_t bj = bjs[block];
        for (lno_t j = bj; j < bj + blockSize; ++j) {
          entries.push_back(j);
          values.push_back(rand.rand(10));
          // values.push_back(1);
        }
      }
    }
  }

  while (rowmap.size() < numRows + 1) {
    rowmap.push_back(entries.size());
  }

  return CrsMat("", numRows, numCols, values.size(), values.data(),
                rowmap.data(), entries.data());
}

struct Coordinate {
  int i;
  int j;
  Coordinate(int _i, int _j) : i(_i), j(_j) {}
  // sort by i then j
  static bool by_ij(const Coordinate &a, const Coordinate &b) {
    if (a.i < b.i) {
      return true;
    } else if (a.i > b.i) {
      return false;
    } else {
      return a.j < b.j;
    }
  }
};
struct Entry {
  Coordinate c;
  double e;
  Entry(int i, int j, double _e) : c(i, j), e(_e) {}
  static bool by_ij(const Entry &a, const Entry &b) {
    return Coordinate::by_ij(a.c, b.c);
  }
};

// expand a pattern into a blocked CrsMatrix
template <typename Matrix,
          std::enable_if_t<is_crs_matrix<Matrix>::value, bool> = true>
Matrix expand_matrix(std::vector<Coordinate> pattern, const int m, const int k,
                     const int blockSize, const int seed = 0) {
  typedef typename Matrix::value_type Scalar;
  typedef typename Matrix::ordinal_type Ordinal;
  typedef typename Matrix::non_const_size_type Offset;
  typedef Kokkos::View<const Offset *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedRowmap;
  typedef Kokkos::View<const Ordinal *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedEntries;
  typedef Kokkos::View<const Scalar *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedValues;

  srand(seed);

  auto gen_rand = []() -> double { return rand() % 10; };

  // check rows and columns
  for (const Coordinate &c : pattern) {
    if (c.i >= m) {
      KokkosKernels::Impl::throw_runtime_exception("i exceeded matrix rows");
    }
    if (c.j >= k) {
      KokkosKernels::Impl::throw_runtime_exception("j exceeded matrix cols");
    }
  }

  // order the blocks
  std::sort(pattern.begin(), pattern.end(), Coordinate::by_ij);

  // create coo entries for each block
  std::vector<Entry> entries;
  for (const Coordinate &c : pattern) {
    for (int i = 0; i < blockSize; ++i) {
      for (int j = 0; j < blockSize; ++j) {
        entries.push_back(
            Entry(c.i * blockSize + i, c.j * blockSize + j, gen_rand()));
      }
    }
  }

  std::sort(entries.begin(), entries.end(), Entry::by_ij);

  std::vector<Offset> rowMap;
  std::vector<Ordinal> colInd;
  std::vector<Scalar> val;

  for (Entry &e : entries) {
    while (rowMap.size() < size_t(e.c.i + 1)) {  // catch empty rows
      rowMap.push_back(colInd.size());
    }
    colInd.push_back(e.c.j);
    val.push_back(e.e);
  }
  // possibly empty rows at end of matrix
  while (rowMap.size() <= size_t(m * blockSize)) {
    rowMap.push_back(colInd.size());
  }

  typename Matrix::row_map_type::non_const_type sparseRowMap("", rowMap.size());
  Kokkos::deep_copy(sparseRowMap,
                    UnmanagedRowmap(rowMap.data(), rowMap.size()));
  typename Matrix::index_type::non_const_type sparseCols("", colInd.size());
  Kokkos::deep_copy(sparseCols, UnmanagedEntries(colInd.data(), colInd.size()));
  typename Matrix::values_type::non_const_type sparseVals("", val.size());
  Kokkos::deep_copy(sparseVals, UnmanagedValues(val.data(), val.size()));

  Matrix mat("crs", m * blockSize, k * blockSize, sparseVals.size(), sparseVals,
             sparseRowMap, sparseCols);
  return mat;
}

template <
    typename Matrix,
    std::enable_if_t<KokkosSparse::Experimental::is_bsr_matrix<Matrix>::value,
                     bool> = true>
Matrix expand_matrix(std::vector<Coordinate> pattern, const int m, const int k,
                     const int blockSize, const int seed = 0) {
  typedef typename Matrix::value_type Scalar;
  typedef typename Matrix::ordinal_type Ordinal;
  typedef typename Matrix::non_const_size_type Offset;
  typedef Kokkos::View<const Offset *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedRowmap;
  typedef Kokkos::View<const Ordinal *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedEntries;
  typedef Kokkos::View<const Scalar *, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      UnmanagedValues;

  srand(seed);

  auto gen_rand = []() -> double { return rand() % 10; };

  // determine the number of rows and columns
  // check rows and columns
  for (const Coordinate &c : pattern) {
    if (c.i >= m) {
      KokkosKernels::Impl::throw_runtime_exception("i exceeded matrix rows");
    }
    if (c.j >= k) {
      KokkosKernels::Impl::throw_runtime_exception("j exceeded matrix cols");
    }
  }

  // order the blocks
  std::sort(pattern.begin(), pattern.end(), Coordinate::by_ij);

  // create values in order of the blocks (storage order for BSR)
  std::vector<Scalar> val(pattern.size() * blockSize * blockSize);
  for (typename std::vector<Scalar>::size_type idx = 0; idx < val.size();
       ++idx) {
    val[idx] = gen_rand();
  }

  /* create the BsrMatrix adjacency info
     use the sorted pattern. val is already in the correct storage order
  */
  std::vector<Offset> rowMap;
  std::vector<Ordinal> colInd;

  for (Coordinate &e : pattern) {
    while (rowMap.size() < size_t(e.i + 1)) {  // catch empty rows
      rowMap.push_back(colInd.size());
    }
    colInd.push_back(e.j);
  }
  // possibly empty rows at end of matrix
  while (rowMap.size() <= size_t(m)) {
    rowMap.push_back(colInd.size());
  }

  typename Matrix::row_map_type::non_const_type sparseRowMap("", rowMap.size());
  Kokkos::deep_copy(sparseRowMap,
                    UnmanagedRowmap(rowMap.data(), rowMap.size()));
  typename Matrix::index_type::non_const_type sparseCols("", colInd.size());
  Kokkos::deep_copy(sparseCols, UnmanagedEntries(colInd.data(), colInd.size()));
  typename Matrix::values_type::non_const_type sparseVals("", val.size());
  Kokkos::deep_copy(sparseVals, UnmanagedValues(val.data(), val.size()));
  Kokkos::fence();

  Matrix mat("bsr", m, k, sparseVals.size(), sparseVals, sparseRowMap,
             sparseCols, blockSize);
  return mat;
}

/* a_scalar_t: the matrix type
   x_scalar_t: the x-vector type
   y_scalar_t: the y-vector type

   blockSize: the size of the dense blocks in the matrix
   pattern: the non-zero locations of the blocks
   m,n: the multiplication dimensions (in terms of blockSize)
   k: number of vectors in the multivector
   y[m*blockSize x k] = A[m*blockSize x n*blockSize] * x[n*blockSize x k]

   Compare the BsrMatrix spmv against a KokkosSparse::spmv on the same operands.
   The controls are used in the BsrMatrix SpMV invocation

*/
template <typename a_scalar_t, typename x_scalar_t, typename y_scalar_t,
          typename lno_t, typename size_type, typename Layout, typename Device>
void test_spmv_bsrmatrix_controls_pattern(
    const KokkosKernels::Experimental::Controls &controls,
    const std::vector<Coordinate> &pattern, const int m, const int n,
    lno_t blockSize, lno_t k, y_scalar_t alpha, y_scalar_t beta) {
  // get the widest passed scalar type
  // typedef typename std::conditional<sizeof(a_scalar_t) >= sizeof(x_scalar_t),
  //                                   a_scalar_t, x_scalar_t>::type wider_t;
  // typedef typename std::conditional<sizeof(wider_t) >= sizeof(y_scalar_t),
  //                                   wider_t, y_scalar_t>::type widest_t;

  typedef typename KokkosSparse::CrsMatrix<a_scalar_t, lno_t, Device, void,
                                           size_type>
      crs_mat_t;
  typedef
      typename KokkosSparse::Experimental::BsrMatrix<a_scalar_t, lno_t, Device,
                                                     void, size_type>
          bsr_mat_t;
  typedef Kokkos::View<x_scalar_t **, Layout, Device> x_view_t;
  typedef Kokkos::View<y_scalar_t **, Layout, Device> y_view_t;

  using DeviceRangePolicy = Kokkos::RangePolicy<Device>;

  crs_mat_t crs = expand_matrix<crs_mat_t>(pattern, m, n, blockSize);
  bsr_mat_t bsr = expand_matrix<bsr_mat_t>(pattern, m, n, blockSize);

  // only tue if the original matrix is a multiple of block size, and all blocks
  // are dense
  EXPECT_TRUE(bsr.nnz() * bsr.blockDim() * bsr.blockDim() == crs.nnz());
  EXPECT_TRUE(bsr.numRows() * bsr.blockDim() == crs.numRows());
  EXPECT_TRUE(bsr.numCols() * bsr.blockDim() == crs.numCols());

  // expected operands
  x_view_t exp_x("exp_x", n * blockSize, k);
  y_view_t exp_y("exp_y", m * blockSize, k);

  // test operands
  y_view_t test_y("test_y", m * blockSize, k);
  x_view_t test_x("test_x", n * blockSize, k);

  // fill expected with random values
  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);
  Kokkos::fill_random(exp_x, rand_pool,
                      randomUpperBound<typename x_view_t::value_type>(10));
  Kokkos::fill_random(exp_y, rand_pool,
                      randomUpperBound<typename y_view_t::value_type>(10));

#if 0
  // fill inputs with 1, for help debugging
  Kokkos::parallel_for("fill",
    Kokkos::MDRangePolicy<Device, Kokkos::Rank<2>>({0,0}, {hi_x.extent(0), hi_x.extent(1)}),
    KOKKOS_LAMBDA (unsigned i, unsigned j) { 
        hi_x(i,j) = 1 + (i == 0 && j == 0); 
    }
  );
#endif

  // copy expected operands to test operands
  Kokkos::deep_copy(test_x, exp_x);
  Kokkos::deep_copy(test_y, exp_y);
  Kokkos::fence();

  // generate expected y vector
  // some error about Blas implementation
  KokkosSparse::spmv("N", alpha, crs, exp_x, beta, exp_y);
  Kokkos::fence();

  // invoke tensor-core spmv
  KokkosSparse::spmv(controls, "N", alpha, bsr, test_x, beta, test_y);
  Kokkos::fence();

  // test each vector
  for (lno_t ki = 0; ki < k; ++ki) {
    auto exp_y_i  = Kokkos::subview(exp_y, Kokkos::ALL(), ki);
    auto test_y_i = Kokkos::subview(test_y, Kokkos::ALL(), ki);

    // count errors
    int num_errors = 0;
    // Kokkos::ArithTraits<half> in CUDA 9 is float on the host
    // for CUDA 9, Kokkos half is actually float. However, the tensor core SpMV
    // uses CUDA's half type, not Kokkos, so we still need a reduced precision
    // test.
    double eps =
        KOKKOSKERNELS_IMPL_FP16_EPSILON * KOKKOSKERNELS_IMPL_FP16_RADIX;
    Kokkos::parallel_reduce("KokkosSparse::Test::spmv_tc",
                            DeviceRangePolicy(0, exp_y_i.extent(0)),
                            Test::fSPMV<decltype(exp_y_i), decltype(test_y_i)>(
                                exp_y_i, test_y_i, eps),
                            num_errors);
    // explicit cast to double since no overload for half::operator<<
    if (num_errors > 0)
      std::cout << "KokkosSparse::Test::spmv_tc: " << num_errors
                << " errors of " << exp_y_i.extent_int(0) << " for mv " << ki
                << " (alpha="
                << double(Kokkos::ArithTraits<y_scalar_t>::abs(alpha))
                << ", beta="
                << double(Kokkos::ArithTraits<y_scalar_t>::abs(beta))
                << ", mode = N"
                << ")\n";
    EXPECT_TRUE(num_errors == 0);
  }
}

/* test a particular pattern with all supported controls
 */
template <typename a_scalar_t, typename x_scalar_t, typename y_scalar_t,
          typename lno_t, typename size_type, typename Layout, typename Device>
void test_spmv_bsrmatrix_pattern(const std::vector<Coordinate> &pattern,
                                 const int m, const int n, lno_t blockSize,
                                 lno_t k, y_scalar_t alpha, y_scalar_t beta) {
  {
    KokkosKernels::Experimental::Controls controls;
    controls.setParameter("algorithm", "experimental_bsr_tc");
    test_spmv_bsrmatrix_controls_pattern<a_scalar_t, x_scalar_t, y_scalar_t,
                                         lno_t, size_type, Layout, Device>(
        controls, pattern, m, n, blockSize, k, alpha, beta);
  }

#if defined(KOKKOS_ARCH_AMPERE)
  {
    KokkosKernels::Experimental::Controls controls;
    controls.setParameter("algorithm", "experimental_bsr_tc");
    controls.setParameter("tc_precision", "double");
    test_spmv_bsrmatrix_controls_pattern<a_scalar_t, x_scalar_t, y_scalar_t,
                                         lno_t, size_type, Layout, Device>(
        controls, pattern, m, n, blockSize, k, alpha, beta);
  }
#endif
}

/* test a bunch of different matrices
 */
template <typename a_scalar_t, typename x_scalar_t, typename y_scalar_t,
          typename lno_t, typename size_type, typename Layout, typename Device>
void test_spmv_bsrmatrix(lno_t blockSize, lno_t k, y_scalar_t alpha,
                         y_scalar_t beta) {
  KokkosKernels::Experimental::Controls controls;
  controls.setParameter("algorithm", "experimental_bsr_tc");

  // 1x1 full
  {
    int m                           = 1;
    int n                           = 1;
    std::vector<Coordinate> pattern = {Coordinate(0, 0)};
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 1x1 empty
  {
    int m                           = 1;
    int n                           = 1;
    std::vector<Coordinate> pattern = {};
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 2x2 top-left
  {
    int m                           = 2;
    int n                           = 2;
    std::vector<Coordinate> pattern = {Coordinate(0, 0)};
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 2x2 bottom right
  {
    int m                           = 2;
    int n                           = 2;
    std::vector<Coordinate> pattern = {Coordinate(1, 1)};
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 2x3 bottom right
  {
    int m                           = 2;
    int n                           = 3;
    std::vector<Coordinate> pattern = {Coordinate(1, 2)};
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 2x10 long bottom row
  {
    int m = 2;
    int n = 10;
    std::vector<Coordinate> pattern;
    for (int j = 0; j < n; ++j) {
      pattern.push_back(Coordinate(1, j));
    }
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }

  // 10x10 column 1 + diagonal
  {
    int m = 10;
    int n = 10;
    std::vector<Coordinate> pattern;
    for (int i = 0; i < n; ++i) {
      pattern.push_back(Coordinate(i, 1));
      if (i != 1) {
        pattern.push_back(Coordinate(i, i));
      }
    }
    test_spmv_bsrmatrix_pattern<a_scalar_t, x_scalar_t, y_scalar_t, lno_t,
                                size_type, Layout, Device>(
        pattern, m, n, blockSize, k, alpha, beta);
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                          \
  TEST_F(TestCategory,                                                         \
         sparse##_##spmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {       \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 3, 200, 10, true); \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 3, 100, 10, true); \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 20, 100, 5, true); \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 3, 20, 10,       \
                                               false);                         \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 3, 100, 10,      \
                                               false);                         \
    test_spmv<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 10000 * 2, 100, 5,       \
                                               false);                         \
    test_spmv_controls<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 10000 * 20,     \
                                                        100, 5);               \
  }

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                    \
  TEST_F(                                                                           \
      TestCategory,                                                                 \
      sparse##_##spmv_mv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        1000, 1000 * 3, 200, 10, true, 1);                                          \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        1000, 1000 * 3, 100, 10, true, 5);                                          \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        1000, 1000 * 2, 100, 5, true, 10);                                          \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        50000, 50000 * 3, 20, 10, false, 1);                                        \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        50000, 50000 * 3, 100, 10, false, 1);                                       \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(                  \
        10000, 10000 * 2, 100, 5, false, 5);                                        \
    test_spmv_mv_heavy<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(            \
        200, 200 * 10, 60, 4, 30);                                                  \
  }

#define EXECUTE_TEST_STRUCT(SCALAR, ORDINAL, OFFSET, DEVICE)                   \
  TEST_F(                                                                      \
      TestCategory,                                                            \
      sparse##_##spmv_struct##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {   \
    test_spmv_struct_1D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 1, 1);            \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(25, 21, 3, 3);        \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(20, 25, 3, 3);        \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(22, 22, 3, 3);        \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(20, 20, 20, 3, 3, 3); \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(22, 22, 22, 3, 3, 3); \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(25, 10, 20, 3, 3, 3); \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 20, 25, 3, 3, 3); \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 24, 20, 3, 3, 3); \
  }

#define EXECUTE_TEST_MV_STRUCT(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                    \
  TEST_F(                                                                                  \
      TestCategory,                                                                        \
      sparse##_##spmv_mv_struct##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_spmv_mv_struct_1D<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(               \
        10, 1);                                                                            \
    test_spmv_mv_struct_1D<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(               \
        10, 2);                                                                            \
  }

/* Tensor Core SpMV
  blocksize, k, alpha, beta
*/
#define EXECUTE_TEST_TC(ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET, LAYOUT,                                           \
                        DEVICE)                                                                                       \
  TEST_F(                                                                                                             \
      TestCategory,                                                                                                   \
      sparse##_##spmv_tensor_core##_##ASCALAR##_##XSCALAR##_##YSCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    /* easy case with different alphas and betas*/                                                                    \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(16, 16, 0, 0);                                                        \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(16, 16, 1, 0);                                                        \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(16, 16, 0, 1);                                                        \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(16, 16, 1, 1);                                                        \
    /* easy case with a real alpha/beta */                                                                            \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(16, 16, 1.25, -2.73);                                                 \
    /* smaller block size with k < and > block size*/                                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(7, 6, 1.25, -2.73);                                                   \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(7, 7, 1.25, -2.73);                                                   \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(7, 8, 1.25, -2.73);                                                   \
    /* smaller block size with k < and > block size*/                                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(15, 14, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(15, 15, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(15, 16, 1.25, -2.73);                                                 \
    /* larger block size with k < and > block size*/                                                                  \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(17, 16, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(17, 17, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(17, 18, 1.25, -2.73);                                                 \
    /* larger block size with k < and > block size*/                                                                  \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(32, 31, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(32, 32, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(32, 33, 1.25, -2.73);                                                 \
    /* more than one team per block*/                                                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(33, 13, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(33, 27, 1.25, -2.73);                                                 \
    test_spmv_bsrmatrix<ASCALAR, XSCALAR, YSCALAR, ORDINAL, OFFSET,                                                   \
                        Kokkos::LAYOUT, DEVICE>(33, 41, 1.25, -2.73);                                                 \
  }

// minimal conditions for tensor core SpMV test
// BsrMatrix spmv is only supported on CUDA for the time being
#if defined(KOKKOS_ENABLE_CUDA) && defined(TEST_CUDA_SPARSE_CPP) && \
    (defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_AMPERE))

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&       \
        defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T) && \
        defined(KOKKOSKERNELS_INST_FLOAT) &&         \
        defined(KOKKOSKERNELS_INST_LAYOUTLEFT) ||    \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&             \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
// EXECUTE_TEST_TC(kokkos_half,  kokkos_half, float,   int, size_t, LayoutLeft,
// TestExecSpace) EXECUTE_TEST_TC(kokkos_half,  float,       float,   int,
// size_t, LayoutLeft,  TestExecSpace) EXECUTE_TEST_TC(float, kokkos_half,
// float,   int, size_t, LayoutLeft,  TestExecSpace)
EXECUTE_TEST_TC(float, float, float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&       \
        defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T) && \
        defined(KOKKOSKERNELS_INST_DOUBLE) &&        \
        defined(KOKKOSKERNELS_INST_LAYOUTLEFT) ||    \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&             \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
// EXECUTE_TEST_TC(kokkos_half,  kokkos_half, double,   int, size_t, LayoutLeft,
// TestExecSpace) EXECUTE_TEST_TC(kokkos_half,  double,       double,   int,
// size_t, LayoutLeft,  TestExecSpace) EXECUTE_TEST_TC(double, kokkos_half,
// double,   int, size_t, LayoutLeft,  TestExecSpace)
EXECUTE_TEST_TC(double, double, double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&       \
        defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T) && \
        defined(KOKKOSKERNELS_INST_FLOAT) &&         \
        defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&             \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
// EXECUTE_TEST_TC(kokkos_half,  kokkos_half, float,   int, size_t, LayoutRight,
// TestExecSpace) EXECUTE_TEST_TC(kokkos_half,  float,       float,   int,
// size_t, LayoutRight,  TestExecSpace) EXECUTE_TEST_TC(float, kokkos_half,
// float,   int, size_t, LayoutRight,  TestExecSpace)
EXECUTE_TEST_TC(float, float, float, int, size_t, LayoutRight, TestExecSpace)
#endif

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&       \
        defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T) && \
        defined(KOKKOSKERNELS_INST_DOUBLE) &&        \
        defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&             \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
// EXECUTE_TEST_TC(kokkos_half,  kokkos_half, double,   int, size_t,
// LayoutRight,  TestExecSpace) EXECUTE_TEST_TC(kokkos_half,  double, double,
// int, size_t, LayoutRight,  TestExecSpace) EXECUTE_TEST_TC(double,
// kokkos_half, double,   int, size_t, LayoutRight,  TestExecSpace)
EXECUTE_TEST_TC(double, double, double, int, size_t, LayoutRight, TestExecSpace)
#endif

#endif  // tensor core SpMV tests

#undef EXECUTE_TEST_TC

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_ISSUE_101(TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestExecSpace)
EXECUTE_TEST_STRUCT(double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestExecSpace)
EXECUTE_TEST_STRUCT(double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&         \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, int, TestExecSpace)
EXECUTE_TEST_STRUCT(float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, int, TestExecSpace)
EXECUTE_TEST_STRUCT(float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(float, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
EXECUTE_TEST_STRUCT(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&  \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(double, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&      \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int64_t, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(double, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&         \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(double, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&      \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int64_t, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(double, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&  \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(float, int, int, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&      \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int64_t, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(float, int64_t, int, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(float, int, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&      \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int64_t, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(float, int64_t, size_t, LayoutLeft, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&             \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int, int, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&             \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int64_t, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int64_t, int, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&             \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int, size_t, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&             \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int64_t, size_t, LayoutLeft,
                TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_double, int64_t, size_t, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int, int, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int64_t, int, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int64_t, int, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int, size_t, LayoutLeft, TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int, size_t, LayoutLeft,
                       TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_LAYOUTLEFT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int64_t, size_t, LayoutLeft,
                TestExecSpace)
EXECUTE_TEST_MV_STRUCT(kokkos_complex_float, int64_t, size_t, LayoutLeft,
                       TestExecSpace)
#endif
#endif  // defined(KOKKOSKERNELS_INST_LAYOUTLEFT)

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&         \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(double, int64_t, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||     \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&    \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&    \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&           \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&     \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||  \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&            \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(float, int64_t, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||            \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||         \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                   \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_double, int64_t, size_t, LayoutRight,
                TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT)) ||           \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int64_t, int, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int, size_t, LayoutRight, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) ||        \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&                  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_MV(kokkos_complex_float, int64_t, size_t, LayoutRight,
                TestExecSpace)
#endif

#undef EXECUTE_TEST
#undef EXECUTE_TEST_STRUCT
#undef EXECUTE_TEST_MV
#undef EXECUTE_TEST_MV_STRUCT
