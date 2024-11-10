//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_Test_Structured_Matrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include <KokkosKernels_Utils.hpp>
#include <KokkosKernels_NaN.hpp>

#include "KokkosKernels_default_types.hpp"

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;
typedef Kokkos::Experimental::half_t kokkos_half;

namespace Test {

// Functor checking that the results of SPMV
// are consistent with a reference sequential
// implementation of the same operation.
//
// Inputs:
// - _ex_y      the expected result calculated
//              from the reference implementation
// - _y         the result from optimized SPMV being
//              tested for correctness
// - _eps       the tolerance required to accept the
//              results as correct
// - _max_val   the largest possible value that can
//              be stored as an intermediate result
//              during the computation
//
//  The criteria to assess correctness is
//     abs(_ex_y - _y) / _max_val < tol
//
//  Note: _max_val in the case of SPMV can be computed
//  as follows. Find the max number of entries per
//  row in the matrix (max_row_length), also find the
//  largest value that can be stored in the matrix, x
//  and y vectors (max_mat, max_x and max_y).
//
//     _max_val = beta*max_y
//                + alpha*max_row_length*max_mat*max_x
template <class VectorType0, class VectorType1>
struct fSPMV {
  using value_type = int;
  using AT         = Kokkos::ArithTraits<typename VectorType1::non_const_value_type>;
  using ATM        = Kokkos::ArithTraits<typename AT::mag_type>;
  using mag_type   = typename AT::mag_type;

  VectorType0 expected_y;
  VectorType1 y;
  mag_type eps;
  mag_type max_val;

  fSPMV(const VectorType0 &_ex_y, const VectorType1 &_y, const mag_type _eps, const mag_type _max_val = ATM::one())
      : expected_y(_ex_y), y(_y), eps(AT::abs(_eps)), max_val(AT::abs(_max_val)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &err) const {
    const mag_type error = AT::abs(expected_y(i) - y(i));

    // only one is NaN or error is too large
    if ((Kokkos::isnan(AT::abs(expected_y(i))) ^ Kokkos::isnan(AT::abs(y(i)))) || (error > eps * max_val)) {
      err++;
      Kokkos::printf("expected_y(%d)=%f, y(%d)=%f err=%e, max_error=%e\n", i, AT::abs(expected_y(i)), i, AT::abs(y(i)),
                     error, eps * max_val);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, value_type &err) const {
    const mag_type error = AT::abs(expected_y(i, j) - y(i, j));

    if (error > eps * max_val) {
      err++;
      Kokkos::printf("expected_y(%d,%d)=%f, y(%d,%d)=%f err=%e, max_error=%e\n", i, j, AT::abs(expected_y(i, j)), i, j,
                     AT::abs(y(i, j)), error, eps * max_val);
    }
  }
};

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void sequential_spmv(crsMat_t input_mat, x_vector_type x, y_vector_type y,
                     typename y_vector_type::non_const_value_type alpha,
                     typename y_vector_type::non_const_value_type beta, const std::string &mode = "N") {
  using graph_t          = typename crsMat_t::StaticCrsGraphType;
  using size_type_view_t = typename graph_t::row_map_type;
  using lno_view_t       = typename graph_t::entries_type;
  using scalar_view_t    = typename crsMat_t::values_type::non_const_type;
  using y_scalar_t       = typename y_vector_type::non_const_value_type;

  using size_type = typename size_type_view_t::non_const_value_type;
  using lno_t     = typename lno_view_t::non_const_value_type;
  using scalar_t  = typename scalar_view_t::non_const_value_type;
  using KAT       = Kokkos::ArithTraits<scalar_t>;

  typename scalar_view_t::HostMirror h_values = Kokkos::create_mirror_view(input_mat.values);
  Kokkos::deep_copy(h_values, input_mat.values);

  typename lno_view_t::HostMirror h_entries = Kokkos::create_mirror_view(input_mat.graph.entries);
  Kokkos::deep_copy(h_entries, input_mat.graph.entries);

  typename size_type_view_t::HostMirror h_rowmap = Kokkos::create_mirror_view(input_mat.graph.row_map);
  Kokkos::deep_copy(h_rowmap, input_mat.graph.row_map);
  Kokkos::fence();

  typename y_vector_type::HostMirror h_y = Kokkos::create_mirror_view(y);
  typename x_vector_type::HostMirror h_x = Kokkos::create_mirror_view(x);

  KokkosKernels::Impl::safe_device_to_host_deep_copy(x.extent(0), x, h_x);
  KokkosKernels::Impl::safe_device_to_host_deep_copy(y.extent(0), y, h_y);
  Kokkos::fence();

  lno_t nr = input_mat.numRows();

  // first, scale y by beta
  for (size_t i = 0; i < h_y.extent(0); i++) {
    if (beta == y_scalar_t(0)) {
      h_y(i) = y_scalar_t(0);
    } else {
      h_y(i) *= beta;
    }
  }

  // then go through the matrix and accumulate the matrix-vector product
  for (lno_t row = 0; row < nr; ++row) {
    for (size_type j = h_rowmap(row); j < h_rowmap(row + 1); ++j) {
      lno_t col    = h_entries(j);
      scalar_t val = h_values(j);
      if (mode == "N")
        h_y(row) += alpha * val * h_x(col);
      else if (mode == "C")
        h_y(row) += alpha * KAT::conj(val) * h_x(col);
      else if (mode == "T")
        h_y(col) += alpha * val * h_x(row);
      else if (mode == "H")
        h_y(col) += alpha * KAT::conj(val) * h_x(row);
    }
  }
  KokkosKernels::Impl::safe_host_to_device_deep_copy(y.extent(0), h_y, y);
  Kokkos::fence();
}

template <typename handle_t, typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv(handle_t *handle, crsMat_t input_mat, x_vector_type x, y_vector_type y,
                typename y_vector_type::non_const_value_type alpha, typename y_vector_type::non_const_value_type beta,
                const std::string &mode,
                typename Kokkos::ArithTraits<typename crsMat_t::value_type>::mag_type max_val) {
  EXPECT_TRUE(mode.size() == 1);

  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  const y_value_mag_type eps = 10 * Kokkos::ArithTraits<y_value_mag_type>::eps();

  y_vector_type actual_y("actual_y", y.extent(0));
  y_vector_type expected_y("expected_y", y.extent(0));
  Kokkos::deep_copy(expected_y, y);
  Kokkos::deep_copy(actual_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta, mode);
  bool threw = false;
  std::string msg;
  try {
    KokkosSparse::spmv(handle, mode.data(), alpha, input_mat, x, beta, actual_y);
    Kokkos::fence();
  } catch (std::exception &e) {
    threw = true;
    msg   = e.what();
  }
  ASSERT_FALSE(threw) << "KokkosSparse::Test::spmv 1D, mode " << mode << ": threw exception:\n" << msg << '\n';

  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosSparse::Test::spmv", my_exec_space(0, actual_y.extent(0)),
                          fSPMV(expected_y, actual_y, eps, max_val), num_errors);
  if (num_errors > 0)
    printf("KokkosSparse::Test::spmv: %i errors of %i with params: %lf %lf\n", num_errors, y.extent_int(0),
           y_value_trait::abs(alpha), y_value_trait::abs(beta));
  EXPECT_TRUE(num_errors == 0);
}

template <typename Handle, typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv(Handle *handle, crsMat_t input_mat, x_vector_type x, y_vector_type y, y_vector_type expected_y,
                   typename y_vector_type::non_const_value_type alpha,
                   typename y_vector_type::non_const_value_type beta, int numMV, const std::string &mode,
                   typename Kokkos::ArithTraits<typename crsMat_t::value_type>::mag_type max_val) {
  EXPECT_TRUE(mode.size() == 1);

  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const y_value_mag_type eps = 10 * Kokkos::ArithTraits<y_value_mag_type>::eps();

  Kokkos::deep_copy(expected_y, y);

  Kokkos::fence();

  bool threw = false;
  std::string msg;
  try {
    KokkosSparse::spmv(handle, mode.data(), alpha, input_mat, x, beta, y);
    Kokkos::fence();
  } catch (std::exception &e) {
    threw = true;
    msg   = e.what();
  }
  ASSERT_FALSE(threw) << "KokkosSparse::Test::spmv 2D, mode " << mode << ": threw exception:\n" << msg << '\n';

  for (int i = 0; i < numMV; ++i) {
    auto x_i = Kokkos::subview(x, Kokkos::ALL(), i);

    auto y_i = Kokkos::subview(expected_y, Kokkos::ALL(), i);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta, mode);

    auto y_spmv    = Kokkos::subview(y, Kokkos::ALL(), i);
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosSparse::Test::spmv_mv", my_exec_space(0, y_i.extent(0)),
                            fSPMV(y_i, y_spmv, eps, max_val), num_errors);
    if (num_errors > 0)
      std::cout << "KokkosSparse::Test::spmv_mv: " << num_errors << " errors of " << y_i.extent_int(0) << " for mv "
                << i << " (alpha=" << alpha << ", beta=" << beta << ", mode = " << mode << ")\n";
    EXPECT_TRUE(num_errors == 0);
  }
}

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_struct(const crsMat_t input_mat, const int stencil_type,
                       const Kokkos::View<typename crsMat_t::non_const_ordinal_type *, Kokkos::HostSpace> structure,
                       x_vector_type x, y_vector_type y, typename y_vector_type::non_const_value_type alpha,
                       typename y_vector_type::non_const_value_type beta,
                       typename Kokkos::ArithTraits<typename crsMat_t::value_type>::mag_type max_val) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const double eps = Kokkos::ArithTraits<y_value_mag_type>::eps();
  const size_t nr  = input_mat.numRows();
  y_vector_type expected_y("expected", nr);
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  sequential_spmv(input_mat, x, expected_y, alpha, beta);
  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, alpha, input_mat, x, beta, y);

  int num_errors = 0;
  Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv_struct", my_exec_space(0, y.extent(0)),
                          fSPMV<y_vector_type, y_vector_type>(expected_y, y, eps, max_val), num_errors);
  if (num_errors > 0) {
    printf(
        "KokkosKernels::UnitTests::spmv_struct: %i errors of %i with params: "
        "%d %lf %lf\n",
        num_errors, y.extent_int(0), stencil_type, y_value_trait::abs(alpha), y_value_trait::abs(beta));
  }
  EXPECT_TRUE(num_errors == 0);
}  // check_spmv_struct

template <typename crsMat_t, typename x_vector_type, typename y_vector_type>
void check_spmv_mv_struct(const crsMat_t input_mat, const int stencil_type,
                          const Kokkos::View<typename crsMat_t::non_const_ordinal_type *, Kokkos::HostSpace> structure,
                          x_vector_type x, y_vector_type y, y_vector_type expected_y,
                          typename y_vector_type::non_const_value_type alpha,
                          typename y_vector_type::non_const_value_type beta, int numMV,
                          typename Kokkos::ArithTraits<typename crsMat_t::value_type>::mag_type max_val) {
  using ExecSpace        = typename crsMat_t::execution_space;
  using my_exec_space    = Kokkos::RangePolicy<ExecSpace>;
  using y_value_type     = typename y_vector_type::non_const_value_type;
  using y_value_trait    = Kokkos::ArithTraits<y_value_type>;
  using y_value_mag_type = typename y_value_trait::mag_type;

  // y is the quantity being tested here,
  // so let us use y_value_type to determine
  // the appropriate tolerance precision.
  const double eps = Kokkos::ArithTraits<y_value_mag_type>::eps();
  Kokkos::deep_copy(expected_y, y);
  Kokkos::fence();

  KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, alpha, input_mat, x, beta, y);

  for (int vectorIdx = 0; vectorIdx < numMV; ++vectorIdx) {
    auto x_i = Kokkos::subview(x, Kokkos::ALL(), vectorIdx);
    auto y_i = Kokkos::subview(expected_y, Kokkos::ALL(), vectorIdx);
    Kokkos::fence();

    sequential_spmv(input_mat, x_i, y_i, alpha, beta);

    auto y_spmv    = Kokkos::subview(y, Kokkos::ALL(), vectorIdx);
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosKernels::UnitTests::spmv_mv_struct", my_exec_space(0, y.extent(0)),
                            fSPMV<decltype(y_i), decltype(y_spmv)>(y_i, y_spmv, eps, max_val), num_errors);
    if (num_errors > 0)
      printf(
          "KokkosKernels::UnitTests::spmv_mv_struct: %i errors of %i with "
          "params: %d %lf %lf, in vector %i\n",
          num_errors, y.extent_int(0), stencil_type, y_value_trait::abs(alpha), y_value_trait::abs(beta), vectorIdx);
    EXPECT_TRUE(num_errors == 0);
  }
}  // check_spmv_mv_struct

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

template <typename scalar_t, typename lno_t, typename size_type, typename Device>
void test_spmv(KokkosSparse::SPMVAlgorithm algo, lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance,
               bool heavy) {
  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using mag_t         = typename Kokkos::ArithTraits<scalar_t>::mag_type;
  using handle_t      = KokkosSparse::SPMVHandle<Device, crsMat_t, x_vector_type, y_vector_type>;
  using y_policy      = Kokkos::RangePolicy<typename y_vector_type::execution_space>;

  constexpr mag_t max_x   = static_cast<mag_t>(1);
  constexpr mag_t max_y   = static_cast<mag_t>(1);
  constexpr mag_t max_val = static_cast<mag_t>(1);

  lno_t numCols = numRows;

  crsMat_t input_mat =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numCols, nnz, row_size_variance, bandwidth);
  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  const lno_t max_nnz_per_row = numRows ? (nnz / numRows + row_size_variance) : 0;

  // Create vectors with and without nans
  x_vector_type input_x("x", nc);
  x_vector_type input_xt("x", nr);
  y_vector_type input_y("y", nr), input_y_nans("y_nans", nr);
  y_vector_type input_yt("y", nc), input_yt_nans("y_nans", nc);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(input_x, rand_pool, randomUpperBound<scalar_t>(max_x));
  Kokkos::fill_random(input_y, rand_pool, randomUpperBound<scalar_t>(max_y));
  Kokkos::fill_random(input_xt, rand_pool, randomUpperBound<scalar_t>(max_x));
  Kokkos::fill_random(input_yt, rand_pool, randomUpperBound<scalar_t>(max_y));

  // sprinkle in some nans
  Kokkos::deep_copy(input_y_nans, input_y);
  Kokkos::deep_copy(input_yt_nans, input_yt);
  Kokkos::parallel_for(
      y_policy(0, input_y_nans.extent(0)), KOKKOS_LAMBDA(const size_t i) {
        if (0 == (i % 19)) {
          input_y_nans(i) = KokkosKernels::Impl::quiet_NaN<scalar_t>();
        }
      });
  Kokkos::parallel_for(
      y_policy(0, input_yt_nans.extent(0)), KOKKOS_LAMBDA(const size_t i) {
        if (0 == (i % 23)) {
          input_yt_nans(i) = KokkosKernels::Impl::quiet_NaN<scalar_t>();
        }
      });

  // We also need to bound the values
  // in the matrix to bound the cancellations
  // coming from arithmetic operations.
  Kokkos::fill_random(input_mat.values, rand_pool, randomUpperBound<scalar_t>(max_val));

  std::vector<const char *> nonTransModes = {"N"};
  std::vector<const char *> transModes    = {"T"};
  std::vector<double> testAlphaBeta       = {0.0, 1.0};
  if (heavy) {
    nonTransModes.push_back("C");
    transModes.push_back("H");
    testAlphaBeta.push_back(-1.0);
    testAlphaBeta.push_back(2.5);
  }

  // This handle can be reused for all following calls, since the matrix does
  // not change
  handle_t handle(algo);

  for (auto mode : nonTransModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        mag_t max_error = beta * max_y + alpha * max_nnz_per_row * max_val * max_x;
        Test::check_spmv(&handle, input_mat, input_x, input_y, alpha, beta, mode, max_error);
        if (0 == beta) {
          Test::check_spmv(&handle, input_mat, input_x, input_y_nans, alpha, beta, mode, max_error);
        }
      }
    }
  }
  for (auto mode : transModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        // hoping the transpose won't have a long column...
        mag_t max_error = beta * max_y + alpha * max_nnz_per_row * max_val * max_x;
        Test::check_spmv(&handle, input_mat, input_xt, input_yt, alpha, beta, mode, max_error);
        if (0 == beta) {
          Test::check_spmv(&handle, input_mat, input_x, input_yt_nans, alpha, beta, mode, max_error);
        }
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename Device>
void test_spmv_algorithms(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, bool heavy) {
  using namespace KokkosSparse;
  // Here, SPMV_MERGE_PATH will test a TPL's algorithm for imbalanced matrices
  // if available (like cuSPARSE ALG2). SPMV_NATIVE_MERGE_PATH will always call
  // the KokkosKernels implmentation of merge path.
  for (SPMVAlgorithm algo : {SPMV_DEFAULT, SPMV_NATIVE, SPMV_MERGE_PATH, SPMV_NATIVE_MERGE_PATH}) {
    test_spmv<scalar_t, lno_t, size_type, Device>(algo, numRows, nnz, bandwidth, row_size_variance, heavy);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, class Device>
void test_spmv_mv(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, bool heavy, int numMV) {
  using mag_t = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  constexpr mag_t max_x   = static_cast<mag_t>(1);
  constexpr mag_t max_y   = static_cast<mag_t>(1);
  constexpr mag_t max_val = static_cast<mag_t>(1);

  lno_t numCols = numRows;

  using crsMat_t  = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using ViewTypeX = Kokkos::View<scalar_t **, layout, Device>;
  using ViewTypeY = Kokkos::View<scalar_t **, layout, Device>;
  using handle_t  = KokkosSparse::SPMVHandle<Device, crsMat_t, ViewTypeX, ViewTypeY>;

  ViewTypeX b_x("A", numCols, numMV);
  ViewTypeY b_y("B", numRows, numMV);
  ViewTypeY b_y_copy("B", numRows, numMV);

  ViewTypeX b_xt("A", numRows, numMV);
  ViewTypeY b_yt("B", numCols, numMV);
  ViewTypeY b_yt_copy("B", numCols, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);
  Kokkos::fill_random(b_x, rand_pool, randomUpperBound<scalar_t>(max_x));
  Kokkos::fill_random(b_y, rand_pool, randomUpperBound<scalar_t>(max_y));
  Kokkos::fill_random(b_xt, rand_pool, randomUpperBound<scalar_t>(max_x));
  Kokkos::fill_random(b_yt, rand_pool, randomUpperBound<scalar_t>(max_y));

  crsMat_t input_mat =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numCols, nnz, row_size_variance, bandwidth);

  const lno_t max_nnz_per_row = numRows ? (nnz / numRows + row_size_variance) : 0;

  // We also need to bound the values
  // in the matrix to bound the cancellations
  // coming from arithmetic operations.
  Kokkos::fill_random(input_mat.values, rand_pool, randomUpperBound<scalar_t>(max_val));

  Kokkos::deep_copy(b_y_copy, b_y);
  Kokkos::deep_copy(b_yt_copy, b_yt);

  std::vector<const char *> nonTransModes = {"N"};
  std::vector<const char *> transModes    = {"T"};
  std::vector<double> testAlphaBeta       = {0.0, 1.0};
  if (heavy) {
    nonTransModes.push_back("C");
    transModes.push_back("H");
    testAlphaBeta.push_back(-1.0);
    testAlphaBeta.push_back(2.5);
  }
  handle_t handle;
  for (auto mode : nonTransModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        mag_t max_error = beta * max_y + alpha * max_nnz_per_row * max_val * max_x;
        Test::check_spmv_mv(&handle, input_mat, b_x, b_y, b_y_copy, alpha, beta, numMV, mode, max_error);
      }
    }
  }
  for (auto mode : transModes) {
    for (double alpha : testAlphaBeta) {
      for (double beta : testAlphaBeta) {
        // hoping the transpose won't have a long column...
        mag_t max_error = beta * max_y + alpha * max_nnz_per_row * max_val * max_x;
        Test::check_spmv_mv(&handle, input_mat, b_xt, b_yt, b_yt_copy, alpha, beta, numMV, mode, max_error);
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout_x, typename layout_y, class Device>
void test_spmv_mv_heavy(lno_t numRows, lno_t numCols, size_type nnz, lno_t bandwidth, lno_t row_size_variance,
                        int numMV) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) || defined(KOKKOS_ARCH_A64FX)
  if (std::is_same<scalar_t, Kokkos::complex<double>>::value) {
    std::cerr << "TEST SKIPPED: See "
                 "https://github.com/kokkos/kokkos-kernels/issues/1331 for details."
              << std::endl;
    return;
  }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL || KOKKOS_ARCH_A64FX
  using crsMat_t  = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using ViewTypeX = Kokkos::View<scalar_t **, layout_x, Device>;
  using ViewTypeY = Kokkos::View<scalar_t **, layout_y, Device>;
  using mag_t     = typename Kokkos::ArithTraits<scalar_t>::mag_type;
  using handle_t  = KokkosSparse::SPMVHandle<Device, crsMat_t, ViewTypeX, ViewTypeY>;

  constexpr mag_t max_x   = static_cast<mag_t>(10);
  constexpr mag_t max_y   = static_cast<mag_t>(10);
  constexpr mag_t max_val = static_cast<mag_t>(10);

  crsMat_t input_mat =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numCols, nnz, row_size_variance, bandwidth);
  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  const lno_t max_nnz_per_row = numRows ? (nnz / numRows + row_size_variance) : 0;

  for (int nv = 1; nv <= numMV; nv++) {
    ViewTypeX b_x("A", numCols, nv);
    ViewTypeY b_y("B", numRows, nv);
    ViewTypeY b_y_copy("B", numRows, nv);

    ViewTypeX b_xt("A", numRows, nv);
    ViewTypeY b_yt("B", numCols, nv);
    ViewTypeY b_yt_copy("B", numCols, nv);

    Kokkos::fill_random(b_x, rand_pool, scalar_t(10));
    Kokkos::fill_random(b_y, rand_pool, scalar_t(10));
    Kokkos::fill_random(b_xt, rand_pool, scalar_t(10));
    Kokkos::fill_random(b_yt, rand_pool, scalar_t(10));
    Kokkos::fill_random(input_mat.values, rand_pool, scalar_t(10));

    Kokkos::deep_copy(b_y_copy, b_y);
    Kokkos::deep_copy(b_yt_copy, b_yt);

    handle_t handle;

    Test::check_spmv_mv(&handle, input_mat, b_x, b_y, b_y_copy, 1.0, 0.0, nv, "N", max_nnz_per_row * max_val * max_x);
    Test::check_spmv_mv(&handle, input_mat, b_x, b_y, b_y_copy, 0.0, 1.0, nv, "N", max_y);
    Test::check_spmv_mv(&handle, input_mat, b_x, b_y, b_y_copy, 1.0, 1.0, nv, "N",
                        max_y + max_nnz_per_row * max_val * max_x);
    Test::check_spmv_mv(&handle, input_mat, b_xt, b_yt, b_yt_copy, 1.0, 0.0, nv, "T",
                        max_nnz_per_row * max_val * max_x);
    Test::check_spmv_mv(&handle, input_mat, b_xt, b_yt, b_yt_copy, 0.0, 1.0, nv, "T", max_y);
    // Testing all modes together, since matrix is square
    std::vector<const char *> modes   = {"N", "C", "T", "H"};
    std::vector<double> testAlphaBeta = {0.0, 1.0, -1.0, 2.5};
    for (auto mode : modes) {
      for (double alpha : testAlphaBeta) {
        for (double beta : testAlphaBeta) {
          mag_t max_error = beta * max_y + alpha * max_nnz_per_row * max_val * max_x;
          if (*mode == 'N' || *mode == 'C') {
            Test::check_spmv_mv(&handle, input_mat, b_x, b_y, b_y_copy, alpha, beta, nv, mode, max_error);
          } else {
            Test::check_spmv_mv(&handle, input_mat, b_xt, b_yt, b_yt_copy, alpha, beta, nv, mode, max_error);
          }
        }
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_1D(lno_t nx, lno_t leftBC, lno_t rightBC) {
  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using mag_t         = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  constexpr mag_t max_x   = static_cast<mag_t>(1);
  constexpr mag_t max_y   = static_cast<mag_t>(1);
  constexpr mag_t max_val = static_cast<mag_t>(2);

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t *[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 1);
  mat_structure(0, 0) = nx;
  if (leftBC == 1) {
    mat_structure(0, 1) = 1;
  }
  if (rightBC == 1) {
    mat_structure(0, 2) = 1;
  }

  crsMat_t input_mat = Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(input_x, rand_pool, max_x);
  Kokkos::fill_random(output_y, rand_pool, max_y);

  const mag_t max_error = max_y + 3 * max_val * max_x;

  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 0.0, max_error);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 0.0, 1.0, max_error);
  Test::check_spmv_struct(input_mat, 1, structure, input_x, output_y, 1.0, 1.0, max_error);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_2D(lno_t nx, lno_t ny, lno_t horizontalBC, lno_t verticalBC) {
  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using mag_t         = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  constexpr mag_t max_x = static_cast<mag_t>(1);
  constexpr mag_t max_y = static_cast<mag_t>(1);

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 2);
  structure(0) = nx;
  structure(1) = ny;
  Kokkos::View<lno_t *[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 2);
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

  crsMat_t input_mat_FD = Test::generate_structured_matrix2D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE = Test::generate_structured_matrix2D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(input_x, rand_pool, max_x);
  Kokkos::fill_random(output_y, rand_pool, max_y);

  {
    constexpr mag_t max_val   = static_cast<mag_t>(4);
    constexpr mag_t max_error = max_y + 5 * max_val * max_x;
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 0.0, max_error);
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0, 1.0, max_error);
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 1.0, max_error);
  }

  {
    constexpr mag_t max_val   = static_cast<mag_t>(8);
    constexpr mag_t max_error = max_y + 9 * max_val * max_x;
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 0.0, max_error);
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0, 1.0, max_error);
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 1.0, max_error);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spmv_struct_3D(lno_t nx, lno_t ny, lno_t nz, lno_t horizontal1BC, lno_t horizontal2BC, lno_t verticalBC) {
  using crsMat_t      = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using scalar_view_t = typename crsMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;
  using mag_t         = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  constexpr mag_t max_x = static_cast<mag_t>(1);
  constexpr mag_t max_y = static_cast<mag_t>(1);

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 3);
  structure(0) = nx;
  structure(1) = ny;
  structure(2) = nz;
  Kokkos::View<lno_t *[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 3);
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

  crsMat_t input_mat_FD = Test::generate_structured_matrix3D<crsMat_t>("FD", mat_structure);
  crsMat_t input_mat_FE = Test::generate_structured_matrix3D<crsMat_t>("FE", mat_structure);

  lno_t nr = input_mat_FD.numRows();
  lno_t nc = input_mat_FD.numCols();

  x_vector_type input_x("x", nc);
  y_vector_type output_y("y", nr);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(input_x, rand_pool, max_x);
  Kokkos::fill_random(output_y, rand_pool, max_y);

  {
    constexpr mag_t max_val   = static_cast<mag_t>(6);
    constexpr mag_t max_error = max_y + 7 * max_val * max_x;
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 0.0, max_error);
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 0.0, 1.0, max_error);
    Test::check_spmv_struct(input_mat_FD, 1, structure, input_x, output_y, 1.0, 1.0, max_error);
  }

  {
    constexpr mag_t max_val   = static_cast<mag_t>(26);
    constexpr mag_t max_error = max_y + 27 * max_val * max_x;
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 0.0, max_error);
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 0.0, 1.0, max_error);
    Test::check_spmv_struct(input_mat_FE, 2, structure, input_x, output_y, 1.0, 1.0, max_error);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename layout, class Device>
void test_spmv_mv_struct_1D(lno_t nx, int numMV) {
  using crsMat_t           = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type>;
  using x_multivector_type = Kokkos::View<scalar_t **, layout, Device>;
  using y_multivector_type = Kokkos::View<scalar_t **, layout, Device>;
  using mag_t              = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  constexpr mag_t max_x = static_cast<mag_t>(1);
  constexpr mag_t max_y = static_cast<mag_t>(1);

  Kokkos::View<lno_t *, Kokkos::HostSpace> structure("Spmv Structure", 1);
  structure(0) = nx;
  Kokkos::View<lno_t *[3], Kokkos::HostSpace> mat_structure("Matrix Structure", 1);
  mat_structure(0, 0) = nx;
  mat_structure(0, 1) = 1;
  mat_structure(0, 2) = 1;

  crsMat_t input_mat = Test::generate_structured_matrix1D<crsMat_t>(mat_structure);

  lno_t nr = input_mat.numRows();
  lno_t nc = input_mat.numCols();

  x_multivector_type input_x("x", nc, numMV);
  y_multivector_type output_y("y", nr, numMV);
  y_multivector_type output_y_copy("y_copy", nr, numMV);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(input_x, rand_pool, max_x);
  Kokkos::fill_random(output_y, rand_pool, max_y);

  constexpr mag_t max_error = 5;

  Kokkos::deep_copy(output_y_copy, output_y);

  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 1.0, 0.0, numMV, max_error);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 0.0, 1.0, numMV, max_error);
  Test::check_spmv_mv_struct(input_mat, 1, structure, input_x, output_y, output_y_copy, 1.0, 1.0, numMV, max_error);
}

// call it if ordinal int and, scalar float and double are instantiated.
template <class DeviceType>
void test_github_issue_101() {
  typedef KokkosSparse::CrsMatrix<float, int, DeviceType> float_matrix_type;
  typedef KokkosSparse::CrsMatrix<double, int, DeviceType> double_matrix_type;
  static_assert(std::is_same<typename float_matrix_type::StaticCrsGraphType,
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

    double_matrix_type A_d("A_d", G, numCols);
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
        static_cast<double>(1.0) + static_cast<double>(EPS_f) / static_cast<double>(2.0);
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

    float_matrix_type A_f("A_f", G, numCols);
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
        static_cast<double>(1.0) + static_cast<double>(EPS_f) / static_cast<double>(2.0);
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

template <class scalar_t, class lno_t, class size_type, class layout_t, class DeviceType>
void test_spmv_all_interfaces_light() {
  // Using a small matrix, run through the various SpMV interfaces and
  // make sure they produce the correct results.
  using execution_space = typename DeviceType::execution_space;
  using mag_t           = typename Kokkos::ArithTraits<scalar_t>::mag_type;
  using crsMat_t        = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, DeviceType, void, size_type>;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  const lno_t m      = 111;
  const lno_t n      = 99;
  const mag_t maxVal = 10.0;
  const mag_t eps    = 10.0 * Kokkos::ArithTraits<mag_t>::eps();
  size_type nnz      = 600;
  crsMat_t A         = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, n, nnz, 2, lno_t(n * 0.7));
  // note: A's values are in range [0, 50)
  const mag_t maxError = (nnz / m) * 50.0 * maxVal;
  using multivector_t  = Kokkos::View<scalar_t **, layout_t, DeviceType>;
  using vector_t       = Kokkos::View<scalar_t *, layout_t, DeviceType>;
  using range1D_t      = Kokkos::RangePolicy<execution_space>;
  using range2D_t      = Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>;
  using v_handle_t     = KokkosSparse::SPMVHandle<DeviceType, crsMat_t, vector_t, vector_t>;
  using mv_handle_t    = KokkosSparse::SPMVHandle<DeviceType, crsMat_t, multivector_t, multivector_t>;
  multivector_t x_mv("x_mv", n, 3);
  vector_t x("x", n);
  // Randomize x (it won't be modified after that)
  Kokkos::fill_random(x_mv, rand_pool, randomUpperBound<scalar_t>(maxVal));
  Kokkos::fill_random(x, rand_pool, randomUpperBound<scalar_t>(maxVal));
  multivector_t y_mv("y_mv", m, 3);
  vector_t y("y", m);
  // Compute the correct y = Ax once
  multivector_t ygold_mv("ygold_mv", m, 3);
  vector_t ygold("ygold", m);
  for (lno_t i = 0; i < 3; i++)
    Test::sequential_spmv(A, Kokkos::subview(x_mv, Kokkos::ALL(), i), Kokkos::subview(ygold_mv, Kokkos::ALL(), i), 1.0,
                          0.0);
  Test::sequential_spmv(A, x, ygold, 1.0, 0.0);
  auto clear_y = [&]() { Kokkos::deep_copy(y_mv, scalar_t(0)); };
  auto verify  = [&]() {
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosSparse::Test::spmv", range1D_t(0, m),
                             Test::fSPMV<vector_t, vector_t>(ygold, y, eps, maxError), num_errors);
    EXPECT_EQ(num_errors, 0);
  };
  auto verify_mv = [&]() {
    int num_errors = 0;
    Kokkos::parallel_reduce("KokkosSparse::Test::spmv", range2D_t({0, 0}, {m, 3}),
                            Test::fSPMV<multivector_t, multivector_t>(ygold_mv, y_mv, eps, maxError), num_errors);
    EXPECT_EQ(num_errors, 0);
  };
  // Now run through the interfaces and check results each time.
  execution_space space;
  std::vector<execution_space> space_partitions;
  if (space.concurrency() > 1) {
    space_partitions = Kokkos::Experimental::partition_space(space, 1, 1);
    space            = space_partitions[1];
  }

  v_handle_t v_handle;
  mv_handle_t mv_handle;

  // space and handle
  spmv(space, &v_handle, "N", 1.0, A, x, 0.0, y);
  space.fence();
  verify();
  clear_y();
  spmv(space, &mv_handle, "N", 1.0, A, x_mv, 0.0, y_mv);
  space.fence();
  verify_mv();
  clear_y();
  // handle
  spmv(&v_handle, "N", 1.0, A, x, 0.0, y);
  verify();
  clear_y();
  spmv(&mv_handle, "N", 1.0, A, x_mv, 0.0, y_mv);
  verify_mv();
  clear_y();
  // space
  spmv(space, "N", 1.0, A, x, 0.0, y);
  space.fence();
  verify();
  clear_y();
  spmv(space, "N", 1.0, A, x_mv, 0.0, y_mv);
  space.fence();
  verify_mv();
  clear_y();
  // neither
  spmv("N", 1.0, A, x, 0.0, y);
  verify();
  clear_y();
  spmv("N", 1.0, A, x_mv, 0.0, y_mv);
  verify_mv();
  clear_y();
}

#define EXECUTE_TEST_ISSUE_101(DEVICE) \
  TEST_F(TestCategory, sparse##_##spmv_issue_101##_##OFFSET##_##DEVICE) { test_github_issue_101<DEVICE>(); }

#define EXECUTE_TEST_FN(SCALAR, ORDINAL, OFFSET, DEVICE)                                     \
  TEST_F(TestCategory, sparse##_##spmv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {       \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 3, 200, 10, true);    \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 3, 100, 10, true);    \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 20, 100, 5, true);    \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 3, 20, 10, false);  \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(50000, 50000 * 3, 100, 10, false); \
    test_spmv_algorithms<SCALAR, ORDINAL, OFFSET, DEVICE>(10000, 10000 * 2, 100, 5, false);  \
  }

#define EXECUTE_TEST_INTERFACES(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                               \
  TEST_F(TestCategory, sparse_spmv_interfaces_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_spmv_all_interfaces_light<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>();                 \
  }

#define EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                                                   \
  TEST_F(TestCategory, sparse##_##spmv_mv##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) {               \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(1001, 1001 * 3, 200, 10, true, 1);               \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(999, 999 * 3, 100, 10, true, 5);                 \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(1003, 1003 * 2, 100, 5, true, 10);               \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(50007, 50007 * 3, 20, 10, false, 1);             \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(50002, 50002 * 3, 100, 10, false, 1);            \
    test_spmv_mv<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(10000, 10000 * 2, 100, 5, false, 5);             \
    test_spmv_mv_heavy<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, Kokkos::LAYOUT, DEVICE>(204, 201, 204 * 10, 60, 4, \
                                                                                        30);                       \
    test_spmv_mv_heavy<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, Kokkos::LAYOUT, DEVICE>(2, 3, 5, 3, 1, 10);        \
  }

#define EXECUTE_TEST_MV_MIXED_LAYOUT(SCALAR, ORDINAL, OFFSET, DEVICE)                                               \
  TEST_F(TestCategory, sparse##_##spmv_mv_mixed_layout##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) {   \
    test_spmv_mv_heavy<SCALAR, ORDINAL, OFFSET, Kokkos::LayoutRight, Kokkos::LayoutLeft, DEVICE>(99, 101, 100 * 15, \
                                                                                                 40, 4, 20);        \
  }

#define EXECUTE_TEST_STRUCT(SCALAR, ORDINAL, OFFSET, DEVICE)                                  \
  TEST_F(TestCategory, sparse##_##spmv_struct##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_spmv_struct_1D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 1, 1);                           \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(25, 21, 3, 3);                       \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(20, 25, 3, 3);                       \
    test_spmv_struct_2D<SCALAR, ORDINAL, OFFSET, DEVICE>(22, 22, 3, 3);                       \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(20, 20, 20, 3, 3, 3);                \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(22, 22, 22, 3, 3, 3);                \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(25, 10, 20, 3, 3, 3);                \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 20, 25, 3, 3, 3);                \
    test_spmv_struct_3D<SCALAR, ORDINAL, OFFSET, DEVICE>(10, 24, 20, 3, 3, 3);                \
  }

#define EXECUTE_TEST_MV_STRUCT(SCALAR, ORDINAL, OFFSET, LAYOUT, DEVICE)                                     \
  TEST_F(TestCategory, sparse##_##spmv_mv_struct##_##SCALAR##_##ORDINAL##_##OFFSET##_##LAYOUT##_##DEVICE) { \
    test_spmv_mv_struct_1D<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(10, 1);                         \
    test_spmv_mv_struct_1D<SCALAR, ORDINAL, OFFSET, Kokkos::LAYOUT, DEVICE>(10, 2);                         \
  }

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST_ISSUE_101(TestDevice)
#endif

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_FN(SCALAR, ORDINAL, OFFSET, TestDevice)              \
  EXECUTE_TEST_STRUCT(SCALAR, ORDINAL, OFFSET, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)       \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestDevice)        \
  EXECUTE_TEST_MV_STRUCT(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestDevice) \
  EXECUTE_TEST_INTERFACES(SCALAR, ORDINAL, OFFSET, LayoutLeft, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST

#endif  // defined(KOKKOSKERNELS_INST_LAYOUTLEFT)

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV(SCALAR, ORDINAL, OFFSET, LayoutRight, TestDevice) \
  EXECUTE_TEST_INTERFACES(SCALAR, ORDINAL, OFFSET, LayoutRight, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#endif

// Test that requires mixing LayoutLeft and LayoutRight (never an ETI'd
// combination)
#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
  EXECUTE_TEST_MV_MIXED_LAYOUT(SCALAR, ORDINAL, OFFSET, TestDevice)

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#endif

#undef EXECUTE_TEST_FN
#undef EXECUTE_TEST_STRUCT
#undef EXECUTE_TEST_MV
#undef EXECUTE_TEST_MV_STRUCT
