// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosSparse_spmv.hpp>

namespace Test {

template <class AMatrix, class XVector, class YVector>
void serial_sell_spmv(const typename YVector::non_const_value_type alpha, const AMatrix& A, const XVector& x,
                      const typename YVector::non_const_value_type beta, const YVector& y) {
  // Copy vectors to host
  typename XVector::host_mirror_type x_h = Kokkos::create_mirror_view(x);
  typename YVector::host_mirror_type y_h = Kokkos::create_mirror_view(y);

  Kokkos::deep_copy(x_h, x);
  Kokkos::deep_copy(y_h, y);

  // Copy matrix data to host
  typename AMatrix::entries_type::host_mirror_type entries_h = Kokkos::create_mirror_view(A.entries);
  typename AMatrix::values_type::host_mirror_type values_h   = Kokkos::create_mirror_view(A.values);

  Kokkos::deep_copy(entries_h, A.entries);
  Kokkos::deep_copy(values_h, A.values);

  typename YVector::non_const_value_type sum = 0;
  const int row_length                       = A.sell_nnz / A.num_rows_per_slice;
  for (int rowIdx = 0; rowIdx < A.num_rows; ++rowIdx) {
    y_h(rowIdx) = beta * y_h(rowIdx);

    sum = 0;
    for (int colIdx = 0; colIdx < row_length; ++colIdx) {
      if (entries_h(colIdx * A.num_rows_per_slice + rowIdx) > -1) {
        sum +=
            values_h(colIdx * A.num_rows_per_slice + rowIdx) * x_h(entries_h(colIdx * A.num_rows_per_slice + rowIdx));
      }
    }
    y_h(rowIdx) += alpha * sum;
  }

  Kokkos::deep_copy(y, y_h);
}

template <typename scalar_t, typename ordinal_t, typename size_type, typename Device>
void test_spmv_sell_analytic() {
  using sellMat_t     = KokkosSparse::Experimental::SellMatrix<scalar_t, ordinal_t, Device, void, size_type>;
  using scalar_view_t = typename sellMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;

  using ordinal_type = typename sellMat_t::ordinal_type;
  using value_type   = typename sellMat_t::value_type;

  using KAT = KokkosKernels::ArithTraits<value_type>;

  using offsets_type = typename sellMat_t::offsets_type;
  using entries_type = typename sellMat_t::entries_type;
  using values_type  = typename sellMat_t::values_type;

  //       0  1  2  3  4  5  6
  //     [ 2 -1               ]             [ 0  1 -1]            [ 2.0 -1.0  0.0]
  //     [-1  2 -1            ]             [ 0  1  2]            [-1.0  2.0 -1.0]
  //     [   -1  2 -1         ]             [ 1  2  3]            [-1.0  2.0 -1.0]
  // A = [      -1  2 -1      ]   colinds = [ 2  3  4]   values = [-1.0  2.0 -1.0]
  //     [         -1  2 -1   ]             [ 3  4  5]            [-1.0  2.0 -1.0]
  //     [            -1  2 -1]             [ 4  5  6]            [-1.0  2.0 -1.0]
  //     [               -1  2]             [ 5  6 -1]            [-1.0  2.0  0.0]

  constexpr ordinal_type nrows = 7, ncols = 7, rows_per_slice = 7, nslices = 1;
  constexpr size_type nnz = 19, sell_nnz = 21;

  offsets_type slice_offsets("slice offsets", nslices + 1);
  typename offsets_type::host_mirror_type slice_offsets_h = Kokkos::create_mirror_view(slice_offsets);
  slice_offsets_h(1)                                      = sell_nnz;
  Kokkos::deep_copy(slice_offsets, slice_offsets_h);

  entries_type colinds("column indices", sell_nnz);
  typename entries_type::host_mirror_type colinds_h = Kokkos::create_mirror_view(colinds);
  colinds_h(0)                                      = 0;
  colinds_h(1)                                      = 0;
  colinds_h(2)                                      = 1;
  colinds_h(3)                                      = 2;
  colinds_h(4)                                      = 3;
  colinds_h(5)                                      = 4;
  colinds_h(6)                                      = 5;
  colinds_h(7)                                      = 1;
  colinds_h(8)                                      = 1;
  colinds_h(9)                                      = 2;
  colinds_h(10)                                     = 3;
  colinds_h(11)                                     = 4;
  colinds_h(12)                                     = 5;
  colinds_h(13)                                     = 6;
  colinds_h(14)                                     = -1;
  colinds_h(15)                                     = 2;
  colinds_h(16)                                     = 3;
  colinds_h(17)                                     = 4;
  colinds_h(18)                                     = 5;
  colinds_h(19)                                     = 6;
  colinds_h(20)                                     = -1;
  Kokkos::deep_copy(colinds, colinds_h);

  values_type values("values", sell_nnz);
  typename values_type::host_mirror_type values_h = Kokkos::create_mirror_view(values);
  const value_type zero                           = KAT::zero();
  const value_type one                            = KAT::one();
  const value_type two                            = KAT::one() + KAT::one();
  values_h(0)                                     = two;
  values_h(1)                                     = -one;
  values_h(2)                                     = -one;
  values_h(3)                                     = -one;
  values_h(4)                                     = -one;
  values_h(5)                                     = -one;
  values_h(6)                                     = -one;
  values_h(7)                                     = -one;
  values_h(8)                                     = two;
  values_h(9)                                     = two;
  values_h(10)                                    = two;
  values_h(11)                                    = two;
  values_h(12)                                    = two;
  values_h(13)                                    = two;
  values_h(14)                                    = zero;
  values_h(15)                                    = -one;
  values_h(16)                                    = -one;
  values_h(17)                                    = -one;
  values_h(18)                                    = -one;
  values_h(19)                                    = -one;
  values_h(20)                                    = zero;
  Kokkos::deep_copy(values, values_h);

  sellMat_t A(nrows, ncols, nnz, sell_nnz, rows_per_slice, slice_offsets, colinds, values);

  // Create and initialize x vector
  x_vector_type x("x vector", ncols);
  typename x_vector_type::host_mirror_type x_h = Kokkos::create_mirror_view(x);
  for (ordinal_type idx = 0; idx < ncols; ++idx) {
    x_h(idx) = idx * idx;
  }
  Kokkos::deep_copy(x, x_h);

  // Create and initialize y vector
  y_vector_type y("y vector", nrows);
  typename y_vector_type::host_mirror_type y_h = Kokkos::create_mirror_view(y);

  // Compute matrix vector product
  KokkosSparse::Experimental::spmv("N", 1.0, A, x, 1.0, y);

  // Copy results back to the host
  Kokkos::deep_copy(y_h, y);

  constexpr int row_length = 3;
  Test::EXPECT_NEAR_KK_REL(y_h(0), scalar_t(-1), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(1), scalar_t(-2), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(2), scalar_t(-2), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(3), scalar_t(-2), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(4), scalar_t(-2), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(5), scalar_t(-2), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  Test::EXPECT_NEAR_KK_REL(y_h(6), scalar_t(47), 10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
}

template <typename scalar_t, typename ordinal_t, typename size_type, typename Device>
void test_spmv_sell(ordinal_t num_rows, ordinal_t num_cols, ordinal_t row_length, ordinal_t variance,
                    ordinal_t bandwidth, scalar_t alpha, scalar_t beta) {
  using sellMat_t     = KokkosSparse::Experimental::SellMatrix<scalar_t, ordinal_t, Device, void, size_type>;
  using scalar_view_t = typename sellMat_t::values_type::non_const_type;
  using x_vector_type = scalar_view_t;
  using y_vector_type = scalar_view_t;

  sellMat_t matA = KokkosSparse::Impl::kk_generate_sell_sparse_matrix<sellMat_t>(num_rows, num_cols, row_length,
                                                                                 variance, bandwidth);

  x_vector_type x("x vector", num_cols);
  y_vector_type y("y vector", num_rows);
  y_vector_type y_ref("y vector ref", num_rows);

  {
    typename x_vector_type::host_mirror_type x_h = Kokkos::create_mirror_view(x);
    typename y_vector_type::host_mirror_type y_h = Kokkos::create_mirror_view(y);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> values_distribution(-10, 10);

    for (int colIdx = 0; colIdx < matA.num_cols; ++colIdx) {
      x_h(colIdx) = values_distribution(gen);
    }

    for (int rowIdx = 0; rowIdx < matA.num_rows; ++rowIdx) {
      y_h(rowIdx) = values_distribution(gen);
    }

    Kokkos::deep_copy(x, x_h);
    Kokkos::deep_copy(y, y_h);
  }

  Kokkos::deep_copy(y_ref, y);

  // Compute matrix vector product
  KokkosSparse::Experimental::spmv("N", alpha, matA, x, beta, y);

  serial_sell_spmv(alpha, matA, x, beta, y_ref);

  typename y_vector_type::host_mirror_type y_h = Kokkos::create_mirror_view(y);
  Kokkos::deep_copy(y_h, y);

  typename y_vector_type::host_mirror_type y_ref_h = Kokkos::create_mirror_view(y_ref);
  Kokkos::deep_copy(y_ref_h, y_ref);

  for (int rowIdx = 0; rowIdx < num_rows; ++rowIdx) {
    Test::EXPECT_NEAR_KK_REL(y_h(rowIdx), y_ref_h(rowIdx),
                             10 * row_length * KokkosKernels::ArithTraits<scalar_t>::eps());
  }
}

}  // namespace Test

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                         \
  TEST_F(TestCategory, sparse##_##spmv_sell##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    ::Test::test_spmv_sell_analytic<SCALAR, ORDINAL, OFFSET, DEVICE>();                     \
    ::Test::test_spmv_sell<SCALAR, ORDINAL, OFFSET, DEVICE>(5, 5, 3, 0, 2, 2.0, 1.0);       \
    ::Test::test_spmv_sell<SCALAR, ORDINAL, OFFSET, DEVICE>(20, 20, 5, 2, 2, 2.0, 1.0);     \
    ::Test::test_spmv_sell<SCALAR, ORDINAL, OFFSET, DEVICE>(100, 100, 7, 3, 10, 2.0, 1.0);  \
    ::Test::test_spmv_sell<SCALAR, ORDINAL, OFFSET, DEVICE>(100, 50, 7, 2, 3, 2.0, 1.0);    \
    ::Test::test_spmv_sell<SCALAR, ORDINAL, OFFSET, DEVICE>(80, 100, 7, 2, 3, 2.0, 1.0);    \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
