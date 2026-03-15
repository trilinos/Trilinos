// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "KokkosSparse_SellMatrix.hpp"
#include "KokkosKernels_ArithTraits.hpp"

namespace Test {

template <class SellMatrixType>
void SellConstructor() {
  using ordinal_type = typename SellMatrixType::ordinal_type;
  using size_type    = typename SellMatrixType::size_type;
  using value_type   = typename SellMatrixType::value_type;

  using KAT = KokkosKernels::ArithTraits<value_type>;

  using offsets_type = typename SellMatrixType::offsets_type;
  using entries_type = typename SellMatrixType::entries_type;
  using values_type  = typename SellMatrixType::values_type;

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
  colinds_h(7)                                      = 1;
  colinds_h(14)                                     = -1;
  colinds_h(1)                                      = 0;
  colinds_h(8)                                      = 1;
  colinds_h(15)                                     = 2;
  colinds_h(2)                                      = 1;
  colinds_h(9)                                      = 2;
  colinds_h(16)                                     = 3;
  colinds_h(3)                                      = 2;
  colinds_h(10)                                     = 3;
  colinds_h(17)                                     = 4;
  colinds_h(4)                                      = 3;
  colinds_h(11)                                     = 4;
  colinds_h(17)                                     = 5;
  colinds_h(5)                                      = 4;
  colinds_h(12)                                     = 5;
  colinds_h(19)                                     = 6;
  colinds_h(6)                                      = 5;
  colinds_h(13)                                     = 6;
  colinds_h(20)                                     = -1;
  Kokkos::deep_copy(colinds, colinds_h);

  values_type values("values", sell_nnz);
  typename values_type::host_mirror_type values_h = Kokkos::create_mirror_view(values);
  const value_type zero                           = KAT::zero();
  const value_type one                            = KAT::one();
  const value_type two                            = KAT::one() + KAT::one();
  values_h(0)                                     = two;
  values_h(7)                                     = -one;
  values_h(14)                                    = zero;
  values_h(1)                                     = -one;
  values_h(8)                                     = two;
  values_h(15)                                    = -one;
  values_h(2)                                     = -one;
  values_h(9)                                     = two;
  values_h(16)                                    = -one;
  values_h(3)                                     = -one;
  values_h(10)                                    = two;
  values_h(17)                                    = -one;
  values_h(4)                                     = -one;
  values_h(11)                                    = two;
  values_h(18)                                    = -one;
  values_h(5)                                     = -one;
  values_h(12)                                    = two;
  values_h(19)                                    = -one;
  values_h(6)                                     = -one;
  values_h(13)                                    = two;
  values_h(20)                                    = zero;
  Kokkos::deep_copy(values, values_h);

  SellMatrixType A(nrows, ncols, nnz, sell_nnz, rows_per_slice, slice_offsets, colinds, values);

  EXPECT_EQ(A.num_slices, 1);

  // Test for bad input detection

  // Number of rows larger than rows_per_slice
  // This one should be removed once we support multiple slices
  EXPECT_THROW(
      { SellMatrixType A_bad_nrows(nrows + 1, ncols, nnz, sell_nnz, rows_per_slice, slice_offsets, colinds, values); },
      std::invalid_argument);

  // Mismatch between slice_offsets.extent(0) and number of slices
  EXPECT_THROW(
      {
        offsets_type bad_slice_offsets("slice offsets", nslices + 2);
        SellMatrixType A_bad_nrows(nrows, ncols, nnz, sell_nnz, rows_per_slice, bad_slice_offsets, colinds, values);
      },
      std::invalid_argument);

  // nnz is larger than sell_nnz
  EXPECT_THROW(
      {
        SellMatrixType A_bad_nrows(nrows, ncols, nnz + sell_nnz + 1, sell_nnz, rows_per_slice, slice_offsets, colinds,
                                   values);
      },
      std::invalid_argument);

  // Mismatch between sell_nnz and size of colinds and values
  EXPECT_THROW(
      { SellMatrixType A_bad_nrows(nrows, ncols, nnz, sell_nnz - 1, rows_per_slice, slice_offsets, colinds, values); },
      std::invalid_argument);
}
}  // namespace Test

template <typename Scalar, typename Ordinal, typename Offset, class Device>
void testSellMatrix() {
  using sell_matrix_type = KokkosSparse::Experimental::SellMatrix<Scalar, Ordinal, Device, void, Offset>;
  Test::SellConstructor<sell_matrix_type>();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                          \
  TEST_F(TestCategory, sparse##_##sellmatrix##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    testSellMatrix<SCALAR, ORDINAL, OFFSET, DEVICE>();                                       \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
