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
#include "KokkosSparse_mdf.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

namespace Test {

template <typename scalar_type, typename ordinal_type, typename size_type, typename device>
void run_test_mdf() {
  using crs_matrix_type = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, device, void, size_type>;
  using crs_graph_type  = typename crs_matrix_type::StaticCrsGraphType;
  using row_map_type    = typename crs_graph_type::row_map_type::non_const_type;
  using col_ind_type    = typename crs_graph_type::entries_type::non_const_type;
  using values_type     = typename crs_matrix_type::values_type::non_const_type;
  using value_type      = typename crs_matrix_type::value_type;

  const value_type four = static_cast<value_type>(4.0);

  constexpr ordinal_type numRows  = 16;
  constexpr ordinal_type numCols  = 16;
  constexpr size_type numNonZeros = 64;
  row_map_type row_map("row map", numRows + 1);
  col_ind_type col_ind("column indices", numNonZeros);
  values_type values("values", numNonZeros);

  {  // create matrix
    const size_type row_mapRaw[]    = {0, 3, 7, 11, 14, 18, 23, 28, 32, 36, 41, 46, 50, 53, 57, 61, 64};
    const ordinal_type col_indRaw[] = {0,  1,  4, 0,  1,  2, 5,  1,  2,  3,  6,  2,  3,  7,  0,  4,
                                       5,  8,  1, 4,  5,  6, 9,  2,  5,  6,  7,  10, 3,  6,  7,  11,
                                       4,  8,  9, 12, 5,  8, 9,  10, 13, 6,  9,  10, 11, 14, 7,  10,
                                       11, 15, 8, 12, 13, 9, 12, 13, 14, 10, 13, 14, 15, 11, 14, 15};
    const value_type values_Raw[]   = {4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, 4,
                                       -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1,
                                       -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1,
                                       4,  -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4};

    typename row_map_type::HostMirror::const_type row_map_host(row_mapRaw, numRows + 1);
    typename col_ind_type::HostMirror::const_type col_ind_host(col_indRaw, numNonZeros);
    typename values_type::HostMirror::const_type values_host(values_Raw, numNonZeros);

    Kokkos::deep_copy(row_map, row_map_host);
    Kokkos::deep_copy(col_ind, col_ind_host);
    Kokkos::deep_copy(values, values_host);
  }

  crs_matrix_type A = crs_matrix_type("A", numRows, numCols, numNonZeros, values, row_map, col_ind);

  KokkosSparse::Experimental::MDF_handle<crs_matrix_type> handle(A);
  handle.set_verbosity(0);
  KokkosSparse::Experimental::mdf_symbolic(A, handle);
  KokkosSparse::Experimental::mdf_numeric(A, handle);

  col_ind_type permutation = handle.get_permutation();

  bool success                                    = true;
  typename col_ind_type::HostMirror permutation_h = Kokkos::create_mirror(permutation);
  Kokkos::deep_copy(permutation_h, permutation);
  const ordinal_type permutation_ref[] = {0, 3, 12, 15, 1, 2, 4, 8, 7, 11, 13, 14, 5, 6, 9, 10};
  printf("MDF ordering: { ");
  for (ordinal_type idx = 0; idx < A.numRows(); ++idx) {
    printf("%d ", static_cast<int>(permutation_h(idx)));
    if (permutation_h(idx) != permutation_ref[idx]) {
      success = false;
    }
  }
  printf("}\n");
  EXPECT_TRUE(success) << "The permutation computed is different from the reference solution!";

  // Check the factors L and U
  handle.sort_factors();
  crs_matrix_type U = handle.getU();
  crs_matrix_type L = handle.getL();

  EXPECT_TRUE(U.numRows() == 16);
  EXPECT_TRUE(U.nnz() == 40);

  {
    auto row_map_U = Kokkos::create_mirror(U.graph.row_map);
    Kokkos::deep_copy(row_map_U, U.graph.row_map);
    auto entries_U = Kokkos::create_mirror(U.graph.entries);
    Kokkos::deep_copy(entries_U, U.graph.entries);
    auto values_U = Kokkos::create_mirror(U.values);
    Kokkos::deep_copy(values_U, U.values);

    const size_type row_map_U_ref[17]    = {0, 3, 6, 9, 12, 15, 17, 20, 22, 25, 27, 30, 32, 35, 37, 39, 40};
    const ordinal_type entries_U_ref[40] = {0, 4,  6, 1, 5,  8, 2,  7,  10, 3,  9,  11, 4,  5,  12, 5,  13, 6,  7,  12,
                                            7, 14, 8, 9, 13, 9, 15, 10, 11, 14, 11, 15, 12, 13, 14, 13, 15, 14, 15, 15};

    const scalar_type val0             = static_cast<scalar_type>(15. / 4.);
    const scalar_type val1             = static_cast<scalar_type>(val0 - 1 / val0);
    const scalar_type val2             = static_cast<scalar_type>(4 - 2 / val0);
    const scalar_type val3             = static_cast<scalar_type>(4 - 1 / val0 - 1 / val1 - 1 / val2);
    const scalar_type val4             = static_cast<scalar_type>(4 - 2 / val1 - 2 / val3);
    const scalar_type values_U_ref[40] = {4,  -1,   -1,   4,    -1,   -1, 4,    -1,   -1,   4,    -1, -1,   val0, -1,
                                          -1, val1, -1,   val0, -1,   -1, val1, -1,   val0, -1,   -1, val1, -1,   val0,
                                          -1, -1,   val1, -1,   val2, -1, -1,   val3, -1,   val3, -1, val4};

    for (int idx = 0; idx < 17; ++idx) {
      EXPECT_TRUE(row_map_U_ref[idx] == row_map_U(idx)) << "rowmap_U(" << idx << ") is wrong!";
    }
    for (int idx = 0; idx < 40; ++idx) {
      EXPECT_TRUE(entries_U_ref[idx] == entries_U(idx)) << "entries_U(" << idx << ") is wrong!";
      EXPECT_NEAR_KK(values_U_ref[idx], values_U(idx), 10 * Kokkos::ArithTraits<scalar_type>::eps(),
                     "An entry in U.values is wrong!");
    }

    auto row_map_L = Kokkos::create_mirror(L.graph.row_map);
    Kokkos::deep_copy(row_map_L, L.graph.row_map);
    auto entries_L = Kokkos::create_mirror(L.graph.entries);
    Kokkos::deep_copy(entries_L, L.graph.entries);
    auto values_L = Kokkos::create_mirror(L.values);
    Kokkos::deep_copy(values_L, L.values);

    const size_type row_map_L_ref[17]    = {0, 1, 2, 3, 4, 6, 9, 11, 14, 16, 19, 21, 24, 27, 31, 35, 40};
    const ordinal_type entries_L_ref[40] = {0,  1, 2,  3,  0, 4, 1,  4, 5, 0,  6,  2, 6,  7,  1,  8, 3,  8,  9,  2,
                                            10, 3, 10, 11, 4, 6, 12, 5, 8, 12, 13, 7, 10, 12, 14, 9, 11, 13, 14, 15};
    const scalar_type values_L_ref[40]   = {
        1, 1,         1,         1,         -1 / four, 1,         -1 / four, -1 / val0, 1,         -1 / four,
        1, -1 / four, -1 / val0, 1,         -1 / four, 1,         -1 / four, -1 / val0, 1,         -1 / four,
        1, -1 / four, -1 / val0, 1,         -1 / val0, -1 / val0, 1,         -1 / val1, -1 / val0, -1 / val2,
        1, -1 / val1, -1 / val0, -1 / val2, 1,         -1 / val1, -1 / val1, -1 / val3, -1 / val3, 1};

    for (int idx = 0; idx < 17; ++idx) {
      EXPECT_TRUE(row_map_L_ref[idx] == row_map_L(idx)) << "rowmap_L(" << idx << ")=" << row_map_L(idx) << " is wrong!";
    }
    for (int idx = 0; idx < 40; ++idx) {
      EXPECT_TRUE(entries_L_ref[idx] == entries_L(idx))
          << "entries_L(" << idx << ")=" << entries_L(idx) << " is wrong, entries_L_ref[" << idx
          << "]=" << entries_L_ref[idx] << "!";
      EXPECT_NEAR_KK(values_L_ref[idx], values_L(idx), 10 * Kokkos::ArithTraits<scalar_type>::eps(),
                     "An entry in L.values is wrong!");
    }
  }
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_mdf() {
  Test::run_test_mdf<scalar_t, lno_t, size_type, device>();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                   \
  TEST_F(TestCategory, sparse##_##mdf##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_mdf<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
