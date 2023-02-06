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

namespace Test {

template <typename scalar_type, typename ordinal_type, typename size_type,
          typename device>
void run_test_mdf() {
  using crs_matrix_type = KokkosSparse::CrsMatrix<scalar_type, ordinal_type,
                                                  device, void, size_type>;
  using crs_graph_type  = typename crs_matrix_type::StaticCrsGraphType;
  using row_map_type    = typename crs_graph_type::row_map_type::non_const_type;
  using col_ind_type    = typename crs_graph_type::entries_type::non_const_type;
  using values_type     = typename crs_matrix_type::values_type::non_const_type;
  using value_type      = typename crs_matrix_type::value_type;

  constexpr ordinal_type numRows  = 16;
  constexpr ordinal_type numCols  = 16;
  constexpr size_type numNonZeros = 64;
  row_map_type row_map("row map", numRows + 1);
  col_ind_type col_ind("column indices", numNonZeros);
  values_type values("values", numNonZeros);

  {  // create matrix
    const size_type row_mapRaw[]    = {0,  3,  7,  11, 14, 18, 23, 28, 32,
                                    36, 41, 46, 50, 53, 57, 61, 64};
    const ordinal_type col_indRaw[] = {
        0,  1,  4, 0,  1,  2, 5,  1,  2,  3,  6,  2,  3,  7,  0,  4,
        5,  8,  1, 4,  5,  6, 9,  2,  5,  6,  7,  10, 3,  6,  7,  11,
        4,  8,  9, 12, 5,  8, 9,  10, 13, 6,  9,  10, 11, 14, 7,  10,
        11, 15, 8, 12, 13, 9, 12, 13, 14, 10, 13, 14, 15, 11, 14, 15};
    const value_type values_Raw[] = {
        4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, 4,
        -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1,
        -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1, 4,  -1, -1, -1, -1,
        4,  -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4,  -1, -1, -1, 4};

    typename row_map_type::HostMirror::const_type row_map_host(row_mapRaw,
                                                               numRows + 1);
    typename col_ind_type::HostMirror::const_type col_ind_host(col_indRaw,
                                                               numNonZeros);
    typename values_type::HostMirror::const_type values_host(values_Raw,
                                                             numNonZeros);

    Kokkos::deep_copy(row_map, row_map_host);
    Kokkos::deep_copy(col_ind, col_ind_host);
    Kokkos::deep_copy(values, values_host);
  }

  crs_matrix_type A = crs_matrix_type("A", numRows, numCols, numNonZeros,
                                      values, row_map, col_ind);

  KokkosSparse::Experimental::MDF_handle<crs_matrix_type> handle(A);
  handle.set_verbosity(0);
  mdf_symbolic_phase(A, handle);
  mdf_numeric_phase(A, handle);

  col_ind_type permutation = handle.get_permutation();

  bool success = true;
  typename col_ind_type::HostMirror permutation_h =
      Kokkos::create_mirror(permutation);
  Kokkos::deep_copy(permutation_h, permutation);
  const ordinal_type permutation_ref[] = {0, 3,  12, 15, 1, 2, 4, 8,
                                          7, 11, 13, 14, 5, 6, 9, 10};
  printf("MDF ordering: { ");
  for (ordinal_type idx = 0; idx < A.numRows(); ++idx) {
    ;
    printf("%d ", static_cast<int>(permutation_h(idx)));
    if (permutation_h(idx) != permutation_ref[idx]) {
      success = false;
    }
  }
  printf("}\n");

  EXPECT_TRUE(success)
      << "The permutation computed is different from the reference solution!";

  handle.sort_factors();
  crs_matrix_type U = handle.getU();
  crs_matrix_type L = handle.getL();
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_mdf() {
  Test::run_test_mdf<scalar_t, lno_t, size_type, device>();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)     \
  TEST_F(TestCategory,                                                  \
         sparse##_##mdf##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_mdf<SCALAR, ORDINAL, OFFSET, DEVICE>();                        \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX
