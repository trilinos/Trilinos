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

#include <string>
#include <stdexcept>

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_par_ilut.hpp"

#include <gtest/gtest.h>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

namespace Test {

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
std::vector<std::vector<scalar_t>> decompress_matrix(
    Kokkos::View<size_type*, device>& row_map,
    Kokkos::View<lno_t*, device>& entries,
    Kokkos::View<scalar_t*, device>& values) {
  const size_type nrows = row_map.size() - 1;
  std::vector<std::vector<scalar_t>> result;
  result.resize(nrows);
  for (auto& row : result) {
    row.resize(nrows, 0.0);
  }

  auto hrow_map = Kokkos::create_mirror_view(row_map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);
  Kokkos::deep_copy(hrow_map, row_map);
  Kokkos::deep_copy(hentries, entries);
  Kokkos::deep_copy(hvalues, values);

  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    const size_type row_nnz_begin = hrow_map(row_idx);
    const size_type row_nnz_end   = hrow_map(row_idx + 1);
    for (size_type row_nnz = row_nnz_begin; row_nnz < row_nnz_end; ++row_nnz) {
      const lno_t col_idx      = hentries(row_nnz);
      const scalar_t value     = hvalues(row_nnz);
      result[row_idx][col_idx] = value;
    }
  }

  return result;
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void check_matrix(const std::string& name,
                  Kokkos::View<size_type*, device>& row_map,
                  Kokkos::View<lno_t*, device>& entries,
                  Kokkos::View<scalar_t*, device>& values,
                  const std::vector<std::vector<scalar_t>>& expected) {
  const auto decompressed_mtx = decompress_matrix(row_map, entries, values);

  const size_type nrows = row_map.size() - 1;
  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    for (size_type col_idx = 0; col_idx < nrows; ++col_idx) {
      EXPECT_NEAR(expected[row_idx][col_idx],
                  decompressed_mtx[row_idx][col_idx], 0.01)
          << "Failed check is: " << name << "[" << row_idx << "][" << col_idx
          << "]";
    }
  }
}

template <typename scalar_t>
void print_matrix(const std::vector<std::vector<scalar_t>>& matrix) {
  for (const auto& row : matrix) {
    for (const auto& item : row) {
      std::printf("%.2f ", item);
    }
    std::cout << std::endl;
  }
}

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void run_test_par_ilut() {
  using RowMapType   = Kokkos::View<size_type*, device>;
  using EntriesType  = Kokkos::View<lno_t*, device>;
  using ValuesType   = Kokkos::View<scalar_t*, device>;
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, typename device::execution_space,
      typename device::memory_space, typename device::memory_space>;

  // Simple test fixture A
  std::vector<std::vector<scalar_t>> A = {{1., 6., 4., 7.},
                                          {2., -5., 0., 8.},
                                          {0.5, -3., 6., 0.},
                                          {0.2, -0.5, -9., 0.}};

  const scalar_t ZERO = scalar_t(0);

  const size_type nrows = A.size();

  // Count A nnz's
  size_type nnz = 0;
  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    for (size_type col_idx = 0; col_idx < nrows; ++col_idx) {
      if (A[row_idx][col_idx] != ZERO) {
        ++nnz;
      }
    }
  }

  // Allocate device CRS views for A
  RowMapType row_map("row_map", nrows + 1);
  EntriesType entries("entries", nnz);
  ValuesType values("values", nnz);

  // Create host mirror views for CRS A
  auto hrow_map = Kokkos::create_mirror_view(row_map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);

  // Compress A into CRS (host views)
  size_type curr_nnz = 0;
  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    for (size_type col_idx = 0; col_idx < nrows; ++col_idx) {
      if (A[row_idx][col_idx] != ZERO) {
        hentries(curr_nnz) = col_idx;
        hvalues(curr_nnz)  = A[row_idx][col_idx];
        ++curr_nnz;
      }
      hrow_map(row_idx + 1) = curr_nnz;
    }
  }

  // Copy host A CRS views to device A CRS views
  Kokkos::deep_copy(row_map, hrow_map);
  Kokkos::deep_copy(entries, hentries);
  Kokkos::deep_copy(values, hvalues);

  // Make kernel handle
  KernelHandle kh;

  kh.create_par_ilut_handle(nrows);

  auto par_ilut_handle = kh.get_par_ilut_handle();

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", nrows + 1);
  RowMapType U_row_map("U_row_map", nrows + 1);

  // Initial L/U approximations for A
  par_ilut_symbolic(&kh, row_map, entries, L_row_map, U_row_map);

  const size_type nnzL = par_ilut_handle->get_nnzL();
  const size_type nnzU = par_ilut_handle->get_nnzU();

  EXPECT_EQ(nnzL, 10);
  EXPECT_EQ(nnzU, 8);

  EntriesType L_entries("L_entries", nnzL);
  ValuesType L_values("L_values", nnzL);
  EntriesType U_entries("U_entries", nnzU);
  ValuesType U_values("U_values", nnzU);

  par_ilut_numeric(&kh, row_map, entries, values, L_row_map, L_entries,
                   L_values, U_row_map, U_entries, U_values,
#ifdef KOKKOS_ENABLE_SERIAL
                   true /*deterministic*/
#else
                   false /*cannot ask for determinism*/
#endif
  );

  // Use this to check LU
  // std::vector<std::vector<scalar_t> > expected_LU = {
  //   {1.0, 6.0, 4.0, 7.0},
  //   {2.0, 7.0, 8.0, 22.0},
  //   {0.5, 18.0, 8.0, -20.5},
  //   {0.2, 3.7, -53.2, -1.60}
  // };

  // check_matrix("LU numeric", L_row_map, L_entries, L_values, expected_LU);

  // Use these fixtures to test add_candidates
  // std::vector<std::vector<scalar_t> > expected_L_candidates = {
  //   {1., 0., 0., 0.},
  //   {2., 1., 0., 0.},
  //   {0.50, -3., 1., 0.},
  //   {0.20, -0.50, -9., 1.}
  // };

  // check_matrix("L numeric", L_row_map, L_entries, L_values,
  // expected_L_candidates);

  // std::vector<std::vector<scalar_t> > expected_U_candidates = {
  //   {1., 6., 4., 7.},
  //   {0., -5., -8., 8.},
  //   {0., 0., 6., 20.50},
  //   {0., 0., 0., 1.}
  // };

  // check_matrix("U numeric", U_row_map, U_entries, U_values,
  // expected_U_candidates);

  // Use these fixtures to test compute_l_u_factors
  // std::vector<std::vector<scalar_t> > expected_L_candidates = {
  //   {1., 0., 0., 0.},
  //   {2., 1., 0., 0.},
  //   {0.50, 0.35, 1., 0.},
  //   {0.20, 0.10, -1.32, 1.}
  // };

  // check_matrix("L numeric", L_row_map, L_entries, L_values,
  // expected_L_candidates);

  // std::vector<std::vector<scalar_t> > expected_U_candidates = {
  //   {1., 6., 4., 7.},
  //   {0., -17., -8., -6.},
  //   {0., 0., 6.82, -1.38},
  //   {0., 0., 0., -2.62}
  // };

  // check_matrix("U numeric", U_row_map, U_entries, U_values,
  // expected_U_candidates);

  // Serial is required for deterministic mode and the checks below cannot
  // reliably pass without determinism.
#ifdef KOKKOS_ENABLE_SERIAL

  // Use these fixtures to test full numeric
  std::vector<std::vector<scalar_t>> expected_L_candidates = {
      {1., 0., 0., 0.},
      {2., 1., 0., 0.},
      {0.50, 0.35, 1., 0.},
      {0., 0., -1.32, 1.}};

  check_matrix("L numeric", L_row_map, L_entries, L_values,
               expected_L_candidates);

  std::vector<std::vector<scalar_t>> expected_U_candidates = {
      {1., 6., 4., 7.},
      {0., -17., -8., -6.},
      {0., 0., 6.82, 0.},
      {0., 0., 0., 0.}  // [3] = 0 for full alg, -2.62 for post-threshold only
  };

  check_matrix("U numeric", U_row_map, U_entries, U_values,
               expected_U_candidates);

  // Checking

  kh.destroy_par_ilut_handle();
#endif
}

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_par_ilut() {
  Test::run_test_par_ilut<scalar_t, lno_t, size_type, device>();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)          \
  TEST_F(TestCategory,                                                       \
         sparse##_##par_ilut##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_par_ilut<SCALAR, ORDINAL, OFFSET, DEVICE>();                        \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX
