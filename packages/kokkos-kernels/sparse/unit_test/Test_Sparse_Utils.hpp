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

#ifndef TEST_SPARSE_UTILS_HPP
#define TEST_SPARSE_UTILS_HPP

#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_SortCrs.hpp"

namespace Test {

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"), crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector_mv(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"), crsMat.numRows(), x_vector.extent(1));
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat_actual, crsMat_t output_mat_reference) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows_actual    = output_mat_actual.numRows();
  size_t ncols_actual    = output_mat_actual.numCols();
  size_t nentries_actual = output_mat_actual.graph.entries.extent(0);
  size_t nvals_actual    = output_mat_actual.values.extent(0);

  size_t nrows_reference    = output_mat_reference.numRows();
  size_t ncols_reference    = output_mat_reference.numCols();
  size_t nentries_reference = output_mat_reference.graph.entries.extent(0);
  size_t nvals_reference    = output_mat_reference.values.extent(0);

  if (nrows_actual != nrows_reference || ncols_actual != ncols_reference) {
    std::cout << "dimensions (actual):" << nrows_actual << 'x' << ncols_actual
              << ", dimensions (reference): " << nrows_reference << 'x' << ncols_reference << '\n';
    return false;
  }
  if (nentries_actual != nentries_reference) {
    std::cout << "nentries_actual:" << nentries_actual << " nentries_reference:" << nentries_reference << std::endl;
    return false;
  }
  if (nvals_actual != nvals_reference) {
    std::cout << "nvals_actual:" << nvals_actual << " nvals_reference:" << nvals_reference << std::endl;
    return false;
  }

  bool is_identical = true;
  // Special case: a matrix with 0 rows can have a rowmap of length 0 or 1.
  // Treat these as equivalent.
  bool zero_row_equivalent = false;
  if (nrows_reference == 0) {
    auto rm1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_actual.graph.row_map);
    auto rm2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_reference.graph.row_map);
    if (rm1.extent_int(0) == 0 && rm2.extent_int(0) == 1) {
      // Make sure the one element of rm2 is 0
      zero_row_equivalent = !rm2(0);
    } else if (rm1.extent_int(0) == 1 && rm2.extent_int(0) == 0) {
      // Make sure the one element of rm1 is 0
      zero_row_equivalent = !rm1(0);
    }
  }
  if (!zero_row_equivalent) {
    is_identical =
        KokkosKernels::Impl::kk_is_identical_view<typename graph_t::row_map_type, typename graph_t::row_map_type,
                                                  typename lno_view_t::value_type, typename device::execution_space>(
            output_mat_actual.graph.row_map, output_mat_reference.graph.row_map, 0);
  }

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    std::cout << "Actual rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.row_map, true);
    std::cout << "Correct rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.row_map, true);
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
                                                typename device::execution_space>(
          output_mat_actual.graph.entries, output_mat_reference.graph.entries, 0);

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.entries);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.entries);
    return false;
  }

  typedef typename Kokkos::ArithTraits<typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type, float>::value ? 3.7e-3 : 1e-7;

  is_identical = KokkosKernels::Impl::kk_is_relatively_identical_view<scalar_view_t, scalar_view_t, eps_type,
                                                                      typename device::execution_space>(
      output_mat_actual.values, output_mat_reference.values, eps);

  if (!is_identical) {
    std::cout << "values are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.values);

    return false;
  }
  return true;
}
}  // namespace Test

#endif  // TEST_SPARSE_UTILS_HPP
