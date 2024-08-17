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

/// \file KokkosSparse_spiluk.hpp
/// \brief Parallel Minimum Discarded Fill method
/// \author Luc Berger-Vergiat
/// \date March 2022
///
/// This file provides KokkosSparse::mdf_symbolic, KokkosSparse::mdf_symbolic
/// and KokkosSparse::mdf_ordering.  These functions perform a
/// local (no MPI) sparse MDF(0) on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_MDF_HPP_
#define KOKKOSSPARSE_MDF_HPP_

#include <Kokkos_UnorderedMap.hpp>
#include "KokkosSparse_mdf_handle.hpp"
#include "KokkosSparse_mdf_impl.hpp"

namespace KokkosSparse {
namespace Experimental {

template <class crs_matrix_type, class MDF_handle>
void mdf_symbolic(const crs_matrix_type& A, MDF_handle& handle) {
  using size_type = typename crs_matrix_type::size_type;

  using execution_space        = typename crs_matrix_type::execution_space;
  using team_range_policy_type = Kokkos::TeamPolicy<execution_space>;

  // Symbolic phase:
  // compute transpose of A for easy access to columns of A
  // allocate temporaries
  // allocate L and U
  size_type nnzL = 0, nnzU = 0;
  team_range_policy_type setupPolicy(A.numRows(), Kokkos::AUTO);
  KokkosSparse::Impl::MDF_count_lower<crs_matrix_type> compute_nnzL(A, handle.permutation, handle.permutation_inv);
  Kokkos::parallel_reduce(setupPolicy, compute_nnzL, nnzL);
  nnzU = A.nnz() - nnzL + A.numRows();
  handle.allocate_data(nnzL, nnzU);

  if (handle.verbosity > 0) {
    printf("MDF symbolic:  nnzL = %d, nnzU = %d\n", static_cast<int>(nnzL), static_cast<int>(nnzU));
  }

  return;
}  // mdf_symbolic

template <class view_t, class ordinal_t = size_t>
void mdf_print_joined_view(const view_t& dev_view, const char* sep,
                           ordinal_t max_count = Kokkos::ArithTraits<ordinal_t>::max()) {
  const auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev_view);

  max_count = max_count > (ordinal_t)host_view.extent(0) ? (ordinal_t)host_view.extent(0) : max_count;
  for (ordinal_t i = 0; i < max_count; ++i) {
    if (i) printf("%s", sep);
    printf("%g", static_cast<double>(host_view[i]));
  }
}

template <class crs_matrix_type, class MDF_handle>
void mdf_numeric(const crs_matrix_type& A, MDF_handle& handle) {
  using col_ind_type    = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using scalar_mag_type = typename KokkosSparse::Impl::MDF_types<crs_matrix_type>::scalar_mag_type;
  using values_mag_type = typename KokkosSparse::Impl::MDF_types<crs_matrix_type>::values_mag_type;
  using ordinal_type    = typename crs_matrix_type::ordinal_type;
  using value_mag_type  = typename values_mag_type::value_type;

  using device_type            = typename crs_matrix_type::device_type;
  using execution_space        = typename crs_matrix_type::execution_space;
  using range_policy_type      = Kokkos::RangePolicy<ordinal_type, execution_space>;
  using team_range_policy_type = Kokkos::TeamPolicy<execution_space>;

  using permutation_set_type = Kokkos::UnorderedMap<ordinal_type, void, device_type>;

  // Numerical phase:
  // loop over rows
  //   compute discarded fill of each row
  //   selected pivot based on MDF
  //   factorize pivot row of A
  const int verbosity_level = handle.verbosity;
  crs_matrix_type Atmp      = crs_matrix_type("A fill", A);
  crs_matrix_type At        = KokkosSparse::Impl::transpose_matrix<crs_matrix_type>(A);
  KokkosSparse::sort_crs_matrix<crs_matrix_type>(At);
  values_mag_type discarded_fill("discarded fill", A.numRows());
  col_ind_type deficiency("deficiency", A.numRows());
  ordinal_type update_list_len = 0;
  col_ind_type update_list("update list", A.numRows());
  col_ind_type factored("factored rows", A.numRows());
  Kokkos::deep_copy(discarded_fill, Kokkos::ArithTraits<value_mag_type>::max());
  Kokkos::deep_copy(deficiency, Kokkos::ArithTraits<ordinal_type>::max());
  permutation_set_type permutation_set(A.numRows());

  KokkosSparse::Impl::MDF_discarded_fill_norm<crs_matrix_type, true> MDF_df_norm(
      Atmp, At, 0, handle.permutation, permutation_set, discarded_fill, deficiency, verbosity_level);
  Kokkos::parallel_for("MDF: initial fill computation",
                       team_range_policy_type(Atmp.numRows(), Kokkos::AUTO, Kokkos::AUTO), MDF_df_norm);

  for (ordinal_type factorization_step = 0; factorization_step < A.numRows(); ++factorization_step) {
    if (verbosity_level > 0) {
      printf("\n\nFactorization step %d\n", static_cast<int>(factorization_step));
    }

    if (update_list_len > 0) {
      team_range_policy_type updatePolicy(update_list_len, Kokkos::AUTO, Kokkos::AUTO);
      KokkosSparse::Impl::MDF_discarded_fill_norm<crs_matrix_type, false> MDF_update_df_norm(
          Atmp, At, factorization_step, handle.permutation, permutation_set, discarded_fill, deficiency,
          verbosity_level, update_list);
      Kokkos::parallel_for("MDF: updating fill norms", updatePolicy, MDF_update_df_norm);
    }

    if (verbosity_level > 1) {
      if constexpr (std::is_arithmetic_v<scalar_mag_type>) {
        printf("  discarded_fill = {");
        mdf_print_joined_view(discarded_fill, ", ");
        printf("}\n");
      }
      printf("  deficiency = {");
      mdf_print_joined_view(deficiency, ", ");
      printf("}\n");
    }

    ordinal_type selected_row_idx = 0;
    {
      range_policy_type stepPolicy(factorization_step, Atmp.numRows());
      KokkosSparse::Impl::MDF_select_row<crs_matrix_type> MDF_row_selector(
          factorization_step, discarded_fill, deficiency, Atmp.graph.row_map, handle.permutation);
      Kokkos::parallel_reduce("MDF: select pivot", stepPolicy, MDF_row_selector, selected_row_idx);
    }

    ordinal_type selected_row_len = 0;
    {
      // vector overloads required for scans to use vector parallel not yet
      // provided by kokkos (https://github.com/kokkos/kokkos/issues/6259)
      team_range_policy_type updateListPolicy(1, Kokkos::AUTO);
      KokkosSparse::Impl::MDF_compute_list_length<crs_matrix_type> updateList(
          Atmp, At, handle.row_mapL, handle.entriesL, handle.valuesL, handle.row_mapU, handle.entriesU, handle.valuesU,
          handle.permutation, handle.permutation_inv, permutation_set, discarded_fill, factored, selected_row_idx,
          factorization_step, update_list, verbosity_level);
      update_list_len = 0;
      Kokkos::parallel_reduce("MDF: compute update list", updateListPolicy, updateList, update_list_len,
                              selected_row_len);
    }

    if (verbosity_level > 1) {
      printf("  updateList = {");
      mdf_print_joined_view(update_list, ", ", update_list_len);
      printf("}\n  permutation = {");
      mdf_print_joined_view(handle.permutation, ", ");
      printf("}\n  permutation_inv = {");
      mdf_print_joined_view(handle.permutation_inv, ", ");
      printf("}\n");
    }
    if (verbosity_level > 0) {
      printf(
          "  Selected row idx %d with length %d. Requires update of %d fill "
          "norms.\n",
          static_cast<int>(selected_row_idx), static_cast<int>(selected_row_len), static_cast<int>(update_list_len));
    }

    // If this was the last row no need to update A and At!
    if (factorization_step < A.numRows() - 1) {
      team_range_policy_type factorizePolicy(selected_row_len, Kokkos::AUTO, Kokkos::AUTO);
      KokkosSparse::Impl::MDF_factorize_row<crs_matrix_type> factorize_row(
          Atmp, At, handle.row_mapL, handle.entriesL, handle.valuesL, handle.row_mapU, handle.entriesU, handle.valuesU,
          handle.permutation, handle.permutation_inv, permutation_set, discarded_fill, factored, selected_row_idx,
          factorization_step, update_list, verbosity_level);
      Kokkos::parallel_for("MDF: factorize row", factorizePolicy, factorize_row);
    }
  }  // Loop over factorization steps

  KokkosSparse::Impl::MDF_reindex_matrix<col_ind_type> reindex_U(handle.permutation_inv, handle.entriesU);
  Kokkos::parallel_for("MDF: re-index U", range_policy_type(0, handle.entriesU.extent(0)), reindex_U);

  KokkosSparse::Impl::MDF_reindex_matrix<col_ind_type> reindex_L(handle.permutation_inv, handle.entriesL);
  Kokkos::parallel_for("MDF: re-index L", range_policy_type(0, handle.entriesL.extent(0)), reindex_L);

  handle.L = KokkosSparse::Impl::transpose_matrix<crs_matrix_type>(handle.L);

  return;
}  // mdf_numeric

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_MDF_HPP_
