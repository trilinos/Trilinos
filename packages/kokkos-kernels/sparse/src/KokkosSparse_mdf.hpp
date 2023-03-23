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

#include "KokkosSparse_mdf_handle.hpp"
#include "KokkosSparse_mdf_impl.hpp"

namespace KokkosSparse {
namespace Experimental {

template <class crs_matrix_type, class MDF_handle>
void mdf_symbolic(crs_matrix_type& A, MDF_handle& handle) {
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;

  using execution_space   = typename crs_matrix_type::execution_space;
  using range_policy_type = Kokkos::RangePolicy<ordinal_type, execution_space>;

  // Symbolic phase:
  // compute transpose of A for easy access to columns of A
  // allocate temporaries
  // allocate L and U
  size_type nnzL = 0, nnzU = 0;
  range_policy_type setupPolicy(0, A.numRows());
  KokkosSparse::Impl::MDF_count_lower<crs_matrix_type> compute_nnzL(
      A, handle.permutation, handle.permutation_inv);
  Kokkos::parallel_reduce(range_policy_type(0, A.numRows()), compute_nnzL,
                          nnzL);
  nnzU = A.nnz() - nnzL + A.numRows();
  handle.allocate_data(nnzL, nnzU);

  if (handle.verbosity > 0) {
    printf("MDF symbolic:  nnzL = %d, nnzU = %d\n", static_cast<int>(nnzL),
           static_cast<int>(nnzU));
  }

  return;
}  // mdf_symbolic

template <class crs_matrix_type, class MDF_handle>
void mdf_numeric(crs_matrix_type& A, MDF_handle& handle) {
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using value_type   = typename crs_matrix_type::value_type;

  using execution_space   = typename crs_matrix_type::execution_space;
  using range_policy_type = Kokkos::RangePolicy<ordinal_type, execution_space>;

  // Numerical phase:
  // loop over rows
  //   compute discarded fill of each row
  //   selected pivot based on MDF
  //   factorize pivot row of A
  const int verbosity_level = handle.verbosity;
  crs_matrix_type Atmp      = crs_matrix_type("A fill", A);
  crs_matrix_type At = KokkosSparse::Impl::transpose_matrix<crs_matrix_type>(A);
  KokkosSparse::sort_crs_matrix<crs_matrix_type>(At);
  values_type discarded_fill("discarded fill", A.numRows());
  col_ind_type deficiency("deficiency", A.numRows());
  col_ind_type update_list_length("update list length", 1);
  typename col_ind_type::HostMirror update_list_length_host =
      Kokkos::create_mirror_view(update_list_length);
  col_ind_type update_list("update list", A.numRows());
  col_ind_type factored("factored rows", A.numRows());
  Kokkos::deep_copy(discarded_fill, Kokkos::ArithTraits<value_type>::max());
  Kokkos::deep_copy(deficiency, Kokkos::ArithTraits<ordinal_type>::max());

  KokkosSparse::Impl::MDF_discarded_fill_norm<crs_matrix_type> MDF_df_norm(
      Atmp, At, 0, handle.permutation, discarded_fill, deficiency,
      verbosity_level);
  Kokkos::parallel_for("MDF: initial fill computation",
                       range_policy_type(0, Atmp.numRows()), MDF_df_norm);

  for (ordinal_type factorization_step = 0; factorization_step < A.numRows();
       ++factorization_step) {
    if (verbosity_level > 0) {
      printf("\n\nFactorization step %d\n\n",
             static_cast<int>(factorization_step));
    }

    Kokkos::deep_copy(update_list_length_host, update_list_length);
    range_policy_type updatePolicy(0, update_list_length_host(0));
    KokkosSparse::Impl::MDF_selective_discarded_fill_norm<crs_matrix_type>
        MDF_update_df_norm(Atmp, At, factorization_step, handle.permutation,
                           update_list, discarded_fill, deficiency,
                           verbosity_level);
    Kokkos::parallel_for("MDF: updating fill norms", updatePolicy,
                         MDF_update_df_norm);

    range_policy_type stepPolicy(factorization_step, Atmp.numRows());
    ordinal_type selected_row_idx = 0;
    KokkosSparse::Impl::MDF_select_row<crs_matrix_type> MDF_row_selector(
        factorization_step, discarded_fill, deficiency, Atmp.graph.row_map,
        handle.permutation);
    Kokkos::parallel_reduce("MDF: select pivot", stepPolicy, MDF_row_selector,
                            selected_row_idx);

    KokkosSparse::Impl::MDF_compute_list_length<crs_matrix_type>
        compute_list_length(selected_row_idx, Atmp, At, handle.permutation,
                            factored, update_list_length, update_list);
    Kokkos::parallel_for("MDF: compute update list", range_policy_type(0, 1),
                         compute_list_length);

    KokkosSparse::Impl::MDF_factorize_row<crs_matrix_type> factorize_row(
        Atmp, At, handle.row_mapL, handle.entriesL, handle.valuesL,
        handle.row_mapU, handle.entriesU, handle.valuesU, handle.permutation,
        handle.permutation_inv, discarded_fill, factored, selected_row_idx,
        factorization_step, verbosity_level);
    Kokkos::parallel_for("MDF: factorize row", range_policy_type(0, 1),
                         factorize_row);

    if (verbosity_level > 0) {
      printf("\n");
    }
  }  // Loop over factorization steps

  KokkosSparse::Impl::MDF_reindex_matrix<col_ind_type> reindex_U(
      handle.permutation_inv, handle.entriesU);
  Kokkos::parallel_for("MDF: re-index U",
                       range_policy_type(0, handle.entriesU.extent(0)),
                       reindex_U);

  KokkosSparse::Impl::MDF_reindex_matrix<col_ind_type> reindex_L(
      handle.permutation_inv, handle.entriesL);
  Kokkos::parallel_for("MDF: re-index L",
                       range_policy_type(0, handle.entriesL.extent(0)),
                       reindex_L);

  handle.L = KokkosSparse::Impl::transpose_matrix<crs_matrix_type>(handle.L);

  return;
}  // mdf_numeric

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_MDF_HPP_
