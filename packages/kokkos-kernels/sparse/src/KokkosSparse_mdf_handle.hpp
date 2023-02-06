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

#ifndef KOKKOSSPARSE_MDF_HANDLE_HPP_
#define KOKKOSSPARSE_MDF_HANDLE_HPP_

#include "KokkosSparse_SortCrs.hpp"
#include "KokkosSparse_Utils.hpp"

namespace KokkosSparse {
namespace Experimental {

template <class matrix_type>
struct MDF_handle {
  using crs_matrix_type = matrix_type;
  using execution_space = typename matrix_type::execution_space;
  using row_map_type    = typename crs_matrix_type::StaticCrsGraphType::
      row_map_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;

  ordinal_type numRows;

  // Views needed to construct L and U
  // at the end of the numerical phase.
  row_map_type row_mapL, row_mapU;
  col_ind_type entriesL, entriesU;
  values_type valuesL, valuesU;

  // Row permutation that defines
  // the MDF ordering or order of
  // elimination during the factorization.
  col_ind_type permutation, permutation_inv;

  int verbosity;

  MDF_handle(const crs_matrix_type A)
      : numRows(A.numRows()),
        permutation(col_ind_type("row permutation", A.numRows())),
        permutation_inv(col_ind_type("inverse row permutation", A.numRows())),
        verbosity(0){};

  void set_verbosity(const int verbosity_level) { verbosity = verbosity_level; }

  void allocate_data(const size_type nnzL, const size_type nnzU) {
    // Allocate L
    row_mapL = row_map_type("row map L", numRows + 1);
    entriesL = col_ind_type("entries L", nnzL);
    valuesL  = values_type("values L", nnzL);

    // Allocate U
    row_mapU = row_map_type("row map U", numRows + 1);
    entriesU = col_ind_type("entries U", nnzU);
    valuesU  = values_type("values U", nnzU);
  }

  col_ind_type get_permutation() { return permutation; }

  void sort_factors() {
    KokkosSparse::sort_crs_matrix<execution_space, row_map_type, col_ind_type,
                                  values_type>(row_mapL, entriesL, valuesL);
    KokkosSparse::sort_crs_matrix<execution_space, row_map_type, col_ind_type,
                                  values_type>(row_mapU, entriesU, valuesU);
  }

  crs_matrix_type getL() {
    return KokkosSparse::Impl::transpose_matrix<crs_matrix_type>(
        crs_matrix_type("L", numRows, numRows, entriesL.extent(0), valuesL,
                        row_mapL, entriesL));
  }

  crs_matrix_type getU() {
    return crs_matrix_type("U", numRows, numRows, entriesU.extent(0), valuesU,
                           row_mapU, entriesU);
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_MDF_HANDLE_HPP_
