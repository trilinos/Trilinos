// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_TESTMATRIXUTILS_HPP
#define KOKKOSKERNELS_TESTMATRIXUTILS_HPP

#include <random>

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_ArithTraits.hpp"
#include "KokkosBatched_Vector.hpp"
// Make this include-able from all subdirectories

namespace Test {

template <typename scalar_t, typename lno_t, typename size_type, typename device, typename crsMat_t>
crsMat_t symmetrize(crsMat_t A) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  auto host_rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto host_entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto host_values  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  lno_t numRows     = A.numRows();
  // symmetrize as input_mat + input_mat^T, to still have a diagonally dominant
  // matrix
  typedef std::map<lno_t, scalar_t> Row;
  std::vector<Row> symRows(numRows);
  for (lno_t r = 0; r < numRows; r++) {
    auto& row = symRows[r];
    for (size_type i = host_rowmap(r); i < host_rowmap(r + 1); i++) {
      lno_t c   = host_entries(i);
      auto& col = symRows[c];
      auto it   = row.find(c);
      if (it == row.end())
        row[c] = host_values(i);
      else
        row[c] += host_values(i);
      it = col.find(r);
      if (it == col.end())
        col[r] = host_values(i);
      else
        col[r] += host_values(i);
    }
  }
  // Count entries
  Kokkos::View<size_type*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_rowmap("Rowmap", numRows + 1);
  size_t accum = 0;
  for (lno_t r = 0; r <= numRows; r++) {
    new_host_rowmap(r) = accum;
    if (r < numRows) accum += symRows[r].size();
  }
  // Allocate new entries/values
  Kokkos::View<lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_entries("Entries", accum);
  Kokkos::View<scalar_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_values("Values", accum);
  for (lno_t r = 0; r < numRows; r++) {
    auto rowIt = symRows[r].begin();
    for (size_type i = new_host_rowmap(r); i < new_host_rowmap(r + 1); i++) {
      new_host_entries(i) = rowIt->first;
      new_host_values(i)  = rowIt->second;
      rowIt++;
    }
  }
  lno_view_t new_rowmap("Rowmap", numRows + 1);
  lno_nnz_view_t new_entries("Entries", accum);
  scalar_view_t new_values("Values", accum);
  Kokkos::deep_copy(new_rowmap, new_host_rowmap);
  Kokkos::deep_copy(new_entries, new_host_entries);
  Kokkos::deep_copy(new_values, new_host_values);
  return crsMat_t("SymA", numRows, numRows, accum, new_values, new_rowmap, new_entries);
}

}  // namespace Test
#endif
