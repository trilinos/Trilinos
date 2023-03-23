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
#ifndef _KOKKOSSPARSE_SORTCRS_HPP
#define _KOKKOSSPARSE_SORTCRS_HPP

#include "Kokkos_Core.hpp"
#include "KokkosKernels_Sorting.hpp"

namespace KokkosSparse {

// ----------------------------------
// BSR matrix/graph sorting utilities
// ----------------------------------

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t,
          typename lno_t = typename entries_t::non_const_value_type>
void sort_bsr_matrix(const lno_t blockdim, const rowmap_t& rowmap,
                     const entries_t& entries, const values_t& values);

template <typename bsrMat_t>
void sort_bsr_matrix(const bsrMat_t& A);

// ----------------------------------
// CRS matrix/graph sorting utilities
// ----------------------------------

// The sort_crs* functions sort the adjacent column list for each row into
// ascending order.

template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries,
                     const values_t& values);

template <typename crsMat_t>
void sort_crs_matrix(const crsMat_t& A);

template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries);

template <typename crsGraph_t>
void sort_crs_graph(const crsGraph_t& G);

// sort_and_merge_matrix produces a new matrix which is equivalent to A but is
// sorted and has no duplicated entries: each (i, j) is unique. Values for
// duplicated entries are summed.
template <typename crsMat_t>
crsMat_t sort_and_merge_matrix(const crsMat_t& A);

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(const crsGraph_t& G);

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in,
                          const entries_t& entries_in, rowmap_t& rowmap_out,
                          entries_t& entries_out);

namespace Impl {

template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t>
struct SortCrsMatrixFunctor {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  using scalar_t  = typename values_t::non_const_value_type;
  using team_mem  = typename Kokkos::TeamPolicy<execution_space>::member_type;
  // The functor owns memory for entriesAux, so it can't have
  // MemoryTraits<Unmanaged>
  using entries_managed_t = Kokkos::View<typename entries_t::data_type,
                                         typename entries_t::device_type>;
  using values_managed_t  = Kokkos::View<typename values_t::data_type,
                                        typename values_t::device_type>;

  SortCrsMatrixFunctor(bool usingRangePol, const rowmap_t& rowmap_,
                       const entries_t& entries_, const values_t& values_)
      : rowmap(rowmap_), entries(entries_), values(values_) {
    if (usingRangePol) {
      entriesAux = entries_managed_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"),
          entries.extent(0));
      valuesAux = values_managed_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values aux"),
          values.extent(0));
    }
    // otherwise, aux arrays won't be allocated (sorting in place)
  }

  KOKKOS_INLINE_FUNCTION void operator()(const lno_t i) const {
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    // Radix sort requires unsigned keys for comparison
    using unsigned_lno_t = typename std::make_unsigned<lno_t>::type;
    KokkosKernels::SerialRadixSort2<lno_t, unsigned_lno_t, scalar_t>(
        (unsigned_lno_t*)entries.data() + rowStart,
        (unsigned_lno_t*)entriesAux.data() + rowStart, values.data() + rowStart,
        valuesAux.data() + rowStart, rowNum);
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_mem t) const {
    size_type i        = t.league_rank();
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    KokkosKernels::TeamBitonicSort2<lno_t, lno_t, scalar_t, team_mem>(
        entries.data() + rowStart, values.data() + rowStart, rowNum, t);
  }

  rowmap_t rowmap;
  entries_t entries;
  entries_managed_t entriesAux;
  values_t values;
  values_managed_t valuesAux;
};

template <typename execution_space, typename rowmap_t, typename entries_t>
struct SortCrsGraphFunctor {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  using team_mem  = typename Kokkos::TeamPolicy<execution_space>::member_type;
  // The functor owns memory for entriesAux, so it can't have
  // MemoryTraits<Unmanaged>
  using entries_managed_t = Kokkos::View<typename entries_t::data_type,
                                         typename entries_t::device_type>;

  SortCrsGraphFunctor(bool usingRangePol, const rowmap_t& rowmap_,
                      const entries_t& entries_)
      : rowmap(rowmap_), entries(entries_) {
    if (usingRangePol) {
      entriesAux = entries_managed_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"),
          entries.extent(0));
    }
    // otherwise, aux arrays won't be allocated (sorting in place)
  }

  KOKKOS_INLINE_FUNCTION void operator()(const lno_t i) const {
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    // Radix sort requires unsigned keys for comparison
    using unsigned_lno_t = typename std::make_unsigned<lno_t>::type;
    KokkosKernels::SerialRadixSort<lno_t, unsigned_lno_t>(
        (unsigned_lno_t*)entries.data() + rowStart,
        (unsigned_lno_t*)entriesAux.data() + rowStart, rowNum);
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_mem t) const {
    size_type i        = t.league_rank();
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    KokkosKernels::TeamBitonicSort<lno_t, lno_t, team_mem>(
        entries.data() + rowStart, rowNum, t);
  }

  rowmap_t rowmap;
  entries_t entries;
  entries_managed_t entriesAux;
};

template <typename rowmap_t, typename entries_t>
struct MergedRowmapFunctor {
  using size_type  = typename rowmap_t::non_const_value_type;
  using lno_t      = typename entries_t::non_const_value_type;
  using c_rowmap_t = typename rowmap_t::const_type;

  // Precondition: entries are sorted within each row
  MergedRowmapFunctor(const rowmap_t& mergedCounts_, const c_rowmap_t& rowmap_,
                      const entries_t& entries_)
      : mergedCounts(mergedCounts_), rowmap(rowmap_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION void operator()(lno_t row, size_type& lnewNNZ) const {
    size_type rowBegin = rowmap(row);
    size_type rowEnd   = rowmap(row + 1);
    if (rowEnd == rowBegin) {
      // Row was empty to begin with
      mergedCounts(row) = 0;
      return;
    }
    // Otherwise, the first entry in the row exists
    lno_t uniqueEntries = 1;
    for (size_type j = rowBegin + 1; j < rowEnd; j++) {
      if (entries(j - 1) != entries(j)) uniqueEntries++;
    }
    mergedCounts(row) = uniqueEntries;
    lnewNNZ += uniqueEntries;
    if (row == lno_t((rowmap.extent(0) - 1) - 1)) mergedCounts(row + 1) = 0;
  }

  rowmap_t mergedCounts;
  c_rowmap_t rowmap;
  entries_t entries;
};

template <typename rowmap_t, typename entries_t, typename values_t>
struct MatrixMergedEntriesFunctor {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  using scalar_t  = typename values_t::non_const_value_type;

  // Precondition: entries are sorted within each row
  MatrixMergedEntriesFunctor(const rowmap_t& rowmap_, const entries_t& entries_,
                             const values_t& values_,
                             const rowmap_t& mergedRowmap_,
                             const entries_t& mergedEntries_,
                             const values_t& mergedValues_)
      : rowmap(rowmap_),
        entries(entries_),
        values(values_),
        mergedRowmap(mergedRowmap_),
        mergedEntries(mergedEntries_),
        mergedValues(mergedValues_) {}

  KOKKOS_INLINE_FUNCTION void operator()(lno_t row) const {
    size_type rowBegin = rowmap(row);
    size_type rowEnd   = rowmap(row + 1);
    if (rowEnd == rowBegin) {
      // Row was empty to begin with, nothing to do
      return;
    }
    // Otherwise, accumulate the value for each column
    scalar_t accumVal   = values(rowBegin);
    lno_t accumCol      = entries(rowBegin);
    size_type insertPos = mergedRowmap(row);
    for (size_type j = rowBegin + 1; j < rowEnd; j++) {
      if (accumCol == entries(j)) {
        // accumulate
        accumVal += values(j);
      } else {
        // write out and reset
        mergedValues(insertPos)  = accumVal;
        mergedEntries(insertPos) = accumCol;
        insertPos++;
        accumVal = values(j);
        accumCol = entries(j);
      }
    }
    // always left with the last unique entry
    mergedValues(insertPos)  = accumVal;
    mergedEntries(insertPos) = accumCol;
  }

  rowmap_t rowmap;
  entries_t entries;
  values_t values;
  rowmap_t mergedRowmap;
  entries_t mergedEntries;
  values_t mergedValues;
};

template <typename rowmap_t, typename entries_t>
struct GraphMergedEntriesFunctor {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;

  // Precondition: entries are sorted within each row
  GraphMergedEntriesFunctor(const rowmap_t& rowmap_, const entries_t& entries_,
                            const rowmap_t& mergedRowmap_,
                            const entries_t& mergedEntries_)
      : rowmap(rowmap_),
        entries(entries_),
        mergedRowmap(mergedRowmap_),
        mergedEntries(mergedEntries_) {}

  KOKKOS_INLINE_FUNCTION void operator()(lno_t row) const {
    size_type rowBegin = rowmap(row);
    size_type rowEnd   = rowmap(row + 1);
    if (rowEnd == rowBegin) {
      // Row was empty to begin with, nothing to do
      return;
    }
    // Otherwise, accumulate the value for each column
    lno_t accumCol      = entries(rowBegin);
    size_type insertPos = mergedRowmap(row);
    for (size_type j = rowBegin + 1; j < rowEnd; j++) {
      if (accumCol != entries(j)) {
        // write out and reset
        mergedEntries(insertPos) = accumCol;
        insertPos++;
        accumCol = entries(j);
      }
    }
    // always left with the last unique entry
    mergedEntries(insertPos) = accumCol;
  }

  rowmap_t rowmap;
  entries_t entries;
  rowmap_t mergedRowmap;
  entries_t mergedEntries;
};

template <typename T>
KOKKOS_INLINE_FUNCTION void kk_swap(T& a, T& b) {
  T t = a;
  a   = b;
  b   = t;
}

template <typename row_map_type, typename entries_type, typename values_type>
struct sort_bsr_functor {
  using lno_t = typename entries_type::non_const_value_type;

  row_map_type rowmap;
  entries_type entries;
  values_type values;
  const lno_t blocksize;

  sort_bsr_functor(row_map_type rowmap_, entries_type entries_,
                   values_type values_, const lno_t blocksize_)
      : rowmap(rowmap_),
        entries(entries_),
        values(values_),
        blocksize(blocksize_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const lno_t i) const {
    const lno_t rowStart = rowmap(i);
    const lno_t rowSize  = rowmap(i + 1) - rowStart;
    auto* e              = entries.data() + rowStart;
    auto* v              = values.data() + rowStart * blocksize;
    bool done            = false;
    while (!done) {
      done = true;
      for (lno_t j = 1; j < rowSize; ++j) {
        const lno_t jp = j - 1;
        if (e[jp] <= e[j]) continue;
        Impl::kk_swap(e[jp], e[j]);
        auto const vb  = v + j * blocksize;
        auto const vbp = v + jp * blocksize;
        for (lno_t k = 0; k < blocksize;
             ++k)  // std::swap_ranges(vb, vb + blocksize, vbp);
          Impl::kk_swap(vb[k], vbp[k]);
        done = false;
      }
    }
  }
};

}  // namespace Impl

// Sort a CRS matrix: within each row, sort entries ascending by column.
// At the same time, permute the values.
template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries,
                     const values_t& values) {
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  bool useRadix = !KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  Impl::SortCrsMatrixFunctor<execution_space, rowmap_t, entries_t, values_t>
      funct(useRadix, rowmap, entries, values);
  if (useRadix) {
    Kokkos::parallel_for("sort_crs_matrix",
                         Kokkos::RangePolicy<execution_space>(0, numRows),
                         funct);
  } else {
    // Try to get teamsize to be largest power of 2 not greater than avg entries
    // per row
    // TODO (probably important for performnce): add thread-level sort also, and
    // use that for small avg degree. But this works for now.
    lno_t idealTeamSize = 1;
    lno_t avgDeg        = (entries.extent(0) + numRows - 1) / numRows;
    while (idealTeamSize < avgDeg / 2) {
      idealTeamSize *= 2;
    }
    team_pol temp(numRows, 1);
    lno_t maxTeamSize = temp.team_size_max(funct, Kokkos::ParallelForTag());
    lno_t teamSize    = std::min(idealTeamSize, maxTeamSize);
    Kokkos::parallel_for("sort_crs_matrix", team_pol(numRows, teamSize), funct);
  }
}

template <typename crsMat_t>
void sort_crs_matrix(const crsMat_t& A) {
  // Note: rowmap_t has const values, but that's OK as sorting doesn't modify it
  using rowmap_t   = typename crsMat_t::row_map_type;
  using entries_t  = typename crsMat_t::index_type::non_const_type;
  using values_t   = typename crsMat_t::values_type::non_const_type;
  using exec_space = typename crsMat_t::execution_space;
  // NOTE: the rowmap of a StaticCrsGraph is const-valued, but the
  // entries and CrsMatrix values are non-const (so sorting them directly
  // is allowed)
  sort_crs_matrix<exec_space, rowmap_t, entries_t, values_t>(
      A.graph.row_map, A.graph.entries, A.values);
}

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t, typename lno_t>
void sort_bsr_matrix(const lno_t blockdim, const rowmap_t& rowmap,
                     const entries_t& entries, const values_t& values) {
  // TODO: this is O(N^2) mock for debugging - do regular implementation based
  // on Radix/Bitonic sort (like CSR) IDEA: maybe we need only one general
  // Radix2/Bitonic2 and CSR sorting may call it with blockSize=1 ?
  lno_t numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  const lno_t blocksize = blockdim * blockdim;

  assert(values.extent(0) == entries.extent(0) * blocksize);
  Impl::sort_bsr_functor<rowmap_t, entries_t, values_t> bsr_sorter(
      rowmap, entries, values, blocksize);
  Kokkos::parallel_for("sort_bsr_matrix",
                       Kokkos::RangePolicy<execution_space>(0, numRows),
                       bsr_sorter);
}

// Sort a BSR matrix (like CRS but single values are replaced with contignous
// blocks)
template <typename bsrMat_t>
void sort_bsr_matrix(const bsrMat_t& A) {
  // NOTE: unlike rowmap, entries and values are non-const, so we can sort them
  // directly
  sort_bsr_matrix<typename bsrMat_t::execution_space,
                  typename bsrMat_t::row_map_type,
                  typename bsrMat_t::index_type::non_const_type,
                  typename bsrMat_t::values_type::non_const_type>(
      A.blockDim(), A.graph.row_map, A.graph.entries, A.values);
}

// Sort a CRS graph: within each row, sort entries ascending by column.
template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries) {
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  bool useRadix = !KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  Impl::SortCrsGraphFunctor<execution_space, rowmap_t, entries_t> funct(
      useRadix, rowmap, entries);
  if (useRadix) {
    Kokkos::parallel_for("sort_crs_graph",
                         Kokkos::RangePolicy<execution_space>(0, numRows),
                         funct);
  } else {
    // Try to get teamsize to be largest power of 2 less than or equal to
    // half the entries per row. 0.5 * #entries is bitonic's parallelism within
    // a row.
    // TODO (probably important for performnce): add thread-level sort also, and
    // use that for small avg degree. But this works for now.
    lno_t idealTeamSize = 1;
    lno_t avgDeg        = (entries.extent(0) + numRows - 1) / numRows;
    while (idealTeamSize < avgDeg / 2) {
      idealTeamSize *= 2;
    }
    team_pol temp(numRows, 1);
    lno_t maxTeamSize = temp.team_size_max(funct, Kokkos::ParallelForTag());
    lno_t teamSize    = std::min(idealTeamSize, maxTeamSize);
    Kokkos::parallel_for("sort_crs_graph", team_pol(numRows, teamSize), funct);
  }
}

template <typename crsGraph_t>
void sort_crs_graph(const crsGraph_t& G) {
  static_assert(
      !std::is_const<typename crsGraph_t::entries_type::value_type>::value,
      "sort_crs_graph requires StaticCrsGraph entries to be non-const.");
  sort_crs_graph<typename crsGraph_t::execution_space,
                 typename crsGraph_t::row_map_type,
                 typename crsGraph_t::entries_type>(G.row_map, G.entries);
}

// Sort the rows of matrix, and merge duplicate entries.
template <typename crsMat_t>
crsMat_t sort_and_merge_matrix(const crsMat_t& A) {
  using c_rowmap_t = typename crsMat_t::row_map_type;
  using rowmap_t   = typename crsMat_t::row_map_type::non_const_type;
  using entries_t  = typename crsMat_t::index_type::non_const_type;
  using values_t   = typename crsMat_t::values_type::non_const_type;
  using size_type  = typename rowmap_t::non_const_value_type;
  using exec_space = typename crsMat_t::execution_space;
  using range_t    = Kokkos::RangePolicy<exec_space>;
  sort_crs_matrix(A);
  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  rowmap_t mergedRowmap(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "SortedMerged rowmap"),
      A.numRows() + 1);
  size_type numCompressedEntries = 0;
  Kokkos::parallel_reduce(range_t(0, A.numRows()),
                          Impl::MergedRowmapFunctor<rowmap_t, entries_t>(
                              mergedRowmap, A.graph.row_map, A.graph.entries),
                          numCompressedEntries);
  // Prefix sum to get rowmap
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<rowmap_t, exec_space>(
      A.numRows() + 1, mergedRowmap);
  entries_t mergedEntries("SortedMerged entries", numCompressedEntries);
  values_t mergedValues("SortedMerged values", numCompressedEntries);
  // Compute merged entries and values
  Kokkos::parallel_for(
      range_t(0, A.numRows()),
      Impl::MatrixMergedEntriesFunctor<c_rowmap_t, entries_t, values_t>(
          A.graph.row_map, A.graph.entries, A.values, mergedRowmap,
          mergedEntries, mergedValues));
  // Finally, construct the new compressed matrix
  return crsMat_t("SortedMerged", A.numRows(), A.numCols(),
                  numCompressedEntries, mergedValues, mergedRowmap,
                  mergedEntries);
}

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in,
                          const entries_t& entries_in, rowmap_t& rowmap_out,
                          entries_t& entries_out) {
  using size_type      = typename rowmap_t::non_const_value_type;
  using lno_t          = typename entries_t::non_const_value_type;
  using range_t        = Kokkos::RangePolicy<exec_space>;
  using const_rowmap_t = typename rowmap_t::const_type;
  lno_t numRows        = rowmap_in.extent(0);
  if (numRows <= 1) {
    // Matrix has zero rows
    rowmap_out  = rowmap_t();
    entries_out = entries_t();
    return;
  }
  numRows--;
  // Sort in place
  sort_crs_graph<exec_space, const_rowmap_t, entries_t>(rowmap_in, entries_in);
  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  rowmap_out = rowmap_t(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "SortedMerged rowmap"),
      numRows + 1);
  size_type numCompressedEntries = 0;
  Kokkos::parallel_reduce(range_t(0, numRows),
                          Impl::MergedRowmapFunctor<rowmap_t, entries_t>(
                              rowmap_out, rowmap_in, entries_in),
                          numCompressedEntries);
  // Prefix sum to get rowmap
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<rowmap_t, exec_space>(
      numRows + 1, rowmap_out);
  entries_out = entries_t("SortedMerged entries", numCompressedEntries);
  // Compute merged entries and values
  Kokkos::parallel_for(
      range_t(0, numRows),
      Impl::GraphMergedEntriesFunctor<const_rowmap_t, entries_t>(
          rowmap_in, entries_in, rowmap_out, entries_out));
}

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(const crsGraph_t& G) {
  using rowmap_t  = typename crsGraph_t::row_map_type::non_const_type;
  using entries_t = typename crsGraph_t::entries_type;
  static_assert(
      !std::is_const<typename entries_t::value_type>::value,
      "sort_and_merge_graph requires StaticCrsGraph entries to be non-const.");
  rowmap_t mergedRowmap;
  entries_t mergedEntries;
  sort_and_merge_graph<typename crsGraph_t::execution_space, rowmap_t,
                       entries_t>(G.row_map, G.entries, mergedRowmap,
                                  mergedEntries);
  return crsGraph_t(mergedEntries, mergedRowmap);
}

}  // namespace KokkosSparse

namespace KokkosKernels {

// ----------------------------------
// BSR matrix/graph sorting utilities
// ----------------------------------

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t,
          typename lno_t = typename entries_t::non_const_value_type>
[[deprecated]] void sort_bsr_matrix(const lno_t blockdim,
                                    const rowmap_t& rowmap,
                                    const entries_t& entries,
                                    const values_t& values) {
  KokkosSparse::sort_bsr_matrix(blockdim, rowmap, entries, values);
}

template <typename bsrMat_t>
[[deprecated]] void sort_bsr_matrix(const bsrMat_t& A) {
  KokkosSparse::sort_bsr_matrix(A);
}

// ----------------------------------
// CRS matrix/graph sorting utilities
// ----------------------------------

// The sort_crs* functions sort the adjacent column list for each row into
// ascending order.

template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t>
[[deprecated]] void sort_crs_matrix(const rowmap_t& rowmap,
                                    const entries_t& entries,
                                    const values_t& values) {
  KokkosSparse::sort_crs_matrix<execution_space, rowmap_t, entries_t>(
      rowmap, entries, values);
}

template <typename crsMat_t>
[[deprecated]] void sort_crs_matrix(const crsMat_t& A) {
  KokkosSparse::sort_crs_matrix(A);
}

template <typename execution_space, typename rowmap_t, typename entries_t>
[[deprecated]] void sort_crs_graph(const rowmap_t& rowmap,
                                   const entries_t& entries) {
  KokkosSparse::sort_crs_graph<execution_space, rowmap_t, entries_t>(rowmap,
                                                                     entries);
}

template <typename crsGraph_t>
[[deprecated]] void sort_crs_graph(const crsGraph_t& G) {
  KokkosSparse::sort_crs_graph(G);
}

// sort_and_merge_matrix produces a new matrix which is equivalent to A but is
// sorted and has no duplicated entries: each (i, j) is unique. Values for
// duplicated entries are summed.
template <typename crsMat_t>
[[deprecated]] crsMat_t sort_and_merge_matrix(const crsMat_t& A) {
  KokkosSparse::sort_and_merge_matrix(A);
}

template <typename crsGraph_t>
[[deprecated]] crsGraph_t sort_and_merge_graph(const crsGraph_t& G) {
  KokkosSparse::sort_and_merge_graph(G);
}

template <typename exec_space, typename rowmap_t, typename entries_t>
[[deprecated]] void sort_and_merge_graph(
    const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
    rowmap_t& rowmap_out, entries_t& entries_out) {
  KokkosSparse::sort_and_merge_graph(rowmap_in, entries_in, rowmap_out,
                                     entries_out);
}

// For backward compatibility: keep the public interface accessible in
// KokkosKernels::Impl::
namespace Impl {
template <typename execution_space, typename rowmap_t, typename entries_t>
[[deprecated]] void sort_crs_graph(const rowmap_t& rowmap,
                                   const entries_t& entries) {
  KokkosKernels::sort_crs_graph<execution_space, rowmap_t, entries_t>(rowmap,
                                                                      entries);
}

template <typename execution_space, typename rowmap_t, typename entries_t,
          typename values_t>
[[deprecated]] void sort_crs_matrix(const rowmap_t& rowmap,
                                    const entries_t& entries,
                                    const values_t& values) {
  KokkosKernels::sort_crs_matrix<execution_space, rowmap_t, entries_t,
                                 values_t>(rowmap, entries, values);
}

template <typename crsMat_t>
[[deprecated]] void sort_crs_matrix(const crsMat_t& A) {
  KokkosKernels::sort_crs_matrix(A);
}

template <typename exec_space, typename rowmap_t, typename entries_t>
[[deprecated]] void sort_and_merge_graph(
    const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
    rowmap_t& rowmap_out, entries_t& entries_out) {
  KokkosKernels::sort_and_merge_graph<exec_space, rowmap_t, entries_t>(
      rowmap_in, entries_in, rowmap_out, entries_out);
}

template <typename crsMat_t>
[[deprecated]] crsMat_t sort_and_merge_matrix(const crsMat_t& A) {
  return KokkosKernels::sort_and_merge_matrix(A);
}

}  // namespace Impl
}  // namespace KokkosKernels

#endif  // _KOKKOSSPARSE_SORTCRS_HPP
