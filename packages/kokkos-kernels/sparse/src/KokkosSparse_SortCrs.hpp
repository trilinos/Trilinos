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
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t,
          typename lno_t = typename entries_t::non_const_value_type>
void sort_bsr_matrix(const lno_t blockdim, const rowmap_t& rowmap, const entries_t& entries, const values_t& values);

// Sort a BRS matrix on the given execution space instance: within each row,
// sort entries ascending by column and permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t,
          typename lno_t = typename entries_t::non_const_value_type>
void sort_bsr_matrix(const execution_space& exec, const lno_t blockdim, const rowmap_t& rowmap,
                     const entries_t& entries, const values_t& values);

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename bsrMat_t>
void sort_bsr_matrix(const bsrMat_t& A);

// Sort a BRS matrix on the given execution space instance: within each row,
// sort entries ascending by column and permute the values accordingly.
template <typename bsrMat_t>
void sort_bsr_matrix(const typename bsrMat_t::execution_space& exec, const bsrMat_t& A);

// ----------------------------------
// CRS matrix/graph sorting utilities
// ----------------------------------

// The sort_crs* functions sort the adjacent column list for each row into
// ascending order. Each version either takes an execution space instance as a
// parameter, or uses the default instance.

// sort_and_merge_matrix produces a new matrix which is equivalent to A but is
// sorted and has no duplicated entries: each (i, j) is unique. Values for
// duplicated entries are summed. Each version either takes an execution space
// instance as a parameter, or uses the default instance. If there are no
// duplicated entries in A, A is sorted and returned (instead of a newly
// allocated matrix).

namespace Impl {

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
struct SortCrsMatrixFunctor {
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t     = typename entries_t::non_const_value_type;
  using scalar_t  = typename values_t::non_const_value_type;
  using team_mem  = typename Kokkos::TeamPolicy<execution_space>::member_type;
  // The functor owns memory for entriesAux, so it can't have
  // MemoryTraits<Unmanaged>
  using entries_managed_t = Kokkos::View<typename entries_t::data_type, typename entries_t::device_type>;
  using values_managed_t  = Kokkos::View<typename values_t::data_type, typename values_t::device_type>;

  SortCrsMatrixFunctor(bool usingRangePol, const rowmap_t& rowmap_, const entries_t& entries_, const values_t& values_)
      : rowmap(rowmap_), entries(entries_), values(values_) {
    if (usingRangePol) {
      entriesAux = entries_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"), entries.extent(0));
      valuesAux  = values_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values aux"), values.extent(0));
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
        (unsigned_lno_t*)entries.data() + rowStart, (unsigned_lno_t*)entriesAux.data() + rowStart,
        values.data() + rowStart, valuesAux.data() + rowStart, rowNum);
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_mem t) const {
    size_type i        = t.league_rank();
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    KokkosKernels::TeamBitonicSort2<lno_t, lno_t, scalar_t, team_mem>(entries.data() + rowStart,
                                                                      values.data() + rowStart, rowNum, t);
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
  using entries_managed_t = Kokkos::View<typename entries_t::data_type, typename entries_t::device_type>;

  SortCrsGraphFunctor(bool usingRangePol, const rowmap_t& rowmap_, const entries_t& entries_)
      : rowmap(rowmap_), entries(entries_) {
    if (usingRangePol) {
      entriesAux = entries_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"), entries.extent(0));
    }
    // otherwise, aux arrays won't be allocated (sorting in place)
  }

  KOKKOS_INLINE_FUNCTION void operator()(const lno_t i) const {
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    // Radix sort requires unsigned keys for comparison
    using unsigned_lno_t = typename std::make_unsigned<lno_t>::type;
    KokkosKernels::SerialRadixSort<lno_t, unsigned_lno_t>((unsigned_lno_t*)entries.data() + rowStart,
                                                          (unsigned_lno_t*)entriesAux.data() + rowStart, rowNum);
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_mem t) const {
    size_type i        = t.league_rank();
    size_type rowStart = rowmap(i);
    size_type rowEnd   = rowmap(i + 1);
    lno_t rowNum       = rowEnd - rowStart;
    KokkosKernels::TeamBitonicSort<lno_t, lno_t, team_mem>(entries.data() + rowStart, rowNum, t);
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
  MergedRowmapFunctor(const rowmap_t& mergedCounts_, const c_rowmap_t& rowmap_, const entries_t& entries_)
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
  MatrixMergedEntriesFunctor(const typename rowmap_t::const_type& rowmap_, const entries_t& entries_,
                             const values_t& values_, const rowmap_t& mergedRowmap_, const entries_t& mergedEntries_,
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

  typename rowmap_t::const_type rowmap;
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
  GraphMergedEntriesFunctor(const typename rowmap_t::const_type& rowmap_, const entries_t& entries_,
                            const rowmap_t& mergedRowmap_, const entries_t& mergedEntries_)
      : rowmap(rowmap_), entries(entries_), mergedRowmap(mergedRowmap_), mergedEntries(mergedEntries_) {}

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

  typename rowmap_t::const_type rowmap;
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

  sort_bsr_functor(row_map_type rowmap_, entries_type entries_, values_type values_, const lno_t blocksize_)
      : rowmap(rowmap_), entries(entries_), values(values_), blocksize(blocksize_) {}

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
        for (lno_t k = 0; k < blocksize; ++k)  // std::swap_ranges(vb, vb + blocksize, vbp);
          Impl::kk_swap(vb[k], vbp[k]);
        done = false;
      }
    }
  }
};

}  // namespace Impl

// Sort a CRS matrix: within each row, sort entries ascending by column.
// At the same time, permute the values.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries,
                     const values_t& values) {
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename rowmap_t::memory_space>::accessible,
                "sort_crs_matrix: rowmap_t is not accessible from the given execution "
                "space");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename entries_t::memory_space>::accessible,
                "sort_crs_matrix: entries_t is not accessible from the given execution "
                "space");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename entries_t::memory_space>::accessible,
                "sort_crs_matrix: values_t is not accessible from the given execution "
                "space");
  static_assert(!std::is_const_v<typename entries_t::value_type>,
                "sort_crs_matrix: entries_t must not be const-valued");
  static_assert(!std::is_const_v<typename values_t::value_type>, "sort_crs_matrix: value_t must not be const-valued");
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  bool useRadix  = !KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows  = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  Impl::SortCrsMatrixFunctor<execution_space, rowmap_t, entries_t, values_t> funct(useRadix, rowmap, entries, values);
  if (useRadix) {
    Kokkos::parallel_for("sort_crs_matrix", Kokkos::RangePolicy<execution_space>(exec, 0, numRows), funct);
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
    team_pol temp(exec, numRows, 1);
    lno_t maxTeamSize = temp.team_size_max(funct, Kokkos::ParallelForTag());
    lno_t teamSize    = std::min(idealTeamSize, maxTeamSize);
    Kokkos::parallel_for("sort_crs_matrix", team_pol(exec, numRows, teamSize), funct);
  }
}

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries, const values_t& values) {
  sort_crs_matrix(execution_space(), rowmap, entries, values);
}

template <typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries, const values_t& values) {
  sort_crs_matrix(typename entries_t::execution_space(), rowmap, entries, values);
}

template <typename crsMat_t>
void sort_crs_matrix(const typename crsMat_t::execution_space& exec, const crsMat_t& A) {
  sort_crs_matrix(exec, A.graph.row_map, A.graph.entries, A.values);
}

template <typename crsMat_t>
void sort_crs_matrix(const crsMat_t& A) {
  sort_crs_matrix(typename crsMat_t::execution_space(), A.graph.row_map, A.graph.entries, A.values);
}

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t, typename lno_t>
void sort_bsr_matrix(const execution_space& exec, const lno_t blockdim, const rowmap_t& rowmap,
                     const entries_t& entries, const values_t& values) {
  // TODO: this is O(N^2) mock for debugging - do regular implementation based
  // on Radix/Bitonic sort (like CSR) IDEA: maybe we need only one general
  // Radix2/Bitonic2 and CSR sorting may call it with blockSize=1 ?
  lno_t numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  const lno_t blocksize = blockdim * blockdim;

  assert(values.extent(0) == entries.extent(0) * blocksize);
  Impl::sort_bsr_functor<rowmap_t, entries_t, values_t> bsr_sorter(rowmap, entries, values, blocksize);
  Kokkos::parallel_for("sort_bsr_matrix", Kokkos::RangePolicy<execution_space>(exec, 0, numRows), bsr_sorter);
}

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t, typename lno_t>
void sort_bsr_matrix(const lno_t blockdim, const rowmap_t& rowmap, const entries_t& entries, const values_t& values) {
  sort_bsr_matrix(execution_space(), blockdim, rowmap, entries, values);
}

// Sort a BSR matrix (like CRS but single values are replaced with contignous
// blocks)
template <typename bsrMat_t>
void sort_bsr_matrix(const typename bsrMat_t::execution_space& exec, const bsrMat_t& A) {
  // NOTE: unlike rowmap, entries and values are non-const, so we can sort them
  // directly
  sort_bsr_matrix<typename bsrMat_t::execution_space, typename bsrMat_t::row_map_type,
                  typename bsrMat_t::index_type::non_const_type, typename bsrMat_t::values_type::non_const_type>(
      exec, A.blockDim(), A.graph.row_map, A.graph.entries, A.values);
}

template <typename bsrMat_t>
void sort_bsr_matrix(const bsrMat_t& A) {
  sort_bsr_matrix(typename bsrMat_t::execution_space(), A);
}

// Sort a CRS graph: within each row, sort entries ascending by column.
template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries) {
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename rowmap_t::memory_space>::accessible,
                "sort_crs_graph: rowmap_t is not accessible from the given execution "
                "space");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename entries_t::memory_space>::accessible,
                "sort_crs_graph: entries_t is not accessible from the given execution "
                "space");
  static_assert(!std::is_const_v<typename entries_t::value_type>, "sort_crs_graph: entries_t must not be const-valued");
  bool useRadix = !KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numRows == 0) return;
  Impl::SortCrsGraphFunctor<execution_space, rowmap_t, entries_t> funct(useRadix, rowmap, entries);
  if (useRadix) {
    Kokkos::parallel_for("sort_crs_graph", Kokkos::RangePolicy<execution_space>(exec, 0, numRows), funct);
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
    team_pol temp(exec, numRows, 1);
    lno_t maxTeamSize = temp.team_size_max(funct, Kokkos::ParallelForTag());
    lno_t teamSize    = std::min(idealTeamSize, maxTeamSize);
    Kokkos::parallel_for("sort_crs_graph", team_pol(exec, numRows, teamSize), funct);
  }
}

template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries) {
  sort_crs_graph(execution_space(), rowmap, entries);
}

// This overload covers 2 cases, while allowing all template args to be deduced:
//  - sort_crs_graph(exec, G)
//  - sort_crs_graph(rowmap, entries)
template <typename Arg1, typename Arg2>
void sort_crs_graph(const Arg1& a1, const Arg2& a2) {
  if constexpr (Kokkos::is_execution_space_v<Arg1>) {
    // a1 is an exec instance, a2 is a graph
    sort_crs_graph(a1, a2.row_map, a2.entries);
  } else if constexpr (Kokkos::is_view_v<Arg1>) {
    // a1 is rowmap, a2 is entries
    sort_crs_graph(typename Arg2::execution_space(), a1, a2);
  } else {
    static_assert(Arg1::doesnthavethisthing,
                  "sort_crs_graph(arg1, arg2): expect either (exec, G) or "
                  "(rowmap, entries)");
  }
}

template <typename crsGraph_t>
void sort_crs_graph(const crsGraph_t& G) {
  sort_crs_graph(typename crsGraph_t::execution_space(), G);
}

template <typename exec_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                           const entries_t& entries_in, const values_t& values_in, rowmap_t& rowmap_out,
                           entries_t& entries_out, values_t& values_out) {
  using nc_rowmap_t = typename rowmap_t::non_const_type;
  using size_type   = typename nc_rowmap_t::value_type;
  using ordinal_t   = typename entries_t::value_type;
  using range_t     = Kokkos::RangePolicy<exec_space>;
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename rowmap_t::memory_space>::accessible,
                "sort_and_merge_matrix: rowmap_t is not accessible from the given "
                "execution space");
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename entries_t::memory_space>::accessible,
                "sort_and_merge_matrix: entries_t is not accessible from the given "
                "execution space");
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename entries_t::memory_space>::accessible,
                "sort_and_merge_matrix: values_t is not accessible from the given "
                "execution space");
  static_assert(!std::is_const_v<typename entries_t::value_type>,
                "sort_and_merge_matrix: entries_t must not be const-valued");
  static_assert(!std::is_const_v<typename values_t::value_type>,
                "sort_and_merge_matrix: value_t must not be const-valued");

  ordinal_t numRows = rowmap_in.extent(0) ? ordinal_t(rowmap_in.extent(0) - 1) : ordinal_t(0);
  size_type nnz     = entries_in.extent(0);

  if (numRows == 0) {
    rowmap_out  = typename rowmap_t::non_const_type("SortedMerged rowmap", rowmap_in.extent(0));
    entries_out = entries_t();
    values_out  = values_t();
    return;
  }

  sort_crs_matrix(exec, rowmap_in, entries_in, values_in);

  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  nc_rowmap_t nc_rowmap_out(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged rowmap"), numRows + 1);
  size_type numCompressedEntries = 0;
  Kokkos::parallel_reduce(range_t(exec, 0, numRows),
                          Impl::MergedRowmapFunctor<nc_rowmap_t, entries_t>(nc_rowmap_out, rowmap_in, entries_in),
                          numCompressedEntries);
  if (nnz == numCompressedEntries) {
    // No merges to do, so just return A. Save the time of allocating and
    // filling a copy.
    if constexpr (std::is_const_v<typename rowmap_t::value_type>) {
      rowmap_out = rowmap_in;
    } else {
      // rowmap_t is non-const, so we can't directly assign rowmap_in to
      // rowmap_out. Forced to deep copy it to maintain const-correctness.
      Kokkos::deep_copy(exec, nc_rowmap_out, rowmap_in);
      rowmap_out = nc_rowmap_out;
    }
    entries_out = entries_in;
    values_out  = values_in;
    return;
  }
  // Have to do the compression. Create a _shallow_ copy of the input
  // to preserve it, in case the input and output views are identical
  // references.
  auto rowmap_orig  = rowmap_in;
  auto entries_orig = entries_in;
  auto values_orig  = values_in;
  // Prefix sum to get rowmap
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<exec_space>(exec, numRows + 1, nc_rowmap_out);
  rowmap_out = nc_rowmap_out;
  entries_out =
      entries_t(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged entries"), numCompressedEntries);
  values_out =
      values_t(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged values"), numCompressedEntries);
  // Compute merged entries and values
  Kokkos::parallel_for(range_t(exec, 0, numRows),
                       Impl::MatrixMergedEntriesFunctor<rowmap_t, entries_t, values_t>(
                           rowmap_orig, entries_orig, values_orig, rowmap_out, entries_out, values_out));
}

// Sort the rows of matrix, and merge duplicate entries.
template <typename crsMat_t>
crsMat_t sort_and_merge_matrix(const typename crsMat_t::execution_space& exec, const crsMat_t& A) {
  using rowmap_t  = typename crsMat_t::row_map_type;
  using entries_t = typename crsMat_t::index_type;
  using values_t  = typename crsMat_t::values_type;

  rowmap_t rowmap_out;
  entries_t entries_out;
  values_t values_out;

  sort_and_merge_matrix(exec, A.graph.row_map, A.graph.entries, A.values, rowmap_out, entries_out, values_out);

  return crsMat_t("SortedMerged", A.numRows(), A.numCols(), values_out.extent(0), values_out, rowmap_out, entries_out);
}

template <typename crsMat_t>
crsMat_t sort_and_merge_matrix(const crsMat_t& A) {
  return sort_and_merge_matrix(typename crsMat_t::execution_space(), A);
}

template <typename exec_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                           const values_t& values_in, rowmap_t& rowmap_out, entries_t& entries_out,
                           values_t& values_out) {
  sort_and_merge_matrix(exec_space(), rowmap_in, entries_in, values_in, rowmap_out, entries_out, values_out);
}

template <typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                           const values_t& values_in, rowmap_t& rowmap_out, entries_t& entries_out,
                           values_t& values_out) {
  sort_and_merge_matrix(typename entries_t::execution_space(), rowmap_in, entries_in, values_in, rowmap_out,
                        entries_out, values_out);
}

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                          const entries_t& entries_in, rowmap_t& rowmap_out, entries_t& entries_out) {
  using size_type   = typename rowmap_t::non_const_value_type;
  using lno_t       = typename entries_t::value_type;
  using range_t     = Kokkos::RangePolicy<exec_space>;
  using nc_rowmap_t = typename rowmap_t::non_const_type;
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename rowmap_t::memory_space>::accessible,
                "sort_and_merge_graph: rowmap_t is not accessible from the given "
                "execution space");
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename entries_t::memory_space>::accessible,
                "sort_and_merge_graph: entries_t is not accessible from the given "
                "execution space");
  static_assert(!std::is_const_v<typename entries_t::value_type>,
                "sort_and_merge_graph: entries_t must not be const-valued");

  lno_t numRows = rowmap_in.extent(0) ? rowmap_in.extent(0) - 1 : 0;
  if (numRows == 0) {
    rowmap_out  = typename rowmap_t::non_const_type("SortedMerged rowmap", rowmap_in.extent(0));
    entries_out = entries_t();
    return;
  }
  // Sort in place
  sort_crs_graph(exec, rowmap_in, entries_in);
  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  nc_rowmap_t nc_rowmap_out(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged rowmap"), numRows + 1);
  size_type numCompressedEntries = 0;
  Kokkos::parallel_reduce(range_t(exec, 0, numRows),
                          Impl::MergedRowmapFunctor<rowmap_t, entries_t>(nc_rowmap_out, rowmap_in, entries_in),
                          numCompressedEntries);
  if (entries_in.extent(0) == size_t(numCompressedEntries)) {
    // No merges to perform, so the output rowmap is unchanged and we can just
    // return the now-sorted entries_in.
    if constexpr (std::is_const_v<typename rowmap_t::value_type>) {
      rowmap_out = rowmap_in;
    } else {
      // rowmap_t is non-const, so we can't directly assign rowmap_in to
      // rowmap_out. Forced to deep copy it to maintain const-correctness.
      Kokkos::deep_copy(exec, nc_rowmap_out, rowmap_in);
      rowmap_out = nc_rowmap_out;
    }
    entries_out = entries_in;
    return;
  }
  // Have to do the compression. Create a _shallow_ copy of the input
  // to preserve it, in case the input and output views are identical
  // references.
  auto rowmap_orig  = rowmap_in;
  auto entries_orig = entries_in;
  // Prefix sum to get rowmap.
  // In the case where the output rowmap is the same as the input, we could just
  // assign "rowmap_out = rowmap_in" except that would break const-correctness.
  // Can skip filling the entries, however.
  KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<exec_space>(exec, numRows + 1, nc_rowmap_out);
  rowmap_out = nc_rowmap_out;
  entries_out =
      entries_t(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged entries"), numCompressedEntries);
  // Compute merged entries and values
  Kokkos::parallel_for(range_t(exec, 0, numRows), Impl::GraphMergedEntriesFunctor<rowmap_t, entries_t>(
                                                      rowmap_orig, entries_orig, rowmap_out, entries_out));
}

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                          rowmap_t& rowmap_out, entries_t& entries_out) {
  return sort_and_merge_graph(exec_space(), rowmap_in, entries_in, rowmap_out, entries_out);
}

template <typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                          rowmap_t& rowmap_out, entries_t& entries_out) {
  return sort_and_merge_graph(typename entries_t::execution_space(), rowmap_in, entries_in, rowmap_out, entries_out);
}

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(const typename crsGraph_t::execution_space& exec, const crsGraph_t& G) {
  using rowmap_t  = typename crsGraph_t::row_map_type::non_const_type;
  using entries_t = typename crsGraph_t::entries_type;
  static_assert(!std::is_const<typename entries_t::value_type>::value,
                "sort_and_merge_graph requires StaticCrsGraph entries to be non-const.");
  rowmap_t mergedRowmap;
  entries_t mergedEntries;
  sort_and_merge_graph(exec, G.row_map, G.entries, mergedRowmap, mergedEntries);
  return crsGraph_t(mergedEntries, mergedRowmap);
}

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(const crsGraph_t& G) {
  return sort_and_merge_graph(typename crsGraph_t::execution_space(), G);
}

}  // namespace KokkosSparse

namespace KokkosKernels {

// ----------------------------------
// BSR matrix/graph sorting utilities
// ----------------------------------

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t,
          typename lno_t = typename entries_t::non_const_value_type>
[[deprecated]] void sort_bsr_matrix(const lno_t blockdim, const rowmap_t& rowmap, const entries_t& entries,
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

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
[[deprecated]] void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries, const values_t& values) {
  KokkosSparse::sort_crs_matrix<execution_space, rowmap_t, entries_t>(rowmap, entries, values);
}

template <typename crsMat_t>
[[deprecated]] void sort_crs_matrix(const crsMat_t& A) {
  KokkosSparse::sort_crs_matrix(A);
}

template <typename execution_space, typename rowmap_t, typename entries_t>
[[deprecated]] void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries) {
  KokkosSparse::sort_crs_graph<execution_space, rowmap_t, entries_t>(rowmap, entries);
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
[[deprecated]] void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                                         rowmap_t& rowmap_out, entries_t& entries_out) {
  KokkosSparse::sort_and_merge_graph(rowmap_in, entries_in, rowmap_out, entries_out);
}

}  // namespace KokkosKernels

#endif  // _KOKKOSSPARSE_SORTCRS_HPP
