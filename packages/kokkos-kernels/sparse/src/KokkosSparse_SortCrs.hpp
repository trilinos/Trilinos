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

#include "KokkosSparse_sort_crs_impl.hpp"
#include "KokkosSparse_Utils.hpp"

namespace KokkosSparse {

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

// Sort a CRS matrix: within each row, sort entries ascending by column.
// At the same time, permute the values.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries,
                     const values_t& values,
                     typename entries_t::non_const_value_type numCols =
                         Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
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
  using Ordinal = typename entries_t::non_const_value_type;
  // This early return condition covers having 0 or 1 entries,
  // which is also implied by having 0 rows or 0 columns.
  // If only 1 entry, the matrix is already sorted.
  if (entries.extent(0) <= size_t(1)) {
    return;
  }
  Ordinal numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if constexpr (!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
    // On CPUs, use a sequential radix sort within each row.
    Kokkos::parallel_for("sort_crs_matrix[CPU,radix]",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, 0, numRows),
                         Impl::MatrixRadixSortFunctor<rowmap_t, entries_t, values_t>(rowmap, entries, values));
  } else {
    // On GPUs:
    //   If the matrix is highly imbalanced, or has long rows AND the dimensions
    //   are not too large to do one large bulk sort, do that. Otherwise, sort
    //   using one Kokkos thread per row.
    Ordinal avgDeg = (entries.extent(0) + numRows - 1) / numRows;
#ifndef KK_DISABLE_BULK_SORT_BY_KEY
    Ordinal maxDeg   = KokkosSparse::Impl::graph_max_degree(exec, rowmap);
    bool useBulkSort = false;
    if (KokkosSparse::Impl::useBulkSortHeuristic<execution_space>(avgDeg, maxDeg)) {
      // Calculate the true number of columns if user didn't pass it in
      if (numCols == Kokkos::ArithTraits<Ordinal>::max()) {
        KokkosKernels::Impl::kk_view_reduce_max(exec, entries.extent(0), entries, numCols);
        numCols++;
      }
      uint64_t maxBulkKey = (uint64_t)numRows * (uint64_t)numCols;
      useBulkSort         = maxBulkKey / numRows == (uint64_t)numCols;
    }
    if (useBulkSort) {
      auto permutation = KokkosSparse::Impl::computeEntryPermutation(exec, rowmap, entries, numCols);
      // Permutations cannot be done in-place
      Kokkos::View<typename values_t::value_type*, execution_space> origValues(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "origValues"), values.extent(0));
      Kokkos::View<typename entries_t::value_type*, execution_space> origEntries(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "origEntries"), entries.extent(0));
      Kokkos::deep_copy(exec, origValues, values);
      Kokkos::deep_copy(exec, origEntries, entries);
      KokkosSparse::Impl::applyPermutation(exec, permutation, origEntries, entries);
      KokkosSparse::Impl::applyPermutation(exec, permutation, origValues, values);
    } else
#else
    (void)numCols;
#endif
    {
      using TeamPol = Kokkos::TeamPolicy<execution_space>;
      // Can't use bulk sort approach as matrix dimensions are too large.
      // Fall back to parallel thread-level sort within each row.
      Ordinal vectorLength = 1;
      while (vectorLength < avgDeg / 2) {
        vectorLength *= 2;
      }
      if (vectorLength > TeamPol ::vector_length_max()) vectorLength = TeamPol ::vector_length_max();
      Impl::MatrixThreadSortFunctor<TeamPol, Ordinal, rowmap_t, entries_t, values_t> funct(numRows, rowmap, entries,
                                                                                           values);
      Ordinal teamSize = TeamPol(exec, 1, 1, vectorLength).team_size_recommended(funct, Kokkos::ParallelForTag());
      Kokkos::parallel_for("sort_crs_matrix[GPU,bitonic]",
                           TeamPol(exec, (numRows + teamSize - 1) / teamSize, teamSize, vectorLength), funct);
    }
  }
}

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries, const values_t& values,
                     typename entries_t::const_value_type numCols =
                         Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  sort_crs_matrix(execution_space(), rowmap, entries, values, numCols);
}

template <typename rowmap_t, typename entries_t, typename values_t>
void sort_crs_matrix(const rowmap_t& rowmap, const entries_t& entries, const values_t& values,
                     typename entries_t::const_value_type numCols =
                         Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  sort_crs_matrix(typename entries_t::execution_space(), rowmap, entries, values, numCols);
}

template <typename crsMat_t>
void sort_crs_matrix(const typename crsMat_t::execution_space& exec, const crsMat_t& A) {
  sort_crs_matrix(exec, A.graph.row_map, A.graph.entries, A.values, A.numCols());
}

template <typename crsMat_t>
void sort_crs_matrix(const crsMat_t& A) {
  sort_crs_matrix(typename crsMat_t::execution_space(), A.graph.row_map, A.graph.entries, A.values, A.numCols());
}

// Sort a BRS matrix: within each row, sort entries ascending by column and
// permute the values accordingly.
template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t,
          typename Ordinal = typename entries_t::non_const_value_type>
void sort_bsr_matrix(const execution_space& exec, Ordinal blockSize, const rowmap_t& rowmap, const entries_t& entries,
                     const values_t& values,
                     typename entries_t::non_const_value_type numCols =
                         Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  static_assert(std::is_same_v<Ordinal, typename entries_t::non_const_value_type>,
                "sort_bsr_matrix: Ordinal type must match nonconst value type of "
                "entries_t (default template parameter)");
  if (entries.extent(0) <= size_t(1)) {
    return;
  }
  Ordinal numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (numCols == Kokkos::ArithTraits<Ordinal>::max()) {
    KokkosKernels::Impl::kk_view_reduce_max(exec, entries.extent(0), entries, numCols);
    numCols++;
  }
  uint64_t maxBulkKey = (uint64_t)numRows * (uint64_t)numCols;
  if (maxBulkKey / numRows != (uint64_t)numCols)
    throw std::invalid_argument(
        "sort_bsr_matrix: implementation requires that numRows * numCols is "
        "representable in uint64_t");
#ifdef KK_DISABLE_BULK_SORT_BY_KEY
  using TeamPol = Kokkos::TeamPolicy<execution_space>;
  using Offset  = typename rowmap_t::non_const_value_type;
  // Temporary workaround: do not use Kokkos::Experimental::sort_by_key, instead
  // sort bulk keys one row at a time
  auto keys = Impl::generateBulkCrsKeys(exec, rowmap, entries, numCols);
  Kokkos::View<Offset*, execution_space> permutation(Kokkos::view_alloc(Kokkos::WithoutInitializing, "permutation"),
                                                     entries.extent(0));
  KokkosKernels::Impl::sequential_fill(exec, permutation);
  Ordinal vectorLength = 1;
  Ordinal avgDeg       = (entries.extent(0) + numRows - 1) / numRows;
  while (vectorLength < avgDeg / 2) {
    vectorLength *= 2;
  }
  if (vectorLength > TeamPol ::vector_length_max()) vectorLength = TeamPol ::vector_length_max();
  Impl::MatrixThreadSortFunctor<TeamPol, Ordinal, rowmap_t, decltype(keys), decltype(permutation)> funct(
      numRows, rowmap, keys, permutation);
  Ordinal teamSize = TeamPol(exec, 1, 1, vectorLength).team_size_recommended(funct, Kokkos::ParallelForTag());
  Kokkos::parallel_for("sort_bulk_keys_by_row[GPU,bitonic]",
                       TeamPol(exec, (numRows + teamSize - 1) / teamSize, teamSize, vectorLength), funct);
#else
  auto permutation = KokkosSparse::Impl::computeEntryPermutation(exec, rowmap, entries, numCols);
#endif
  // Permutations cannot be done in-place
  Kokkos::View<typename values_t::value_type*, execution_space> origValues(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "origValues"), values.extent(0));
  Kokkos::View<typename entries_t::value_type*, execution_space> origEntries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "origEntries"), entries.extent(0));
  Kokkos::deep_copy(exec, origValues, values);
  Kokkos::deep_copy(exec, origEntries, entries);
  KokkosSparse::Impl::applyPermutation(exec, permutation, origEntries, entries);
  KokkosSparse::Impl::applyPermutationBlockValues(exec, permutation, origValues, values, blockSize);
}

template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t, typename Ordinal>
void sort_bsr_matrix(Ordinal blockdim, const rowmap_t& rowmap, const entries_t& entries, const values_t& values,
                     Ordinal numCols = Kokkos::ArithTraits<Ordinal>::max()) {
  sort_bsr_matrix(execution_space(), blockdim, rowmap, entries, values, numCols);
}

// Sort a BSR matrix (like CRS but single values are replaced with contignous
// blocks)
template <typename bsrMat_t>
void sort_bsr_matrix(const typename bsrMat_t::execution_space& exec, const bsrMat_t& A) {
  // NOTE: unlike rowmap, entries and values are non-const, so we can sort them
  // directly
  sort_bsr_matrix<typename bsrMat_t::execution_space, typename bsrMat_t::row_map_type,
                  typename bsrMat_t::index_type::non_const_type, typename bsrMat_t::values_type::non_const_type>(
      exec, A.blockDim(), A.graph.row_map, A.graph.entries, A.values, A.numCols());
}

template <typename bsrMat_t>
void sort_bsr_matrix(const bsrMat_t& A) {
  sort_bsr_matrix(typename bsrMat_t::execution_space(), A);
}

// Sort a CRS graph: within each row, sort entries ascending by column.
template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries,
                    typename entries_t::non_const_value_type numCols =
                        Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  using Ordinal = typename entries_t::non_const_value_type;
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename rowmap_t::memory_space>::accessible,
                "sort_crs_graph: rowmap_t is not accessible from the given execution "
                "space");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename entries_t::memory_space>::accessible,
                "sort_crs_graph: entries_t is not accessible from the given execution "
                "space");
  static_assert(!std::is_const_v<typename entries_t::value_type>, "sort_crs_graph: entries_t must not be const-valued");
  Ordinal numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  if (entries.extent(0) <= size_t(1)) {
    return;
  }
  if constexpr (!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
    // If on CPU, sort each row independently. Don't need to know numCols for
    // this.
    Kokkos::parallel_for("sort_crs_graph[CPU,radix]",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, 0, numRows),
                         Impl::GraphRadixSortFunctor<rowmap_t, entries_t>(rowmap, entries));
  } else {
    // On GPUs:
    //   If the graph is highly imbalanced AND the dimensions are not too large
    //   to do one large bulk sort, do that. Otherwise, sort using one Kokkos
    //   thread per row.
    Ordinal avgDeg = (entries.extent(0) + numRows - 1) / numRows;
#ifndef KK_DISABLE_BULK_SORT_BY_KEY
    Ordinal maxDeg   = KokkosSparse::Impl::graph_max_degree(exec, rowmap);
    bool useBulkSort = false;
    if (KokkosSparse::Impl::useBulkSortHeuristic<execution_space>(avgDeg, maxDeg)) {
      // Calculate the true number of columns if user didn't pass it in
      if (numCols == Kokkos::ArithTraits<Ordinal>::max()) {
        KokkosKernels::Impl::kk_view_reduce_max(exec, entries.extent(0), entries, numCols);
        numCols++;
      }
      uint64_t maxBulkKey = (uint64_t)numRows * (uint64_t)numCols;
      useBulkSort         = maxBulkKey / numRows == (uint64_t)numCols;
    }
    if (useBulkSort) {
      auto keys = KokkosSparse::Impl::generateBulkCrsKeys(exec, rowmap, entries, numCols);
      Kokkos::Experimental::sort_by_key(exec, keys, entries);
    } else
#else
    (void)numCols;
#endif
    {
      using TeamPol = Kokkos::TeamPolicy<execution_space>;
      // Fall back to thread-level sort within each row
      Ordinal vectorLength = 1;
      while (vectorLength < avgDeg / 2) {
        vectorLength *= 2;
      }
      if (vectorLength > TeamPol ::vector_length_max()) vectorLength = TeamPol ::vector_length_max();

      Impl::GraphThreadSortFunctor<TeamPol, Ordinal, rowmap_t, entries_t> funct(numRows, rowmap, entries);
      Ordinal teamSize = TeamPol(exec, 1, 1, vectorLength).team_size_recommended(funct, Kokkos::ParallelForTag());
      Kokkos::parallel_for("sort_crs_graph[GPU,bitonic]",
                           TeamPol(exec, (numRows + teamSize - 1) / teamSize, teamSize, vectorLength), funct);
    }
  }
}

template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries) {
  sort_crs_graph(execution_space(), rowmap, entries);
}

template <typename rowmap_t, typename entries_t>
typename std::enable_if_t<Kokkos::is_view_v<rowmap_t>> sort_crs_graph(
    const rowmap_t& rowmap, const entries_t& entries,
    typename entries_t::const_value_type& numCols =
        Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  sort_crs_graph(typename entries_t::execution_space(), rowmap, entries, numCols);
}

template <typename execution_space, typename crsGraph_t>
typename std::enable_if_t<Kokkos::is_execution_space_v<execution_space>> sort_crs_graph(
    const execution_space& exec, const crsGraph_t& G,
    typename crsGraph_t::entries_type::const_value_type& numCols =
        Kokkos::ArithTraits<typename crsGraph_t::entries_type::non_const_value_type>::max()) {
  sort_crs_graph(exec, G.row_map, G.entries, numCols);
}

template <typename crsGraph_t>
void sort_crs_graph(const crsGraph_t& G,
                    typename crsGraph_t::entries_type::const_value_type& numCols =
                        Kokkos::ArithTraits<typename crsGraph_t::entries_type::non_const_value_type>::max()) {
  sort_crs_graph(typename crsGraph_t::execution_space(), G, numCols);
}

template <typename exec_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                           const entries_t& entries_in, const values_t& values_in, rowmap_t& rowmap_out,
                           entries_t& entries_out, values_t& values_out,
                           typename entries_t::const_value_type& numCols =
                               Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  using nc_rowmap_t = typename rowmap_t::non_const_type;
  using Offset      = typename nc_rowmap_t::value_type;
  using Ordinal     = typename entries_t::value_type;
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

  Ordinal numRows = rowmap_in.extent(0) ? Ordinal(rowmap_in.extent(0) - 1) : Ordinal(0);
  Offset nnz      = entries_in.extent(0);

  if (numRows == 0) {
    rowmap_out  = typename rowmap_t::non_const_type("SortedMerged rowmap", rowmap_in.extent(0));
    entries_out = entries_t();
    values_out  = values_t();
    return;
  }

  sort_crs_matrix(exec, rowmap_in, entries_in, values_in, numCols);

  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  nc_rowmap_t nc_rowmap_out(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged rowmap"), numRows + 1);
  Offset numCompressedEntries = 0;
  Kokkos::parallel_reduce("KokkosSparse::Impl::MergedRowmapFunctor", range_t(exec, 0, numRows),
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
  Kokkos::parallel_for("KokkosSparse::Impl::MatrixMergedEntriesFunctor", range_t(exec, 0, numRows),
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

  sort_and_merge_matrix(exec, A.graph.row_map, A.graph.entries, A.values, rowmap_out, entries_out, values_out,
                        A.numCols());

  return crsMat_t("SortedMerged", A.numRows(), A.numCols(), values_out.extent(0), values_out, rowmap_out, entries_out);
}

template <typename crsMat_t>
crsMat_t sort_and_merge_matrix(const crsMat_t& A) {
  return sort_and_merge_matrix(typename crsMat_t::execution_space(), A);
}

template <typename exec_space, typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                           const values_t& values_in, rowmap_t& rowmap_out, entries_t& entries_out,
                           values_t& values_out,
                           typename entries_t::const_value_type& numCols =
                               Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  sort_and_merge_matrix(exec_space(), rowmap_in, entries_in, values_in, rowmap_out, entries_out, values_out, numCols);
}

template <typename rowmap_t, typename entries_t, typename values_t>
void sort_and_merge_matrix(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                           const values_t& values_in, rowmap_t& rowmap_out, entries_t& entries_out,
                           values_t& values_out,
                           typename entries_t::const_value_type& numCols =
                               Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  sort_and_merge_matrix(typename entries_t::execution_space(), rowmap_in, entries_in, values_in, rowmap_out,
                        entries_out, values_out, numCols);
}

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                          const entries_t& entries_in, rowmap_t& rowmap_out, entries_t& entries_out,
                          typename entries_t::const_value_type& numCols =
                              Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  using Offset      = typename rowmap_t::non_const_value_type;
  using Ordinal     = typename entries_t::value_type;
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

  Ordinal numRows = rowmap_in.extent(0) ? rowmap_in.extent(0) - 1 : 0;
  if (numRows == 0) {
    rowmap_out  = typename rowmap_t::non_const_type("SortedMerged rowmap", rowmap_in.extent(0));
    entries_out = entries_t();
    return;
  }
  // Sort in place
  sort_crs_graph(exec, rowmap_in, entries_in, numCols);
  // Count entries per row into a new rowmap, in terms of merges that can be
  // done
  nc_rowmap_t nc_rowmap_out(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "SortedMerged rowmap"), numRows + 1);
  Offset numCompressedEntries = 0;
  Kokkos::parallel_reduce("KokkosSparse::Impl::MergedRowmapFunctor", range_t(exec, 0, numRows),
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
  Kokkos::parallel_for(
      "KokkosSparse::Impl::GraphMergedEntriesFunctor", range_t(exec, 0, numRows),
      Impl::GraphMergedEntriesFunctor<rowmap_t, entries_t>(rowmap_orig, entries_orig, rowmap_out, entries_out));
}

template <typename exec_space, typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                          rowmap_t& rowmap_out, entries_t& entries_out,
                          typename entries_t::const_value_type& numCols =
                              Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  return sort_and_merge_graph(exec_space(), rowmap_in, entries_in, rowmap_out, entries_out, numCols);
}

template <typename rowmap_t, typename entries_t>
void sort_and_merge_graph(const typename rowmap_t::const_type& rowmap_in, const entries_t& entries_in,
                          rowmap_t& rowmap_out, entries_t& entries_out,
                          typename entries_t::const_value_type& numCols =
                              Kokkos::ArithTraits<typename entries_t::non_const_value_type>::max()) {
  return sort_and_merge_graph(typename entries_t::execution_space(), rowmap_in, entries_in, rowmap_out, entries_out,
                              numCols);
}

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(
    const typename crsGraph_t::execution_space& exec, const crsGraph_t& G,
    typename crsGraph_t::entries_type::const_value_type& numCols =
        Kokkos::ArithTraits<typename crsGraph_t::entries_type::non_const_value_type>::max()) {
  using rowmap_t  = typename crsGraph_t::row_map_type::non_const_type;
  using entries_t = typename crsGraph_t::entries_type;
  static_assert(!std::is_const<typename entries_t::value_type>::value,
                "sort_and_merge_graph requires StaticCrsGraph entries to be non-const.");
  rowmap_t mergedRowmap;
  entries_t mergedEntries;
  sort_and_merge_graph(exec, G.row_map, G.entries, mergedRowmap, mergedEntries, numCols);
  return crsGraph_t(mergedEntries, mergedRowmap);
}

template <typename crsGraph_t>
crsGraph_t sort_and_merge_graph(
    const crsGraph_t& G, typename crsGraph_t::entries_type::const_value_type& numCols =
                             Kokkos::ArithTraits<typename crsGraph_t::entries_type::non_const_value_type>::max()) {
  return sort_and_merge_graph(typename crsGraph_t::execution_space(), G, numCols);
}

}  // namespace KokkosSparse

#endif  // _KOKKOSSPARSE_SORTCRS_HPP
