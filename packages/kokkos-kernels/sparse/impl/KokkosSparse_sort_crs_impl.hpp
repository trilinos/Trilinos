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
#ifndef _KOKKOSSPARSE_SORTCRS_IMPL_HPP
#define _KOKKOSSPARSE_SORTCRS_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "KokkosKernels_Sorting.hpp"

// Workaround for issue with Kokkos::Experimental::sort_by_key, with nvcc and OpenMP enabled
// (Kokkos issue #7036, fixed in 4.4 release)
// Once support for Kokkos < 4.4 is dropped,
// all code inside "ifdef KK_DISABLE_BULK_SORT_BY_KEY" can be deleted.
#if (KOKKOS_VERSION < 40400) && defined(KOKKOS_ENABLE_CUDA)
#define KK_DISABLE_BULK_SORT_BY_KEY
#endif

namespace KokkosSparse {
namespace Impl {

template <typename rowmap_t, typename entries_t, typename values_t>
struct MatrixRadixSortFunctor {
  using Offset          = typename rowmap_t::non_const_value_type;
  using Ordinal         = typename entries_t::non_const_value_type;
  using UnsignedOrdinal = typename std::make_unsigned<Ordinal>::type;
  using Scalar          = typename values_t::non_const_value_type;
  // The functor owns memory for entriesAux, so it can't have
  // MemoryTraits<Unmanaged>
  using entries_managed_t = Kokkos::View<typename entries_t::data_type, typename entries_t::device_type>;
  using values_managed_t  = Kokkos::View<typename values_t::data_type, typename values_t::device_type>;

  MatrixRadixSortFunctor(const rowmap_t& rowmap_, const entries_t& entries_, const values_t& values_)
      : rowmap(rowmap_), entries(entries_), values(values_) {
    entriesAux = entries_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"), entries.extent(0));
    valuesAux  = values_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values aux"), values.extent(0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Ordinal i) const {
    Offset rowStart = rowmap(i);
    Offset rowEnd   = rowmap(i + 1);
    Ordinal rowNum  = rowEnd - rowStart;
    // Radix sort requires unsigned keys for comparison
    KokkosKernels::SerialRadixSort2<Ordinal, UnsignedOrdinal, Scalar>(
        (UnsignedOrdinal*)entries.data() + rowStart, (UnsignedOrdinal*)entriesAux.data() + rowStart,
        values.data() + rowStart, valuesAux.data() + rowStart, rowNum);
  }

  rowmap_t rowmap;
  entries_t entries;
  entries_managed_t entriesAux;
  values_t values;
  values_managed_t valuesAux;
};

template <typename Policy, typename Ordinal, typename rowmap_t, typename entries_t, typename values_t>
struct MatrixThreadSortFunctor {
  using Offset = typename rowmap_t::non_const_value_type;

  MatrixThreadSortFunctor(Ordinal numRows_, const rowmap_t& rowmap_, const entries_t& entries_, const values_t& values_)
      : numRows(numRows_), rowmap(rowmap_), entries(entries_), values(values_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const typename Policy::member_type& t) const {
    Ordinal i = t.league_rank() * t.team_size() + t.team_rank();
    if (i >= numRows) return;
    Offset rowStart = rowmap(i);
    Offset rowEnd   = rowmap(i + 1);
    auto rowEntries = Kokkos::subview(entries, Kokkos::make_pair(rowStart, rowEnd));
    auto rowValues  = Kokkos::subview(values, Kokkos::make_pair(rowStart, rowEnd));
    Kokkos::Experimental::sort_by_key_thread(t, rowEntries, rowValues);
  }

  Ordinal numRows;
  rowmap_t rowmap;
  entries_t entries;
  values_t values;
};

template <typename rowmap_t, typename entries_t>
struct GraphRadixSortFunctor {
  using Offset          = typename rowmap_t::non_const_value_type;
  using Ordinal         = typename entries_t::non_const_value_type;
  using UnsignedOrdinal = typename std::make_unsigned<Ordinal>::type;
  // The functor owns memory for entriesAux, so it can't have
  // MemoryTraits<Unmanaged>
  using entries_managed_t = Kokkos::View<typename entries_t::data_type, typename entries_t::device_type>;

  GraphRadixSortFunctor(const rowmap_t& rowmap_, const entries_t& entries_) : rowmap(rowmap_), entries(entries_) {
    entriesAux = entries_managed_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries aux"), entries.extent(0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Ordinal i) const {
    Offset rowStart = rowmap(i);
    Offset rowEnd   = rowmap(i + 1);
    Ordinal rowNum  = rowEnd - rowStart;
    // Radix sort requires unsigned keys for comparison
    KokkosKernels::SerialRadixSort<Ordinal, UnsignedOrdinal>((UnsignedOrdinal*)entries.data() + rowStart,
                                                             (UnsignedOrdinal*)entriesAux.data() + rowStart, rowNum);
  }

  rowmap_t rowmap;
  entries_t entries;
  entries_managed_t entriesAux;
};

template <typename Policy, typename Ordinal, typename rowmap_t, typename entries_t>
struct GraphThreadSortFunctor {
  using Offset = typename rowmap_t::non_const_value_type;

  GraphThreadSortFunctor(Ordinal numRows_, const rowmap_t& rowmap_, const entries_t& entries_)
      : numRows(numRows_), rowmap(rowmap_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const typename Policy::member_type& t) const {
    Ordinal i = t.league_rank() * t.team_size() + t.team_rank();
    if (i >= numRows) return;
    Offset rowStart = rowmap(i);
    Offset rowEnd   = rowmap(i + 1);
    auto rowEntries = Kokkos::subview(entries, Kokkos::make_pair(rowStart, rowEnd));
    Kokkos::Experimental::sort_thread(t, rowEntries);
  }

  Ordinal numRows;
  rowmap_t rowmap;
  entries_t entries;
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

template <typename Offset, typename Keys, typename Entries>
struct MaxScanFunctor {
  using value_type = uint64_t;

  MaxScanFunctor(uint64_t ncols_, const Keys& keys_, const Entries& entries_)
      : ncols(ncols_), keys(keys_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION
  void init(uint64_t& update) const { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(uint64_t& update, const uint64_t& input) const { update = Kokkos::max(update, input); }

  KOKKOS_INLINE_FUNCTION
  void operator()(Offset i, uint64_t& lmax, bool finalPass) const {
    lmax = Kokkos::max(lmax, keys(i));
    if (finalPass) {
      // lmax is the row containing entry i.
      // The key is equivalent to the entry's linear
      // index if the matrix were dense and row-major.
      keys(i) = lmax * ncols + entries(i);
    }
  }

  uint64_t ncols;
  Keys keys;
  Entries entries;
};

template <typename ExecSpace, typename Rowmap, typename Entries>
Kokkos::View<uint64_t*, ExecSpace> generateBulkCrsKeys(const ExecSpace& exec, const Rowmap& rowmap,
                                                       const Entries& entries,
                                                       typename Entries::non_const_value_type ncols) {
  using Offset    = typename Rowmap::non_const_value_type;
  using Ordinal   = typename Entries::non_const_value_type;
  Ordinal numRows = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
  Kokkos::View<uint64_t*, ExecSpace> keys("keys", entries.extent(0));
  Kokkos::parallel_for(
      "CRS bulk sorting: mark row begins", Kokkos::RangePolicy<ExecSpace>(exec, 0, numRows), KOKKOS_LAMBDA(Ordinal i) {
        Offset rowBegin = rowmap(i);
        // Only mark the beginnings of non-empty rows.
        // Otherwise multiple rows could try to update the same key.
        if (rowmap(i + 1) != rowBegin) {
          keys(rowBegin) = uint64_t(i);
        }
      });
  Kokkos::fence();
  Kokkos::parallel_scan("CRS bulk sorting: compute keys", Kokkos::RangePolicy<ExecSpace>(exec, 0, entries.extent(0)),
                        MaxScanFunctor<Offset, decltype(keys), Entries>(ncols, keys, entries));
  Kokkos::fence();
  return keys;
}

#ifndef KK_DISABLE_BULK_SORT_BY_KEY
template <typename ExecSpace, typename Rowmap, typename Entries>
Kokkos::View<typename Rowmap::non_const_value_type*, ExecSpace> computeEntryPermutation(
    const ExecSpace& exec, const Rowmap& rowmap, const Entries& entries, typename Entries::non_const_value_type ncols) {
  using Offset = typename Rowmap::non_const_value_type;
  auto keys    = generateBulkCrsKeys(exec, rowmap, entries, ncols);
  Kokkos::View<Offset*, ExecSpace> permutation(Kokkos::view_alloc(Kokkos::WithoutInitializing, "permutation"),
                                               entries.extent(0));
  // This initializes permutation as the identity
  KokkosKernels::Impl::sequential_fill(exec, permutation);
  Kokkos::Experimental::sort_by_key(exec, keys, permutation);
  return permutation;
}

// Heuristic for choosing bulk sorting algorithm
template <typename ExecSpace, typename Ordinal>
bool useBulkSortHeuristic(Ordinal avgDeg, Ordinal maxDeg) {
  // Issue 2352: the KokkosSparse::sort_crs_matrix uses Kokkos::Experimental::sort_by_key when this returns true.
  // sort_by_key executes on the host when a thrust-like library is not available, which really kills the performance in
  // a scenario where the bulk sort algorithm would otherwise be appropriate. Additionally, On MI300A, sorting via
  // ROCTHRUST was observed to be ~3x slower than the Kokkos kernels native implementation on some matrices of interest,
  // so on that architecture only always bypass bulk sort.
  // * GPU execution space, SYLC is enabled, but no ONEDPL does not have sort_by_key
  // * GPU execution space, HIP is enabled, but no ROCTHRUST
  // * GPU execution space, HIP is enabled, and GPU is GFX942
  // (Kokkos seems to require thrust when CUDA is enabled)
  if constexpr (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecSpace>()) {
#if (defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ONEDPL_HAS_SORT_BY_KEY)) || \
    (defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_ROCTHRUST)) ||        \
    (defined(KOKKOS_ENABLE_HIP) && defined(KOKKOS_ARCH_AMD_GFX942))
    return false;
#else
    // Use bulk sort if matrix is highly imbalanced,
    // OR the longest rows have many entries.
    return (maxDeg / 10 > avgDeg) || (maxDeg > 1024);
#endif
  } else {
    // Use bulk sort if matrix is highly imbalanced,
    // OR the longest rows have many entries.
    return (maxDeg / 10 > avgDeg) || (maxDeg > 1024);
  }
}
#endif

template <typename ExecSpace, typename Permutation, typename InView, typename OutView>
void applyPermutation(const ExecSpace& exec, const Permutation& permutation, const InView& in, const OutView& out) {
  Kokkos::parallel_for(
      "CRS bulk sorting: permute", Kokkos::RangePolicy<ExecSpace>(exec, 0, in.extent(0)),
      KOKKOS_LAMBDA(size_t i) { out(i) = in(permutation(i)); });
}

template <typename ExecSpace, typename Permutation, typename InView, typename OutView, typename Ordinal>
void applyPermutationBlockValues(const ExecSpace& exec, const Permutation& permutation, const InView& in,
                                 const OutView& out, Ordinal blockSize) {
  uint64_t scalarsPerBlock = (uint64_t)blockSize * blockSize;
  if (in.extent(0) % scalarsPerBlock)
    throw std::invalid_argument(
        "sort_bsr_matrix: matrix values extent not divisible by graph entries "
        "extent");
  Kokkos::parallel_for(
      "BSR bulk sorting: permute", Kokkos::RangePolicy<ExecSpace>(exec, 0, in.extent(0)), KOKKOS_LAMBDA(size_t i) {
        uint64_t blockIndex    = i / scalarsPerBlock;
        uint64_t offsetInBlock = i % scalarsPerBlock;
        out(i)                 = in(permutation(blockIndex) * scalarsPerBlock + offsetInBlock);
      });
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
