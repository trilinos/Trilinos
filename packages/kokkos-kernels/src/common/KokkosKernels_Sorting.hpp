/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef _KOKKOSKERNELS_SORTING_HPP
#define _KOKKOSKERNELS_SORTING_HPP

#include "Kokkos_Core.hpp"
#include "KokkosKernels_SimpleUtils.hpp"  //for kk_exclusive_parallel_prefix_sum
#include "KokkosKernels_ExecSpaceUtils.hpp"  //for kk_is_gpu_exec_space
#include <type_traits>

namespace KokkosKernels {

namespace Impl {
template <typename Value>
struct DefaultComparator {
  KOKKOS_INLINE_FUNCTION bool operator()(const Value lhs,
                                         const Value rhs) const {
    return lhs < rhs;
  }
};
}  // namespace Impl

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

// ----------------------------
// General device-level sorting
// ----------------------------

// Bitonic sort: sorts v according to the comparator object's operator().
// Default comparator is just operator< for v's element type.
template <
    typename View, typename ExecSpace, typename Ordinal,
    typename Comparator = Impl::DefaultComparator<typename View::value_type>>
void bitonicSort(View v, const Comparator& comp = Comparator());

// --------------------------------------------------------
// Serial sorting (callable inside any kernel or host code)
// --------------------------------------------------------

// Radix sort. Not in-place: requires scratch array 'valuesAux' to be the same
// size as values. ValueType must be an unsigned integer type.
template <typename Ordinal, typename ValueType>
KOKKOS_INLINE_FUNCTION void SerialRadixSort(ValueType* values,
                                            ValueType* valuesAux, Ordinal n);

// Same as SerialRadixSort, but also permutes perm[0...n] as it sorts
// values[0...n].
template <typename Ordinal, typename ValueType, typename PermType>
KOKKOS_INLINE_FUNCTION void SerialRadixSort2(ValueType* values,
                                             ValueType* valuesAux,
                                             PermType* perm, PermType* permAux,
                                             Ordinal n);

// -------------------------------------------------------------------
// Team-level parallel sorting (callable inside any TeamPolicy kernel)
// -------------------------------------------------------------------

// Comparison based sorting that uses the entire team (described by mem) to sort
// raw array according to the comparator.
template <typename Ordinal, typename ValueType, typename TeamMember,
          typename Comparator = Impl::DefaultComparator<ValueType>>
KOKKOS_INLINE_FUNCTION void TeamBitonicSort(
    ValueType* values, Ordinal n, const TeamMember mem,
    const Comparator& comp = Comparator());

// Same as SerialRadixSort, but also permutes perm[0...n] as it sorts
// values[0...n].
template <typename Ordinal, typename ValueType, typename PermType,
          typename TeamMember,
          typename Comparator = Impl::DefaultComparator<ValueType>>
KOKKOS_INLINE_FUNCTION void TeamBitonicSort2(
    ValueType* values, PermType* perm, Ordinal n, const TeamMember mem,
    const Comparator& comp = Comparator());

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

// Functor that sorts a view on one team
template <typename View, typename Ordinal, typename TeamMember,
          typename Comparator>
struct BitonicSingleTeamFunctor {
  BitonicSingleTeamFunctor(View& v_, const Comparator& comp_)
      : v(v_), comp(comp_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const {
    KokkosKernels::TeamBitonicSort<Ordinal, typename View::value_type,
                                   TeamMember, Comparator>(
        v.data(), v.extent(0), t, comp);
  };
  View v;
  Comparator comp;
};

// Functor that sorts equally sized chunks on each team
template <typename View, typename Ordinal, typename TeamMember,
          typename Comparator>
struct BitonicChunkFunctor {
  BitonicChunkFunctor(View& v_, const Comparator& comp_, Ordinal chunkSize_)
      : v(v_), comp(comp_), chunkSize(chunkSize_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const {
    Ordinal chunk      = t.league_rank();
    Ordinal chunkStart = chunk * chunkSize;
    Ordinal n          = chunkSize;
    if (chunkStart + n > Ordinal(v.extent(0))) n = v.extent(0) - chunkStart;
    KokkosKernels::TeamBitonicSort<Ordinal, typename View::value_type,
                                   TeamMember, Comparator>(
        v.data() + chunkStart, n, t, comp);
  };
  View v;
  Comparator comp;
  Ordinal chunkSize;
};

// Functor that does just the first phase (brown) of bitonic sort on
// equally-sized chunks
template <typename View, typename Ordinal, typename TeamMember,
          typename Comparator>
struct BitonicPhase1Functor {
  typedef typename View::value_type Value;
  BitonicPhase1Functor(View& v_, const Comparator& comp_, Ordinal boxSize_,
                       Ordinal teamsPerBox_)
      : v(v_), comp(comp_), boxSize(boxSize_), teamsPerBox(teamsPerBox_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const {
    Ordinal box         = t.league_rank() / teamsPerBox;
    Ordinal boxStart    = boxSize * box;
    Ordinal work        = boxSize / teamsPerBox / 2;
    Ordinal workStart   = work * (t.league_rank() % teamsPerBox);
    Ordinal workReflect = boxSize - workStart - 1;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(t, work),
                         [&](const Ordinal i) {
                           Ordinal elem1 = boxStart + workStart + i;
                           Ordinal elem2 = boxStart + workReflect - i;
                           if (elem2 < Ordinal(v.extent(0))) {
                             if (comp(v(elem2), v(elem1))) {
                               Value temp = v(elem1);
                               v(elem1)   = v(elem2);
                               v(elem2)   = temp;
                             }
                           }
                         });
  };
  View v;
  Comparator comp;
  Ordinal boxSize;
  Ordinal teamsPerBox;
};

// Functor that does the second phase (red) of bitonic sort
template <typename View, typename Ordinal, typename TeamMember,
          typename Comparator>
struct BitonicPhase2Functor {
  typedef typename View::value_type Value;
  BitonicPhase2Functor(View& v_, const Comparator& comp_, Ordinal boxSize_,
                       Ordinal teamsPerBox_)
      : v(v_), comp(comp_), boxSize(boxSize_), teamsPerBox(teamsPerBox_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const {
    Ordinal logBoxSize = 1;
    while ((Ordinal(1) << logBoxSize) < boxSize) logBoxSize++;
    Ordinal box       = t.league_rank() / teamsPerBox;
    Ordinal boxStart  = boxSize * box;
    Ordinal work      = boxSize / teamsPerBox / 2;
    Ordinal workStart = boxStart + work * (t.league_rank() % teamsPerBox);
    Ordinal jump      = boxSize / 2;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(t, work),
                         [&](const Ordinal i) {
                           Ordinal elem1 = workStart + i;
                           Ordinal elem2 = workStart + jump + i;
                           if (elem2 < Ordinal(v.extent(0))) {
                             if (comp(v(elem2), v(elem1))) {
                               Value temp = v(elem1);
                               v(elem1)   = v(elem2);
                               v(elem2)   = temp;
                             }
                           }
                         });
    if (teamsPerBox == 1) {
      // This team can finish phase 2 for all the smaller red boxes that follow,
      // since there are no longer cross-team data dependencies
      for (Ordinal subLevel = 1; subLevel < logBoxSize; subLevel++) {
        t.team_barrier();
        Ordinal logSubBoxSize = logBoxSize - subLevel;
        Ordinal subBoxSize    = Ordinal(1) << logSubBoxSize;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(t, work), [&](const Ordinal i) {
              Ordinal globalThread = i + t.league_rank() * work;
              Ordinal subBox       = globalThread >> (logSubBoxSize - 1);
              Ordinal subBoxStart  = subBox << logSubBoxSize;
              Ordinal subBoxOffset =
                  globalThread & ((Ordinal(1) << (logSubBoxSize - 1)) -
                                  1);  // i % (subBoxSize / 2)
              Ordinal elem1 = subBoxStart + subBoxOffset;
              // later phases (pink box): within a block, compare with fixed
              // distance (boxSize / 2) apart
              Ordinal elem2 = elem1 + subBoxSize / 2;
              if (elem2 < Ordinal(v.extent(0))) {
                if (comp(v(elem2), v(elem1))) {
                  Value temp = v(elem1);
                  v(elem1)   = v(elem2);
                  v(elem2)   = temp;
                }
              }
            });
      }
    }
  };
  View v;
  Comparator comp;
  Ordinal boxSize;
  Ordinal teamsPerBox;
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
  bool useRadix  = !Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows  = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
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

// Sort a CRS graph: within each row, sort entries ascending by column.
template <typename execution_space, typename rowmap_t, typename entries_t>
void sort_crs_graph(const rowmap_t& rowmap, const entries_t& entries) {
  using lno_t    = typename entries_t::non_const_value_type;
  using team_pol = Kokkos::TeamPolicy<execution_space>;
  bool useRadix  = !Impl::kk_is_gpu_exec_space<execution_space>();
  lno_t numRows  = rowmap.extent(0) ? rowmap.extent(0) - 1 : 0;
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
  Impl::kk_exclusive_parallel_prefix_sum<rowmap_t, exec_space>(A.numRows() + 1,
                                                               mergedRowmap);
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
  Impl::kk_exclusive_parallel_prefix_sum<rowmap_t, exec_space>(numRows + 1,
                                                               rowmap_out);
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

// Version to be called from host on a single array
// Generally ~2x slower than Kokkos::sort() for large arrays (> 50 M elements),
// but faster for smaller arrays.
//
// This is more general than BinSort: bitonic supports any trivially copyable
// type and an arbitrary device-compatible comparison operator (provided through
// operator() of Comparator) If comparator is void, use operator< (which should
// only be used for primitives)
template <typename View, typename ExecSpace, typename Ordinal,
          typename Comparator>
void bitonicSort(View v, const Comparator& comp) {
  typedef Kokkos::TeamPolicy<ExecSpace> team_policy;
  typedef typename team_policy::member_type team_member;
  Ordinal n = v.extent(0);
  // If n is small, just sort on a single team
  if (n <= Ordinal(1) << 12) {
    Kokkos::parallel_for(
        team_policy(1, Kokkos::AUTO()),
        Impl::BitonicSingleTeamFunctor<View, Ordinal, team_member, Comparator>(
            v, comp));
  } else {
    Ordinal npot = 1;
    while (npot < n) npot <<= 1;
    // Partition the data equally among fixed number of teams
    Ordinal chunkSize = 512;
    Ordinal numTeams  = npot / chunkSize;
    // First, sort within teams
    Kokkos::parallel_for(
        team_policy(numTeams, Kokkos::AUTO()),
        Impl::BitonicChunkFunctor<View, Ordinal, team_member, Comparator>(
            v, comp, chunkSize));
    for (int teamsPerBox = 2; teamsPerBox <= npot / chunkSize;
         teamsPerBox *= 2) {
      Ordinal boxSize = teamsPerBox * chunkSize;
      Kokkos::parallel_for(
          team_policy(numTeams, Kokkos::AUTO()),
          Impl::BitonicPhase1Functor<View, Ordinal, team_member, Comparator>(
              v, comp, boxSize, teamsPerBox));
      for (int boxDiv = 1; teamsPerBox >> boxDiv; boxDiv++) {
        Kokkos::parallel_for(
            team_policy(numTeams, Kokkos::AUTO()),
            Impl::BitonicPhase2Functor<View, Ordinal, team_member, Comparator>(
                v, comp, boxSize >> boxDiv, teamsPerBox >> boxDiv));
      }
    }
  }
}

// Radix sort for integers, on a single thread within a team.
// Pros: few diverging branches, so OK for sorting on a single GPU vector lane.
// Better on CPU cores. Con: requires auxiliary storage, and this version only
// works for integers
template <typename Ordinal, typename ValueType>
KOKKOS_INLINE_FUNCTION void SerialRadixSort(ValueType* values,
                                            ValueType* valuesAux, Ordinal n) {
  static_assert(
      std::is_integral<ValueType>::value && std::is_unsigned<ValueType>::value,
      "radixSort can only be run on unsigned integers.");
  if (n <= 1) return;
  ValueType maxVal = 0;
  for (Ordinal i = 0; i < n; i++) {
    if (maxVal < values[i]) maxVal = values[i];
  }
  // determine how many significant bits the data has
  int passes = 0;
  while (maxVal) {
    maxVal >>= 4;
    passes++;
  }
  // Is the data currently held in values (false) or valuesAux (true)?
  bool inAux = false;
  // sort 4 bits at a time, into 16 buckets
  ValueType mask = 0xF;
  // maskPos counts the low bit index of mask (0, 4, 8, ...)
  Ordinal maskPos = 0;
  for (int p = 0; p < passes; p++) {
    // Count the number of elements in each bucket
    Ordinal count[16] = {0};
    Ordinal offset[17];
    if (!inAux) {
      for (Ordinal i = 0; i < n; i++) {
        count[(values[i] & mask) >> maskPos]++;
      }
    } else {
      for (Ordinal i = 0; i < n; i++) {
        count[(valuesAux[i] & mask) >> maskPos]++;
      }
    }
    offset[0] = 0;
    // get offset as the prefix sum for count
    for (Ordinal i = 0; i < 16; i++) {
      offset[i + 1] = offset[i] + count[i];
    }
    // now for each element in [lo, hi), move it to its offset in the other
    // buffer this branch should be ok because whichBuf is the same on all
    // threads
    if (!inAux) {
      for (Ordinal i = 0; i < n; i++) {
        Ordinal bucket = (values[i] & mask) >> maskPos;
        valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
        count[bucket]--;
      }
    } else {
      for (Ordinal i = 0; i < n; i++) {
        Ordinal bucket = (valuesAux[i] & mask) >> maskPos;
        values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
        count[bucket]--;
      }
    }
    inAux = !inAux;
    mask  = mask << 4;
    maskPos += 4;
  }
  // Move values back into main array if they are currently in aux.
  // This is the case if an odd number of rounds were done.
  if (inAux) {
    for (Ordinal i = 0; i < n; i++) {
      values[i] = valuesAux[i];
    }
  }
}

// Radix sort for integers (no internal parallelism).
// While sorting, also permute "perm" array along with the values.
// Pros: few diverging branches, so good for sorting on a single GPU vector
// lane. Con: requires auxiliary storage, this version only works for integers
// (although float/double is possible)
template <typename Ordinal, typename ValueType, typename PermType>
KOKKOS_INLINE_FUNCTION void SerialRadixSort2(ValueType* values,
                                             ValueType* valuesAux,
                                             PermType* perm, PermType* permAux,
                                             Ordinal n) {
  static_assert(
      std::is_integral<ValueType>::value && std::is_unsigned<ValueType>::value,
      "radixSort can only be run on unsigned integers.");
  if (n <= 1) return;
  ValueType maxVal = 0;
  for (Ordinal i = 0; i < n; i++) {
    if (maxVal < values[i]) maxVal = values[i];
  }
  int passes = 0;
  while (maxVal) {
    maxVal >>= 4;
    passes++;
  }
  // Is the data currently held in values (false) or valuesAux (true)?
  bool inAux = false;
  // sort 4 bits at a time, into 16 buckets
  ValueType mask = 0xF;
  // maskPos counts the low bit index of mask (0, 4, 8, ...)
  Ordinal maskPos = 0;
  for (int p = 0; p < passes; p++) {
    // Count the number of elements in each bucket
    Ordinal count[16] = {0};
    Ordinal offset[17];
    if (!inAux) {
      for (Ordinal i = 0; i < n; i++) {
        count[(values[i] & mask) >> maskPos]++;
      }
    } else {
      for (Ordinal i = 0; i < n; i++) {
        count[(valuesAux[i] & mask) >> maskPos]++;
      }
    }
    offset[0] = 0;
    // get offset as the prefix sum for count
    for (Ordinal i = 0; i < 16; i++) {
      offset[i + 1] = offset[i] + count[i];
    }
    // now for each element in [lo, hi), move it to its offset in the other
    // buffer this branch should be ok because whichBuf is the same on all
    // threads
    if (!inAux) {
      for (Ordinal i = 0; i < n; i++) {
        Ordinal bucket = (values[i] & mask) >> maskPos;
        valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
        permAux[offset[bucket + 1] - count[bucket]]   = perm[i];
        count[bucket]--;
      }
    } else {
      for (Ordinal i = 0; i < n; i++) {
        Ordinal bucket = (valuesAux[i] & mask) >> maskPos;
        values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
        perm[offset[bucket + 1] - count[bucket]]   = permAux[i];
        count[bucket]--;
      }
    }
    inAux = !inAux;
    mask  = mask << 4;
    maskPos += 4;
  }
  // Move values back into main array if they are currently in aux.
  // This is the case if an odd number of rounds were done.
  if (inAux) {
    for (Ordinal i = 0; i < n; i++) {
      values[i] = valuesAux[i];
      perm[i]   = permAux[i];
    }
  }
}

// Bitonic merge sort (requires only comparison operators and
// trivially-copyable) Pros: In-place, plenty of parallelism for GPUs, and
// memory references are coalesced Con: O(n log^2(n)) serial time is bad on CPUs
// Good diagram of the algorithm at https://en.wikipedia.org/wiki/Bitonic_sorter
template <typename Ordinal, typename ValueType, typename TeamMember,
          typename Comparator>
KOKKOS_INLINE_FUNCTION void TeamBitonicSort(ValueType* values, Ordinal n,
                                            const TeamMember mem,
                                            const Comparator& comp) {
  // Algorithm only works on power-of-two input size only.
  // If n is not a power-of-two, will implicitly pretend
  // that values[i] for i >= n is just the max for ValueType, so it never gets
  // swapped
  Ordinal npot   = 1;
  Ordinal levels = 0;
  while (npot < n) {
    levels++;
    npot <<= 1;
  }
  for (Ordinal i = 0; i < levels; i++) {
    for (Ordinal j = 0; j <= i; j++) {
      // n/2 pairs of items are compared in parallel
      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(mem, npot / 2), [=](const Ordinal t) {
            // How big are the brown/pink boxes?
            Ordinal boxSize = Ordinal(2) << (i - j);
            // Which box contains this thread?
            Ordinal boxID     = t >> (i - j);          // t * 2 / boxSize;
            Ordinal boxStart  = boxID << (1 + i - j);  // boxID * boxSize
            Ordinal boxOffset = t - (boxStart >> 1);   // t - boxID * boxSize /
                                                       // 2;
            Ordinal elem1 = boxStart + boxOffset;
            if (j == 0) {
              // first phase (brown box): within a block, compare with the
              // opposite value in the box
              Ordinal elem2 = boxStart + boxSize - 1 - boxOffset;
              if (elem2 < n) {
                // both elements in bounds, so compare them and swap if out of
                // order
                if (comp(values[elem2], values[elem1])) {
                  ValueType temp = values[elem1];
                  values[elem1]  = values[elem2];
                  values[elem2]  = temp;
                }
              }
            } else {
              // later phases (pink box): within a block, compare with fixed
              // distance (boxSize / 2) apart
              Ordinal elem2 = elem1 + boxSize / 2;
              if (elem2 < n) {
                if (comp(values[elem2], values[elem1])) {
                  ValueType temp = values[elem1];
                  values[elem1]  = values[elem2];
                  values[elem2]  = temp;
                }
              }
            }
          });
      mem.team_barrier();
    }
  }
}

// Sort "values", while applying the same swaps to "perm"
template <typename Ordinal, typename ValueType, typename PermType,
          typename TeamMember, typename Comparator>
KOKKOS_INLINE_FUNCTION void TeamBitonicSort2(ValueType* values, PermType* perm,
                                             Ordinal n, const TeamMember mem,
                                             const Comparator& comp) {
  // Algorithm only works on power-of-two input size only.
  // If n is not a power-of-two, will implicitly pretend
  // that values[i] for i >= n is just the max for ValueType, so it never gets
  // swapped
  Ordinal npot   = 1;
  Ordinal levels = 0;
  while (npot < n) {
    levels++;
    npot <<= 1;
  }
  for (Ordinal i = 0; i < levels; i++) {
    for (Ordinal j = 0; j <= i; j++) {
      // n/2 pairs of items are compared in parallel
      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(mem, npot / 2), [=](const Ordinal t) {
            // How big are the brown/pink boxes?
            Ordinal boxSize = Ordinal(2) << (i - j);
            // Which box contains this thread?
            Ordinal boxID     = t >> (i - j);          // t * 2 / boxSize;
            Ordinal boxStart  = boxID << (1 + i - j);  // boxID * boxSize
            Ordinal boxOffset = t - (boxStart >> 1);   // t - boxID * boxSize /
                                                       // 2;
            Ordinal elem1 = boxStart + boxOffset;
            if (j == 0) {
              // first phase (brown box): within a block, compare with the
              // opposite value in the box
              Ordinal elem2 = boxStart + boxSize - 1 - boxOffset;
              if (elem2 < n) {
                // both elements in bounds, so compare them and swap if out of
                // order
                if (comp(values[elem2], values[elem1])) {
                  ValueType temp1 = values[elem1];
                  values[elem1]   = values[elem2];
                  values[elem2]   = temp1;
                  PermType temp2  = perm[elem1];
                  perm[elem1]     = perm[elem2];
                  perm[elem2]     = temp2;
                }
              }
            } else {
              // later phases (pink box): within a block, compare with fixed
              // distance (boxSize / 2) apart
              Ordinal elem2 = elem1 + boxSize / 2;
              if (elem2 < n) {
                if (comp(values[elem2], values[elem1])) {
                  ValueType temp1 = values[elem1];
                  values[elem1]   = values[elem2];
                  values[elem2]   = temp1;
                  PermType temp2  = perm[elem1];
                  perm[elem1]     = perm[elem2];
                  perm[elem2]     = temp2;
                }
              }
            }
          });
      mem.team_barrier();
    }
  }
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

template <
    typename View, typename ExecSpace, typename Ordinal,
    typename Comparator = Impl::DefaultComparator<typename View::value_type>>
[[deprecated]] void bitonicSort(View v, const Comparator& comp = Comparator()) {
  KokkosKernels::bitonicSort<View, ExecSpace, Ordinal, Comparator>(v, comp);
}

template <typename Ordinal, typename ValueType>
[[deprecated]] KOKKOS_INLINE_FUNCTION void SerialRadixSort(ValueType* values,
                                                           ValueType* valuesAux,
                                                           Ordinal n) {
  KokkosKernels::SerialRadixSort<Ordinal, ValueType>(values, valuesAux, n);
}

// Same as SerialRadixSort, but also permutes perm[0...n] as it sorts
// values[0...n].
template <typename Ordinal, typename ValueType, typename PermType>
[[deprecated]] KOKKOS_INLINE_FUNCTION void SerialRadixSort2(
    ValueType* values, ValueType* valuesAux, PermType* perm, PermType* permAux,
    Ordinal n) {
  KokkosKernels::SerialRadixSort2<Ordinal, ValueType, PermType>(
      values, valuesAux, perm, permAux, n);
}

template <typename Ordinal, typename ValueType, typename TeamMember,
          typename Comparator = Impl::DefaultComparator<ValueType>>
[[deprecated]] KOKKOS_INLINE_FUNCTION void TeamBitonicSort(
    ValueType* values, Ordinal n, const TeamMember mem,
    const Comparator& comp = Comparator()) {
  KokkosKernels::TeamBitonicSort<Ordinal, ValueType, TeamMember, Comparator>(
      values, n, mem, comp);
}

// Same as SerialRadixSort, but also permutes perm[0...n] as it sorts
// values[0...n].
template <typename Ordinal, typename ValueType, typename PermType,
          typename TeamMember,
          typename Comparator = Impl::DefaultComparator<ValueType>>
[[deprecated]] KOKKOS_INLINE_FUNCTION void TeamBitonicSort2(
    ValueType* values, PermType* perm, Ordinal n, const TeamMember mem,
    const Comparator& comp = Comparator()) {
  KokkosKernels::TeamBitonicSort2<Ordinal, ValueType, PermType, TeamMember,
                                  Comparator>(values, perm, n, mem, comp);
}
}  // namespace Impl

}  // namespace KokkosKernels

#endif
