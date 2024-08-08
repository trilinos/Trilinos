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

#ifndef _KOKKOS_SPADD_SYMBOLIC_IMPL_HPP
#define _KOKKOS_SPADD_SYMBOLIC_IMPL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosSparse {
namespace Impl {

// Two types are the same (ignoring const)
template <typename T, typename U>
constexpr bool spadd_symbolic_same_type =
    std::is_same_v<typename std::remove_const_t<T>, typename std::remove_const_t<U>>;

// get C rowmap for sorted input
template <typename size_type, typename ordinal_type, typename ARowPtrsT, typename BRowPtrsT, typename AColIndsT,
          typename BColIndsT, typename CRowPtrsT, typename ExecSpace>
struct SortedCountEntriesRange {
  SortedCountEntriesRange(ordinal_type nrows_, const typename ARowPtrsT::const_type& Arowptrs_,
                          const AColIndsT& Acolinds_, const typename BRowPtrsT::const_type& Browptrs_,
                          const BColIndsT& Bcolinds_, const CRowPtrsT& Crowcounts_)
      : nrows(nrows_),
        Arowptrs(Arowptrs_),
        Acolinds(Acolinds_),
        Browptrs(Browptrs_),
        Bcolinds(Bcolinds_),
        Crowcounts(Crowcounts_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    const ordinal_type ORDINAL_MAX = Kokkos::ArithTraits<ordinal_type>::max();

    // count the union of nonzeros in Arow and Brow
    size_type numEntries = 0;
    size_type ai         = 0;
    size_type bi         = 0;
    size_type Arowstart  = Arowptrs(i);
    size_type Arowlen    = Arowptrs(i + 1) - Arowstart;
    size_type Browstart  = Browptrs(i);
    size_type Browlen    = Browptrs(i + 1) - Browstart;
    ordinal_type Acol    = (Arowlen == 0) ? ORDINAL_MAX : Acolinds(Arowstart);
    ordinal_type Bcol    = (Browlen == 0) ? ORDINAL_MAX : Bcolinds(Browstart);
    while (Acol != ORDINAL_MAX || Bcol != ORDINAL_MAX) {
      ordinal_type Ccol = (Acol < Bcol) ? Acol : Bcol;
      numEntries++;
      // Eat all entries in both A and B which have this column
      // This also results in Acol/Bcol being updated to following entries for
      // next loop iter
      while (Acol == Ccol) Acol = (ai == Arowlen) ? ORDINAL_MAX : Acolinds(Arowstart + ai++);
      while (Bcol == Ccol) Bcol = (bi == Browlen) ? ORDINAL_MAX : Bcolinds(Browstart + bi++);
    }
    Crowcounts(i) = numEntries;
  }

  ordinal_type nrows;
  const typename ARowPtrsT::const_type Arowptrs;
  const AColIndsT Acolinds;
  const typename BRowPtrsT::const_type Browptrs;
  const BColIndsT Bcolinds;
  CRowPtrsT Crowcounts;
};

template <typename size_type, typename ordinal_type, typename ARowPtrsT, typename BRowPtrsT, typename AColIndsT,
          typename BColIndsT, typename CRowPtrsT, typename ExecSpace>
struct SortedCountEntriesTeam {
  SortedCountEntriesTeam(ordinal_type nrows_, const typename ARowPtrsT::const_type& Arowptrs_,
                         const AColIndsT& Acolinds_, const typename BRowPtrsT::const_type& Browptrs_,
                         const BColIndsT& Bcolinds_, const CRowPtrsT& Crowcounts_)
      : nrows(nrows_),
        Arowptrs(Arowptrs_),
        Acolinds(Acolinds_),
        Browptrs(Browptrs_),
        Bcolinds(Bcolinds_),
        Crowcounts(Crowcounts_) {}

  using TeamPol = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMem = typename TeamPol::member_type;

  KOKKOS_INLINE_FUNCTION void longRowFallback(const ordinal_type i) const {
    const ordinal_type ORDINAL_MAX = Kokkos::ArithTraits<ordinal_type>::max();

    // count the union of nonzeros in Arow and Brow
    size_type numEntries = 0;
    size_type ai         = 0;
    size_type bi         = 0;
    size_type Arowstart  = Arowptrs(i);
    size_type Arowlen    = Arowptrs(i + 1) - Arowstart;
    size_type Browstart  = Browptrs(i);
    size_type Browlen    = Browptrs(i + 1) - Browstart;
    ordinal_type Acol    = (Arowlen == 0) ? ORDINAL_MAX : Acolinds(Arowstart);
    ordinal_type Bcol    = (Browlen == 0) ? ORDINAL_MAX : Bcolinds(Browstart);
    while (Acol != ORDINAL_MAX || Bcol != ORDINAL_MAX) {
      ordinal_type Ccol = (Acol < Bcol) ? Acol : Bcol;
      numEntries++;
      // Eat all entries in both A and B which have this column
      // This also results in Acol/Bcol being updated to following entries for
      // next loop iter
      while (Acol == Ccol) Acol = (ai == Arowlen) ? ORDINAL_MAX : Acolinds(Arowstart + ai++);
      while (Bcol == Ccol) Bcol = (bi == Browlen) ? ORDINAL_MAX : Bcolinds(Browstart + bi++);
    }
    Crowcounts(i) = numEntries;
  }

  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem t) const {
    ordinal_type i = t.league_rank() * t.team_size() + t.team_rank();
    if (i >= nrows) return;
    ordinal_type* allScratch = (ordinal_type*)t.team_shmem().get_shmem(totalShared);
    ordinal_type* scratch    = allScratch + t.team_rank() * sharedPerThread;
    ordinal_type Arowstart   = Arowptrs(i);
    ordinal_type Arowlen     = Arowptrs(i + 1) - Arowstart;
    ordinal_type Browstart   = Browptrs(i);
    ordinal_type Browlen     = Browptrs(i + 1) - Browstart;
    ordinal_type n           = Arowlen + Browlen;
    if (n > sharedPerThread) {
      // fall back to slow serial method
      Kokkos::single(Kokkos::PerThread(t), [&]() { longRowFallback(i); });
      return;
    }
    if (n == 0) {
      Kokkos::single(Kokkos::PerThread(t), [&]() { Crowcounts(i) = 0; });
      return;
    }
    // Figure out the number of bitonic steps: ceil(log2(n))
    ordinal_type npot   = 1;
    ordinal_type levels = 0;
    while (npot < n) {
      levels++;
      npot <<= 1;
    }
    // Copy A and B entries to scratch
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, Arowlen),
                         [&](ordinal_type j) { scratch[j] = Acolinds(Arowstart + j); });
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, Browlen),
                         [&](ordinal_type j) { scratch[npot - 1 - j] = Bcolinds(Browstart + j); });
    // Fill space between A and B with ORDINAL_MAX,
    // to maintain a valid bitonic sequence of power-of-two length
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, npot - n),
                         [&](ordinal_type j) { scratch[Arowlen + j] = Kokkos::ArithTraits<ordinal_type>::max(); });
    // npot = 2^levels
    for (ordinal_type level = 0; level < levels; level++) {
      // npot/2 pairs of items are compared in parallel
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, npot >> 1), [&](const ordinal_type j) {
        ordinal_type boxSize = npot >> level;
        // Which box contains this thread?
        // box = (j / boxSize), and boxSize =
        // 2^(levels-level), so box = j * 2^(level-levels)
        // = j >> (levels - level)
        ordinal_type boxID = (j * 2) >> (levels - level);
        // boxStart = boxID * boxSize = boxID *
        // 2^(levels-level) = boxID << (levels-level)
        ordinal_type boxStart  = boxID << (levels - level);
        ordinal_type boxOffset = j - boxID * boxSize / 2;
        ordinal_type elem1     = boxStart + boxOffset;
        ordinal_type elem2     = elem1 + (boxSize >> 1);
        if (scratch[elem2] < scratch[elem1]) {
          ordinal_type temp = scratch[elem1];
          scratch[elem1]    = scratch[elem2];
          scratch[elem2]    = temp;
        }
      });
    }
    // Finally, count the number of distinct entries (this is #rising edges + 1)
    ordinal_type risingEdges;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(t, n - 1),
        [&](const ordinal_type j, ordinal_type& lcount) {
          if (scratch[j] != scratch[j + 1]) lcount++;
        },
        risingEdges);
    Kokkos::single(Kokkos::PerThread(t), [&]() { Crowcounts(i) = risingEdges + 1; });
  }

  size_t team_shmem_size(int teamSize) const { return sharedPerThread * sizeof(ordinal_type) * teamSize; }

  ordinal_type nrows;
  const typename ARowPtrsT::const_type Arowptrs;
  const AColIndsT Acolinds;
  const typename BRowPtrsT::const_type Browptrs;
  const BColIndsT Bcolinds;
  CRowPtrsT Crowcounts;
  int sharedPerThread;  // Shared for each thread, measured in
                        // sizeof(ordinal_type)
  int totalShared;      // Shared for whole team, measured in bytes
};

// get upper bound for C entries per row (assumes worst case, that entries in A
// and B on each row are disjoint)
template <typename size_type, typename ordinal_type, typename ARowPtrsT, typename BRowPtrsT, typename CRowPtrsT>
struct UnsortedEntriesUpperBound {
  UnsortedEntriesUpperBound(ordinal_type nrows_, const typename ARowPtrsT::const_type& Arowptrs_,
                            const typename BRowPtrsT::const_type& Browptrs_, const CRowPtrsT& Crowcounts_)
      : nrows(nrows_), Arowptrs(Arowptrs_), Browptrs(Browptrs_), Crowcounts(Crowcounts_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    Crowcounts(i) = (Arowptrs(i + 1) - Arowptrs(i)) + (Browptrs(i + 1) - Browptrs(i));
    if (i == nrows - 1) {
      // last workitem also zeros the one-past-end entry of row counts, so
      // that prefix sum is correct
      Crowcounts(nrows) = 0;
    }
  }
  ordinal_type nrows;
  const typename ARowPtrsT::const_type Arowptrs;
  const typename BRowPtrsT::const_type Browptrs;
  CRowPtrsT Crowcounts;
};

// Unsorted symbolic: new functors:
//  -compute uncompressed C (entries only, no values)
//  -sort uncompressed C entries within row, while permuting A union B
//  permutation array -compress sorted C entries and A,B perm arrays at the same
//  time, which produces Crowcounts value
// Inputs: A, B rowptrs/colinds, C uncompressed rowptrs (and allocated C
// entries) Output: C uncompressed colinds
template <typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
          typename AcolindsT, typename BcolindsT, typename CcolindsT>
struct UnmergedSumFunctor {
  UnmergedSumFunctor(ordinal_type nrows_, const ArowptrsT& Arowptrs_, const AcolindsT& Acolinds_,
                     const BrowptrsT& Browptrs_, const BcolindsT& Bcolinds_, const CrowptrsT& Crowptrs_,
                     const CcolindsT& Ccolinds_, const CcolindsT& ABperm_)
      : nrows(nrows_),
        Arowptrs(Arowptrs_),
        Acolinds(Acolinds_),
        Browptrs(Browptrs_),
        Bcolinds(Bcolinds_),
        Crowptrs(Crowptrs_),
        Ccolinds(Ccolinds_),
        ABperm(ABperm_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    size_type inserted  = 0;
    size_type crowstart = Crowptrs(i);
    size_type arowstart = Arowptrs(i);
    size_type arowlen   = Arowptrs(i + 1) - arowstart;
    size_type browstart = Browptrs(i);
    size_type browlen   = Browptrs(i + 1) - browstart;
    // Insert all A entries, then all B entries
    for (size_type j = 0; j < arowlen; j++) {
      Ccolinds(crowstart + inserted) = Acolinds(arowstart + j);
      ABperm(crowstart + inserted)   = j;
      inserted++;
    }
    for (size_type j = 0; j < browlen; j++) {
      Ccolinds(crowstart + inserted) = Bcolinds(browstart + j);
      // tell A and B permutation values apart by adding arowlen as a bias to B
      // values
      ABperm(crowstart + inserted) = j + arowlen;
      inserted++;
    }
  }
  ordinal_type nrows;
  const ArowptrsT Arowptrs;
  const AcolindsT Acolinds;
  const BrowptrsT Browptrs;
  const BcolindsT Bcolinds;
  const CrowptrsT Crowptrs;
  CcolindsT Ccolinds;
  CcolindsT ABperm;
};

template <typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
          typename CcolindsT, typename OffsetView>
struct MergeEntriesFunctor {
  MergeEntriesFunctor(ordinal_type nrows_, const ArowptrsT& Arowptrs_, const BrowptrsT& Browptrs_,
                      const OffsetView& Crowptrs_, const CrowptrsT& Crowcounts_, const CcolindsT& Ccolinds_,
                      const CcolindsT& ABperm_, const CcolindsT& Apos_, const CcolindsT& Bpos_)
      : nrows(nrows_),
        Arowptrs(Arowptrs_),
        Browptrs(Browptrs_),
        Crowptrs(Crowptrs_),
        Crowcounts(Crowcounts_),
        Ccolinds(Ccolinds_),
        ABperm(ABperm_),
        Apos(Apos_),
        Bpos(Bpos_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    size_type CrowStart = Crowptrs(i);
    size_type CrowEnd   = Crowptrs(i + 1);
    if (CrowEnd == CrowStart) {
      Crowcounts(i) = 0;
      return;
    }
    size_type ArowStart = Arowptrs(i);
    size_type ArowNum   = Arowptrs(i + 1) - ArowStart;
    size_type BrowStart = Browptrs(i);
    ordinal_type CFit   = 0;  // counting through merged C indices (within row)
    for (size_type Cit = CrowStart; Cit < CrowEnd; Cit++) {
      if ((Cit > CrowStart) && (Ccolinds(Cit) != Ccolinds(Cit - 1))) {
        // This is a different column than the previous entry, and is not the
        // first entry. This means that this is the first occurence of a unique
        // column.
        CFit++;
      }
      size_type permVal = ABperm(Cit);
      if (permVal < ArowNum) {
        // Entry belongs to A
        ordinal_type Aindex = permVal;
        // The Aindex'th entry in row i of A will be added into the CFit'th
        // entry in C
        Apos(ArowStart + Aindex) = CFit;
      } else {
        // Entry belongs to B
        ordinal_type Bindex = permVal - ArowNum;
        // The Bindex'th entry in row i of B will be added into the CFit'th
        // entry in C
        Bpos(BrowStart + Bindex) = CFit;
      }
    }
    // At end of the row, know how many entries are in merged C.
    // Right now, CFit is the index of the last Apos/Bpos,
    // so adding one gives the total number of entries.
    Crowcounts(i) = CFit + 1;
  }
  ordinal_type nrows;
  const ArowptrsT Arowptrs;
  const BrowptrsT Browptrs;
  const OffsetView Crowptrs;
  CrowptrsT Crowcounts;
  CcolindsT Ccolinds;
  const CcolindsT ABperm;
  CcolindsT Apos;
  CcolindsT Bpos;
};

// Run SortedCountEntries: non-GPU, always uses the RangePolicy version.
template <typename execution_space, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename clno_row_view_t_>
void runSortedCountEntries(
    const execution_space& exec, const alno_row_view_t_& a_rowmap, const alno_nnz_view_t_& a_entries,
    const blno_row_view_t_& b_rowmap, const blno_nnz_view_t_& b_entries, const clno_row_view_t_& c_rowmap,
    typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;
  using range_type   = Kokkos::RangePolicy<execution_space>;
  auto nrows         = c_rowmap.extent(0) - 1;
  SortedCountEntriesRange<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, alno_nnz_view_t_,
                          blno_nnz_view_t_, clno_row_view_t_, execution_space>
      countEntries(nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap);
  Kokkos::parallel_for("KokkosSparse::SpAdd::Symbolic::InputSorted::CountEntries", range_type(exec, 0, nrows),
                       countEntries);
}

// Run SortedCountEntries: GPU, uses the TeamPolicy or RangePolicy depending
//  on average nz per row (a runtime decision)
template <typename execution_space, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename clno_row_view_t_>
void runSortedCountEntries(
    const execution_space& exec, const alno_row_view_t_& a_rowmap, const alno_nnz_view_t_& a_entries,
    const blno_row_view_t_& b_rowmap, const blno_nnz_view_t_& b_entries, const clno_row_view_t_& c_rowmap,
    typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr) {
  using size_type     = typename KernelHandle::size_type;
  using ordinal_type  = typename KernelHandle::nnz_lno_t;
  using RangePol      = Kokkos::RangePolicy<execution_space>;
  using TeamPol       = Kokkos::TeamPolicy<execution_space>;
  auto nrows          = c_rowmap.extent(0) - 1;
  size_type c_est_nnz = 1.4 * (a_entries.extent(0) + b_entries.extent(0)) / nrows;
  if (c_est_nnz <= 512) {
    // Convert c_est_nnz to a power of 2
    size_type pot_est_nnz = 1;
    while (pot_est_nnz < c_est_nnz) pot_est_nnz *= 2;
    // Estimate max number of uncompressed entries in each row of C
    int vector_length     = 1;
    int vector_length_max = KokkosKernels::Impl::kk_get_max_vector_size<execution_space>();
    while (vector_length * 2 <= vector_length_max && (size_type)vector_length * 2 <= pot_est_nnz) {
      vector_length *= 2;
    }
    SortedCountEntriesTeam<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, alno_nnz_view_t_,
                           blno_nnz_view_t_, clno_row_view_t_, execution_space>
        countEntries(nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap);
    countEntries.sharedPerThread = pot_est_nnz;
    // compute largest possible team size
    TeamPol testPolicy(exec, 1, 1, vector_length);
    testPolicy.set_scratch_size(0, Kokkos::PerThread(pot_est_nnz * sizeof(ordinal_type)));
    int team_size = testPolicy.team_size_recommended(countEntries, Kokkos::ParallelForTag());
    // construct real policy
    int league_size = (nrows + team_size - 1) / team_size;
    TeamPol policy(exec, league_size, team_size, vector_length);
    policy.set_scratch_size(0, Kokkos::PerThread(pot_est_nnz * sizeof(ordinal_type)));
    countEntries.totalShared = countEntries.sharedPerThread * team_size * sizeof(ordinal_type);
    Kokkos::parallel_for("KokkosSparse::SpAdd::Symbolic::InputSorted::CountEntries", policy, countEntries);
  } else {
    SortedCountEntriesRange<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, alno_nnz_view_t_,
                            blno_nnz_view_t_, clno_row_view_t_, execution_space>
        countEntries(nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap);
    Kokkos::parallel_for("KokkosSparse::SpAdd::Symbolic::InputSorted::CountEntries", RangePol(exec, 0, nrows),
                         countEntries);
  }
}

// Symbolic: count entries in each row in C to produce rowmap
// kernel handle has information about whether it is sorted add or not.
template <typename execution_space, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename clno_row_view_t_>
void spadd_symbolic_impl(const execution_space& exec, KernelHandle* handle, const alno_row_view_t_ a_rowmap,
                         const alno_nnz_view_t_ a_entries, const blno_row_view_t_ b_rowmap,
                         const blno_nnz_view_t_ b_entries,
                         clno_row_view_t_ c_rowmap)  // c_rowmap must already be allocated (doesn't
                                                     // need to be initialized)
{
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;
  typedef typename KernelHandle::SPADDHandleType::nnz_lno_view_t ordinal_view_t;
  typedef typename KernelHandle::SPADDHandleType::nnz_row_view_t offset_view_t;
  // Check that A/B/C data types match KernelHandle types, and that C data types
  // are nonconst (doesn't matter if A/B types are const)
  static_assert(spadd_symbolic_same_type<typename alno_row_view_t_::non_const_value_type, size_type>,
                "add_symbolic: A size_type must match KernelHandle size_type (const "
                "doesn't matter)");
  static_assert(spadd_symbolic_same_type<typename blno_row_view_t_::non_const_value_type, size_type>,
                "add_symbolic: B size_type must match KernelHandle size_type (const "
                "doesn't matter)");
  static_assert(spadd_symbolic_same_type<typename clno_row_view_t_::non_const_value_type, size_type>,
                "add_symbolic: C size_type must match KernelHandle size_type)");
  static_assert(std::is_same_v<typename clno_row_view_t_::non_const_value_type, typename clno_row_view_t_::value_type>,
                "add_symbolic: C size_type must not be const");
  static_assert(spadd_symbolic_same_type<typename alno_nnz_view_t_::non_const_value_type, ordinal_type>,
                "add_symbolic: A entry type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");
  static_assert(spadd_symbolic_same_type<typename blno_nnz_view_t_::non_const_value_type, ordinal_type>,
                "add_symbolic: B entry type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");
  static_assert(std::is_same_v<typename clno_row_view_t_::non_const_value_type, typename clno_row_view_t_::value_type>,
                "add_symbolic: C entry type must not be const");
  // symbolic just needs to compute c_rowmap
  // easy for sorted, but for unsorted is easiest to just compute the whole sum
  auto addHandle = handle->get_spadd_handle();
  if (a_rowmap.extent(0) == 0 || a_rowmap.extent(0) == 1) {
    // Have 0 rows, so nothing to do except set #nnz to 0
    addHandle->set_c_nnz(0);
    // If c_rowmap has a single entry, it must be 0
    if (c_rowmap.extent(0)) Kokkos::deep_copy(c_rowmap, (size_type)0);
    addHandle->set_call_symbolic();
    return;
  }
  ordinal_type nrows = a_rowmap.extent(0) - 1;
  typedef Kokkos::RangePolicy<execution_space, ordinal_type> range_type;
  if (addHandle->is_input_sorted()) {
    runSortedCountEntries<execution_space, KernelHandle, alno_row_view_t_, alno_nnz_view_t_, blno_row_view_t_,
                          blno_nnz_view_t_, clno_row_view_t_>(exec, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap);
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<execution_space>(exec, nrows + 1, c_rowmap);
  } else {
    // note: scoping individual parts of the process to free views sooner,
    // minimizing peak memory usage run the unsorted c_rowmap upper bound
    // functor (just adds together A and B entry counts row by row)
    offset_view_t c_rowmap_upperbound(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "C row counts upper bound"),
                                      nrows + 1);
    size_type c_nnz_upperbound = 0;
    {
      UnsortedEntriesUpperBound<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, offset_view_t>
          countEntries(nrows, a_rowmap, b_rowmap, c_rowmap_upperbound);
      Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::CountEntries", range_type(exec, 0, nrows),
                           countEntries);
      KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<execution_space>(exec, nrows + 1, c_rowmap_upperbound,
                                                                             c_nnz_upperbound);
    }
    ordinal_view_t c_entries_uncompressed(
        Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "C entries uncompressed"), c_nnz_upperbound);
    ordinal_view_t ab_perm(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "A and B permuted entry indices"),
                           c_nnz_upperbound);
    // compute the unmerged sum
    UnmergedSumFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, offset_view_t, alno_nnz_view_t_,
                       blno_nnz_view_t_, ordinal_view_t>
        unmergedSum(nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap_upperbound, c_entries_uncompressed,
                    ab_perm);
    Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::UnmergedSum", range_type(exec, 0, nrows),
                         unmergedSum);
    // sort the unmerged sum
    KokkosSparse::sort_crs_matrix<execution_space, offset_view_t, ordinal_view_t, ordinal_view_t>(
        exec, c_rowmap_upperbound, c_entries_uncompressed, ab_perm);
    ordinal_view_t a_pos(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "A entry positions"),
                         a_entries.extent(0));
    ordinal_view_t b_pos(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "B entry positions"),
                         b_entries.extent(0));
    // merge the entries and compute Apos/Bpos, as well as Crowcounts
    {
      MergeEntriesFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, offset_view_t, ordinal_view_t,
                          offset_view_t>
          mergeEntries(nrows, a_rowmap, b_rowmap, c_rowmap_upperbound, c_rowmap, c_entries_uncompressed, ab_perm, a_pos,
                       b_pos);
      Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::MergeEntries", range_type(exec, 0, nrows),
                           mergeEntries);
      // compute actual c_rowmap
      KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<execution_space>(exec, nrows + 1, c_rowmap);
    }
    addHandle->set_a_b_pos(a_pos, b_pos);
  }
  // provide the number of NNZ in C to user through handle
  size_type cmax;
  Kokkos::deep_copy(exec, cmax, Kokkos::subview(c_rowmap, nrows));
  exec.fence("fence before cmax used on host");
  addHandle->set_c_nnz(cmax);
  addHandle->set_call_symbolic();
  addHandle->set_call_numeric(false);
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
