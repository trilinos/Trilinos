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

#ifndef _KOKKOSGRAPH_DISTANCE2_MIS_IMPL_HPP
#define _KOKKOSGRAPH_DISTANCE2_MIS_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Bitset.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_Utils.hpp"
#include <cstdint>

namespace KokkosGraph {
namespace Impl {

template <typename device_t, typename rowmap_t, typename entries_t, typename lno_view_t>
struct D2_MIS_RandomPriority {
  using exec_space     = typename device_t::execution_space;
  using mem_space      = typename device_t::memory_space;
  using bitset_t       = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type      = typename rowmap_t::non_const_value_type;
  using lno_t          = typename entries_t::non_const_value_type;
  // The type of status/priority values.
  using status_t        = typename std::make_unsigned<lno_t>::type;
  using status_view_t   = Kokkos::View<status_t*, mem_space>;
  using range_pol       = Kokkos::RangePolicy<exec_space>;
  using team_pol        = Kokkos::TeamPolicy<exec_space>;
  using team_mem        = typename team_pol::member_type;
  using all_worklists_t = Kokkos::View<lno_t**, Kokkos::LayoutLeft, mem_space>;
  using worklist_t      = Kokkos::View<lno_t*, Kokkos::LayoutLeft, mem_space>;

  // Priority values 0 and max are special, they mean the vertex is
  // in the independent set or eliminated from consideration, respectively.
  // Values in between represent a priority for being added to the set,
  // based on degree and vertex ID as a tiebreak
  //   (higher priority = less preferred to being in the independent set)

  static constexpr status_t IN_SET  = 0;
  static constexpr status_t OUT_SET = ~IN_SET;

  D2_MIS_RandomPriority(const rowmap_t& rowmap_, const entries_t& entries_)
      : rowmap(rowmap_), entries(entries_), numVerts(rowmap.extent(0) - 1) {
    status_t i = numVerts + 1;
    nvBits     = 0;
    while (i) {
      i >>= 1;
      nvBits++;
    }
    // Each value in rowStatus represents the status and priority of each row.
    // Each value in colStatus represents the lowest nonzero priority of any row
    // adjacent to the column.
    //  This counts up monotonically as vertices are eliminated (given status
    //  OUT_SET)
    rowStatus    = status_view_t(Kokkos::ViewAllocateWithoutInitializing("RowStatus"), numVerts);
    colStatus    = status_view_t(Kokkos::ViewAllocateWithoutInitializing("ColStatus"), numVerts);
    allWorklists = Kokkos::View<lno_t**, Kokkos::LayoutLeft, mem_space>(
        Kokkos::ViewAllocateWithoutInitializing("AllWorklists"), numVerts, 3);
  }

  struct RefreshRowStatus {
    RefreshRowStatus(const status_view_t& rowStatus_, const worklist_t& worklist_, lno_t nvBits_, int round)
        : rowStatus(rowStatus_), worklist(worklist_), nvBits(nvBits_) {
      hashedRound = KokkosKernels::Impl::xorshiftHash<status_t>(round);
    }

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const {
      lno_t i = worklist(w);
      // Combine vertex and round to get some pseudorandom priority bits that
      // change each round
      status_t priority =
          KokkosKernels::Impl::xorshiftHash<status_t>(KokkosKernels::Impl::xorshiftHash<status_t>(i) ^ hashedRound);
      // Generate unique status per row, with IN_SET < status < OUT_SET,
      status_t newStatus = (status_t)(i + 1) | (priority << nvBits);
      if (newStatus == OUT_SET) newStatus--;
      rowStatus(i) = newStatus;
    }

    status_view_t rowStatus;
    worklist_t worklist;
    int nvBits;
    uint32_t hashedRound;
  };

  struct RefreshColStatus {
    RefreshColStatus(const status_view_t& colStatus_, const worklist_t& worklist_, const status_view_t& rowStatus_,
                     const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_, lno_t worklistLen_)
        : colStatus(colStatus_),
          worklist(worklist_),
          rowStatus(rowStatus_),
          rowmap(rowmap_),
          entries(entries_),
          nv(nv_),
          worklistLen(worklistLen_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const {
      lno_t i = worklist(w);
      // iterate over {i} union the neighbors of i, to find
      // minimum status.
      status_t s         = rowStatus(i);
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      for (size_type j = rowBegin; j < rowEnd; j++) {
        lno_t nei = entries(j);
        if (nei < nv && nei != i) {
          status_t neiStat = rowStatus(nei);
          if (neiStat < s) s = neiStat;
        }
      }
      if (s == IN_SET) s = OUT_SET;
      colStatus(i) = s;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const team_mem& t) const {
      using MinReducer = Kokkos::Min<status_t>;
      lno_t w          = t.league_rank() * t.team_size() + t.team_rank();
      if (w >= worklistLen) return;
      lno_t i            = worklist(w);
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      lno_t rowLen       = rowEnd - rowBegin;
      // iterate over {i} union the neighbors of i, to find
      // minimum status.
      status_t s;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(t, rowLen + 1),
          [&](lno_t j, status_t& ls) {
            lno_t nei = (j == rowLen) ? i : entries(rowBegin + j);
            if (nei < nv) {
              status_t neiStat = rowStatus(nei);
              if (neiStat < ls) ls = neiStat;
            }
          },
          MinReducer(s));
      Kokkos::single(Kokkos::PerThread(t), [&]() {
        if (s == IN_SET) s = OUT_SET;
        colStatus(i) = s;
      });
    }

    status_view_t colStatus;
    worklist_t worklist;
    status_view_t rowStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
    lno_t worklistLen;
  };

  struct DecideSetFunctor {
    DecideSetFunctor(const status_view_t& rowStatus_, const status_view_t& colStatus_, const rowmap_t& rowmap_,
                     const entries_t& entries_, lno_t nv_, const worklist_t& worklist_, lno_t worklistLen_)
        : rowStatus(rowStatus_),
          colStatus(colStatus_),
          rowmap(rowmap_),
          entries(entries_),
          nv(nv_),
          worklist(worklist_),
          worklistLen(worklistLen_) {}

    // Enum values to be used as flags, so that the team policy version can
    // express the neighbor checking as an OR-reduction
    enum { NEI_OUT_SET = 1, NEI_DIFFERENT_STATUS = 2 };

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const {
      lno_t i = worklist(w);
      // Processing row i.
      status_t s = rowStatus(i);
      if (s == IN_SET || s == OUT_SET) return;
      // s is the status which must be the minimum among all neighbors
      // to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      bool neiOut        = false;
      bool neiMismatchS  = false;
      for (size_type j = rowBegin; j <= rowEnd; j++) {
        lno_t nei = (j == rowEnd) ? i : entries(j);
        if (nei >= nv) continue;
        status_t neiStat = colStatus(nei);
        if (neiStat == OUT_SET) {
          neiOut = true;
          break;
        } else if (neiStat != s) {
          neiMismatchS = true;
        }
      }
      if (neiOut) {
        // In order to make future progress, need to update the
        // col statuses for all neighbors of i.
        rowStatus(i) = OUT_SET;
      } else if (!neiMismatchS) {
        // all neighboring col statuses match s, therefore s is the minimum
        // status among all d2 neighbors
        rowStatus(i) = IN_SET;
      }
    }

    KOKKOS_INLINE_FUNCTION void operator()(const team_mem& t) const {
      using OrReducer = Kokkos::BOr<int>;
      lno_t w         = t.league_rank() * t.team_size() + t.team_rank();
      if (w >= worklistLen) return;
      lno_t i = worklist(w);
      // Processing row i.
      status_t s = rowStatus(i);
      if (s == IN_SET || s == OUT_SET) return;
      // s is the status which must be the minimum among all neighbors
      // to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      lno_t rowLen       = rowEnd - rowBegin;
      int flags          = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(t, rowLen + 1),
          [&](lno_t j, int& lflags) {
            lno_t nei = (j == rowLen) ? i : entries(rowBegin + j);
            if (nei >= nv) return;
            status_t neiStat = colStatus(nei);
            if (neiStat == OUT_SET)
              lflags |= NEI_OUT_SET;
            else if (neiStat != s)
              lflags |= NEI_DIFFERENT_STATUS;
          },
          OrReducer(flags));
      Kokkos::single(Kokkos::PerThread(t), [&]() {
        if (flags & NEI_OUT_SET) {
          // In order to make future progress, need to update the
          // col statuses for all neighbors of i.
          rowStatus(i) = OUT_SET;
        } else if (!(flags & NEI_DIFFERENT_STATUS)) {
          // all neighboring col statuses match s, therefore s is the minimum
          // status among all d2 neighbors
          rowStatus(i) = IN_SET;
        }
      });
    }

    status_view_t rowStatus;
    status_view_t colStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
    worklist_t worklist;
    lno_t worklistLen;
  };

  struct CountInSet {
    CountInSet(const status_view_t& rowStatus_) : rowStatus(rowStatus_) {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet) const {
      if (rowStatus(i) == IN_SET) lNumInSet++;
    }
    status_view_t rowStatus;
  };

  struct CompactInSet {
    CompactInSet(const status_view_t& rowStatus_, const lno_view_t& setList_)
        : rowStatus(rowStatus_), setList(setList_) {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet, bool finalPass) const {
      if (rowStatus(i) == IN_SET) {
        if (finalPass) setList(lNumInSet) = i;
        lNumInSet++;
      }
    }
    status_view_t rowStatus;
    lno_view_t setList;
  };

  struct MaskedWorklist {
    MaskedWorklist(const lno_view_t& mask_, const worklist_t& worklist_) : mask(mask_), worklist(worklist_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInList, bool finalPass) const {
      if (mask(i) < 0) {
        if (finalPass) worklist(lNumInList) = i;
        lNumInList++;
      }
    }
    lno_view_t mask;
    worklist_t worklist;
  };

  struct CompactWorklistFunctor {
    CompactWorklistFunctor(const worklist_t& src_, const worklist_t& dst_, const status_view_t& status_)
        : src(src_), dst(dst_), status(status_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w, lno_t& lNumInSet, bool finalPass) const {
      lno_t i    = src(w);
      status_t s = status(i);
      if (s != IN_SET && s != OUT_SET) {
        // next worklist needs to contain i
        if (finalPass) dst(lNumInSet) = i;
        lNumInSet++;
      }
    }

    worklist_t src;
    worklist_t dst;
    status_view_t status;
  };

  lno_view_t compute() {
    // Initialize first worklist to 0...numVerts
    worklist_t rowWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 0);
    worklist_t colWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 1);
    KokkosKernels::Impl::sequential_fill(rowWorklist);
    KokkosKernels::Impl::sequential_fill(colWorklist);
    worklist_t thirdWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 2);
    auto execSpaceEnum       = KokkosKernels::Impl::kk_get_exec_space_type<exec_space>();
    bool useTeams    = KokkosKernels::Impl::kk_is_gpu_exec_space<exec_space>() && (entries.extent(0) / numVerts >= 16);
    int vectorLength = KokkosKernels::Impl::kk_get_suggested_vector_size(numVerts, entries.extent(0), execSpaceEnum);
    int round        = 0;
    lno_t rowWorkLen = numVerts;
    lno_t colWorkLen = numVerts;
    int refreshColTeamSize = 0;
    int decideSetTeamSize  = 0;
    if (useTeams) {
      team_pol dummyPolicy(1, 1, vectorLength);
      // Compute the recommended team size for RefreshColStatus and
      // DecideSetFunctor (will be constant)
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        refreshColTeamSize = dummyPolicy.team_size_max(refreshCol, Kokkos::ParallelForTag());
      }
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        decideSetTeamSize = dummyPolicy.team_size_max(decideSet, Kokkos::ParallelForTag());
      }
    }
    while (true) {
      // Compute new row statuses
      Kokkos::parallel_for(range_pol(0, rowWorkLen), RefreshRowStatus(rowStatus, rowWorklist, nvBits, round));
      // Compute new col statuses
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        if (useTeams)
          Kokkos::parallel_for(
              team_pol((colWorkLen + refreshColTeamSize - 1) / refreshColTeamSize, refreshColTeamSize, vectorLength),
              refreshCol);
        else
          Kokkos::parallel_for(range_pol(0, colWorkLen), refreshCol);
      }
      // Decide row statuses where enough information is available
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        if (useTeams)
          Kokkos::parallel_for(
              team_pol((rowWorkLen + decideSetTeamSize - 1) / decideSetTeamSize, decideSetTeamSize, vectorLength),
              decideSet);
        else
          Kokkos::parallel_for(range_pol(0, rowWorkLen), decideSet);
      }
      round++;
      // Compact row worklist
      Kokkos::parallel_scan(range_pol(0, rowWorkLen), CompactWorklistFunctor(rowWorklist, thirdWorklist, rowStatus),
                            rowWorkLen);
      if (rowWorkLen == 0) break;
      std::swap(rowWorklist, thirdWorklist);
      // Compact col worklist
      Kokkos::parallel_scan(range_pol(0, colWorkLen), CompactWorklistFunctor(colWorklist, thirdWorklist, colStatus),
                            colWorkLen);
      std::swap(colWorklist, thirdWorklist);
    }
    // now that every vertex has been decided IN_SET/OUT_SET,
    // build a compact list of the vertices which are IN_SET.
    lno_t numInSet = 0;
    Kokkos::parallel_reduce(range_pol(0, numVerts), CountInSet(rowStatus), numInSet);
    lno_view_t setList(Kokkos::ViewAllocateWithoutInitializing("D2MIS"), numInSet);
    Kokkos::parallel_scan(range_pol(0, numVerts), CompactInSet(rowStatus, setList));
    return setList;
  }

  // Compute with an initial mask: vertices with mask value < 0 are completely
  // ignored
  lno_view_t compute(const lno_view_t& mask) {
    // Initialize first worklist to 0...numVerts
    worklist_t rowWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 0);
    worklist_t colWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 1);
    lno_t rowWorkLen       = numVerts;
    lno_t colWorkLen       = numVerts;
    // Row worklist: initially only the non-masked vertices
    Kokkos::parallel_scan(range_pol(0, numVerts), MaskedWorklist(mask, rowWorklist), rowWorkLen);
    KokkosKernels::Impl::sequential_fill(colWorklist);
    // Need to fill rowStatus with OUT_SET initially so that vertices not in the
    // worklist don't affect algorithm
    Kokkos::deep_copy(rowStatus, ~(status_t(0)));
    worklist_t thirdWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 2);
    auto execSpaceEnum       = KokkosKernels::Impl::kk_get_exec_space_type<exec_space>();
    bool useTeams    = KokkosKernels::Impl::kk_is_gpu_exec_space<exec_space>() && (entries.extent(0) / numVerts >= 16);
    int vectorLength = KokkosKernels::Impl::kk_get_suggested_vector_size(numVerts, entries.extent(0), execSpaceEnum);
    int round        = 0;
    int refreshColTeamSize = 0;
    int decideSetTeamSize  = 0;
    if (useTeams) {
      team_pol dummyPolicy(1, 1, vectorLength);
      // Compute the recommended team size for RefreshColStatus and
      // DecideSetFunctor (will be constant)
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        refreshColTeamSize = dummyPolicy.team_size_max(refreshCol, Kokkos::ParallelForTag());
      }
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        decideSetTeamSize = dummyPolicy.team_size_max(decideSet, Kokkos::ParallelForTag());
      }
    }
    while (true) {
      // Compute new row statuses
      Kokkos::parallel_for(range_pol(0, rowWorkLen), RefreshRowStatus(rowStatus, rowWorklist, nvBits, round));
      // Compute new col statuses
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        if (useTeams)
          Kokkos::parallel_for(
              team_pol((colWorkLen + refreshColTeamSize - 1) / refreshColTeamSize, refreshColTeamSize, vectorLength),
              refreshCol);
        else
          Kokkos::parallel_for(range_pol(0, colWorkLen), refreshCol);
      }
      // Decide row statuses where enough information is available
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        if (useTeams)
          Kokkos::parallel_for(
              team_pol((rowWorkLen + decideSetTeamSize - 1) / decideSetTeamSize, decideSetTeamSize, vectorLength),
              decideSet);
        else
          Kokkos::parallel_for(range_pol(0, rowWorkLen), decideSet);
      }
      round++;
      // Compact row worklist
      Kokkos::parallel_scan(range_pol(0, rowWorkLen), CompactWorklistFunctor(rowWorklist, thirdWorklist, rowStatus),
                            rowWorkLen);
      if (rowWorkLen == 0) break;
      std::swap(rowWorklist, thirdWorklist);
      // Compact col worklist
      Kokkos::parallel_scan(range_pol(0, colWorkLen), CompactWorklistFunctor(colWorklist, thirdWorklist, colStatus),
                            colWorkLen);
      std::swap(colWorklist, thirdWorklist);
    }
    // now that every vertex has been decided IN_SET/OUT_SET,
    // build a compact list of the vertices which are IN_SET.
    lno_t numInSet = 0;
    Kokkos::parallel_reduce(range_pol(0, numVerts), CountInSet(rowStatus), numInSet);
    lno_view_t setList(Kokkos::ViewAllocateWithoutInitializing("D2MIS"), numInSet);
    Kokkos::parallel_scan(range_pol(0, numVerts), CompactInSet(rowStatus, setList));
    return setList;
  }

  rowmap_t rowmap;
  entries_t entries;
  lno_t numVerts;
  status_view_t rowStatus;
  status_view_t colStatus;
  all_worklists_t allWorklists;
  // The number of bits required to represent vertex IDs, in the ECL-MIS
  // tiebreak scheme:
  //  ceil(log_2(numVerts + 1))
  int nvBits;
};

template <typename device_t, typename rowmap_t, typename entries_t, typename lno_view_t>
struct D2_MIS_FixedPriority {
  using exec_space     = typename device_t::execution_space;
  using mem_space      = typename device_t::memory_space;
  using bitset_t       = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type      = typename rowmap_t::non_const_value_type;
  using lno_t          = typename entries_t::non_const_value_type;
  // The type of status/priority values.
  using status_t      = typename std::make_unsigned<lno_t>::type;
  using status_view_t = Kokkos::View<status_t*, mem_space>;
  using range_pol     = Kokkos::RangePolicy<exec_space>;

  // Priority values 0 and max are special, they mean the vertex is
  // in the independent set or eliminated from consideration, respectively.
  // Values in between represent a priority for being added to the set,
  // based on degree and vertex ID as a tiebreak
  //   (higher priority = less preferred to being in the independent set)

  static constexpr status_t IN_SET  = 0;
  static constexpr status_t OUT_SET = ~IN_SET;

  D2_MIS_FixedPriority(const rowmap_t& rowmap_, const entries_t& entries_)
      : rowmap(rowmap_),
        entries(entries_),
        numVerts(rowmap.extent(0) - 1),
        colUpdateBitset(numVerts),
        worklist1(Kokkos::view_alloc(Kokkos::WithoutInitializing, "WL1"), numVerts),
        worklist2(Kokkos::view_alloc(Kokkos::WithoutInitializing, "WL2"), numVerts) {
    status_t i = numVerts + 1;
    nvBits     = 0;
    while (i) {
      i >>= 1;
      nvBits++;
    }
    // Each value in rowStatus represents the status and priority of each row.
    // Each value in colStatus represents the lowest nonzero priority of any row
    // adjacent to the column.
    //  This counts up monotonically as vertices are eliminated (given status
    //  OUT_SET)
    rowStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "RowStatus"), numVerts);
    colStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ColStatus"), numVerts);
    KokkosSparse::Impl::graph_min_max_degree<device_t, lno_t, rowmap_t>(rowmap, minDegree, maxDegree);
    // Compute row statuses
    Kokkos::parallel_for(range_pol(0, numVerts),
                         InitRowStatus(rowStatus, rowmap, numVerts, nvBits, minDegree, maxDegree));
    // Compute col statuses
    Kokkos::parallel_for(range_pol(0, numVerts), InitColStatus(colStatus, rowStatus, rowmap, entries, numVerts));
  }

  struct InitRowStatus {
    InitRowStatus(const status_view_t& rowStatus_, const rowmap_t& rowmap_, lno_t nv_, lno_t nvBits_, lno_t minDeg_,
                  lno_t maxDeg_)
        : rowStatus(rowStatus_),
          rowmap(rowmap_),
          nv(nv_),
          nvBits(nvBits_),
          minDeg(minDeg_),
          maxDeg(maxDeg_),
          invDegRange(1.f / (maxDeg - minDeg)) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const {
      // Generate unique status per row, with IN_SET < status < OUT_SET,
      int degBits = sizeof(status_t) * 8 - nvBits;
      if (degBits == 0) {
        // no space to store degree information. Algorithm will still work but
        // will probably produce a lower quality MIS.
        rowStatus(i) = i + 1;
        return;
      }
      status_t maxDegRange = (((status_t)1) << degBits) - 2;
      lno_t deg            = rowmap(i + 1) - rowmap(i);
      float degScore       = (float)(deg - minDeg) * invDegRange;
      rowStatus(i)         = (status_t)(i + 1) + (((status_t)(degScore * maxDegRange)) << nvBits);
    }

    status_view_t rowStatus;
    rowmap_t rowmap;
    lno_t nv;
    int nvBits;
    lno_t minDeg;
    lno_t maxDeg;
    float invDegRange;
  };

  struct InitColStatus {
    InitColStatus(const status_view_t& colStatus_, const status_view_t& rowStatus_, const rowmap_t& rowmap_,
                  const entries_t& entries_, lno_t nv_)
        : colStatus(colStatus_), rowStatus(rowStatus_), rowmap(rowmap_), entries(entries_), nv(nv_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const {
      // iterate over {i} union the neighbors of i, to find
      // minimum status.
      status_t s         = rowStatus(i);
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      for (size_type j = rowBegin; j < rowEnd; j++) {
        lno_t nei = entries(j);
        if (nei != i && nei < nv) {
          status_t neiStat = rowStatus(nei);
          if (neiStat < s) s = neiStat;
        }
      }
      colStatus(i) = s;
    }

    status_view_t colStatus;
    status_view_t rowStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
  };

  struct IterateStatusFunctor {
    IterateStatusFunctor(const status_view_t& rowStatus_, const status_view_t& colStatus_, const rowmap_t& rowmap_,
                         const entries_t& entries_, lno_t nv_, const lno_view_t& worklist_,
                         const bitset_t& colUpdateBitset_)
        : rowStatus(rowStatus_),
          colStatus(colStatus_),
          rowmap(rowmap_),
          entries(entries_),
          nv(nv_),
          worklist(worklist_),
          colUpdateBitset(colUpdateBitset_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const {
      lno_t i = worklist(w);
      // Processing row i.
      status_t s = rowStatus(i);
      // s is the status which must be the minimum among all neighbors
      // to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd   = rowmap(i + 1);
      bool neiOut        = false;
      bool neiMismatchS  = false;
      for (size_type j = rowBegin; j <= rowEnd; j++) {
        lno_t nei = (j == rowEnd) ? i : entries(j);
        if (nei >= nv) continue;
        status_t neiStat = colStatus(nei);
        if (neiStat == OUT_SET) {
          neiOut = true;
          break;
        } else if (neiStat != s) {
          neiMismatchS = true;
        }
      }
      bool statusChanged = neiOut || !neiMismatchS;
      if (neiOut) {
        // In order to make future progress, need to update the
        // col statuses for all neighbors of i which have status s.
        // This will increase the minimum to the next smallest row,
        // so that another nearby vertex can be added to the set.
        rowStatus(i) = OUT_SET;
      } else if (!neiMismatchS) {
        rowStatus(i) = IN_SET;
      }
      if (statusChanged) {
        for (size_type j = rowBegin; j <= rowEnd; j++) {
          lno_t nei = (j == rowEnd) ? i : entries(j);
          if (nei < nv && colStatus(nei) == s) colUpdateBitset.set(nei);
        }
      }
      // else: still undecided
    }

    status_view_t rowStatus;
    status_view_t colStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
    lno_view_t worklist;
    bitset_t colUpdateBitset;
  };

  struct UpdateWorklistFunctor {
    UpdateWorklistFunctor(const status_view_t& rowStatus_, const lno_view_t& oldWorklist_,
                          const lno_view_t& newWorklist_)
        : rowStatus(rowStatus_), oldWorklist(oldWorklist_), newWorklist(newWorklist_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w, lno_t& lcount, bool finalPass) const {
      // processing row i
      lno_t i = oldWorklist(w);
      // Bit i will be set when it's decided IN_SET/OUT_SET.
      // If clear, vertex i needs to be processed still.
      status_t s = rowStatus(i);
      if (s != IN_SET && s != OUT_SET) {
        if (finalPass) newWorklist(lcount) = i;
        lcount++;
      }
    }

    status_view_t rowStatus;
    lno_view_t oldWorklist;
    lno_view_t newWorklist;
  };

  struct ColRefreshWorklist {
    ColRefreshWorklist(const bitset_t& colUpdateBitset_, const lno_view_t& refreshList_)
        : colUpdateBitset(colUpdateBitset_), refreshList(refreshList_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lindex, bool finalPass) const {
      if (colUpdateBitset.test(i)) {
        if (finalPass) {
          refreshList(lindex) = i;
          colUpdateBitset.reset(i);
        }
        lindex++;
      }
    }

    bitset_t colUpdateBitset;
    lno_view_t refreshList;
  };

  struct RefreshColStatus {
    RefreshColStatus(const lno_view_t& worklist_, const status_view_t& rowStatus_, const status_view_t& colStatus_,
                     const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_)
        : worklist(worklist_),
          rowStatus(rowStatus_),
          colStatus(colStatus_),
          rowmap(rowmap_),
          entries(entries_),
          nv(nv_) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const {
      lno_t col           = worklist(w);
      status_t minNeiStat = OUT_SET;
      size_type rowBegin  = rowmap(col);
      size_type rowEnd    = rowmap(col + 1);
      for (size_type j = rowBegin; j <= rowEnd; j++) {
        lno_t nei = (j == rowEnd) ? col : entries(j);
        if (nei >= nv) continue;
        status_t neiStat = rowStatus(nei);
        if (neiStat < minNeiStat) minNeiStat = neiStat;
      }
      if (minNeiStat == IN_SET) minNeiStat = OUT_SET;
      colStatus(col) = minNeiStat;
    }

    lno_view_t worklist;
    status_view_t rowStatus;
    status_view_t colStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
  };

  struct InitWorklistFunctor {
    InitWorklistFunctor(const lno_view_t& worklist_) : worklist(worklist_) {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const { worklist(i) = i; }
    lno_view_t worklist;
  };

  struct CountInSet {
    CountInSet(const status_view_t& rowStatus_) : rowStatus(rowStatus_) {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet) const {
      if (rowStatus(i) == IN_SET) lNumInSet++;
    }
    status_view_t rowStatus;
  };

  struct CompactInSet {
    CompactInSet(const status_view_t& rowStatus_, const lno_view_t& setList_)
        : rowStatus(rowStatus_), setList(setList_) {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet, bool finalPass) const {
      if (rowStatus(i) == IN_SET) {
        if (finalPass) setList(lNumInSet) = i;
        lNumInSet++;
      }
    }
    status_view_t rowStatus;
    lno_view_t setList;
  };

  lno_view_t compute() {
    // Initialize first worklist to 0...numVerts
    Kokkos::parallel_for(range_pol(0, numVerts), InitWorklistFunctor(worklist1));
    lno_t workRemain = numVerts;
    while (workRemain) {
      // do another iteration
      Kokkos::parallel_for(range_pol(0, workRemain), IterateStatusFunctor(rowStatus, colStatus, rowmap, entries,
                                                                          numVerts, worklist1, colUpdateBitset));
      // And refresh the column statuses using the other worklist.
      lno_t colsToRefresh;
      Kokkos::parallel_scan(range_pol(0, numVerts), ColRefreshWorklist(colUpdateBitset, worklist2), colsToRefresh);
      Kokkos::parallel_for(range_pol(0, colsToRefresh),
                           RefreshColStatus(worklist2, rowStatus, colStatus, rowmap, entries, numVerts));
      // then build the next worklist with a scan. Also get the length of the
      // next worklist.
      lno_t newWorkRemain = 0;
      Kokkos::parallel_scan(range_pol(0, workRemain), UpdateWorklistFunctor(rowStatus, worklist1, worklist2),
                            newWorkRemain);
      // Finally, flip the worklists
      std::swap(worklist1, worklist2);
      workRemain = newWorkRemain;
    }
    // now that every vertex has been decided IN_SET/OUT_SET,
    // build a compact list of the vertices which are IN_SET.
    lno_t numInSet = 0;
    Kokkos::parallel_reduce(range_pol(0, numVerts), CountInSet(rowStatus), numInSet);
    lno_view_t setList(Kokkos::view_alloc(Kokkos::WithoutInitializing, "D2MIS"), numInSet);
    Kokkos::parallel_scan(range_pol(0, numVerts), CompactInSet(rowStatus, setList));
    return setList;
  }

  rowmap_t rowmap;
  entries_t entries;
  lno_t numVerts;
  status_view_t rowStatus;
  status_view_t colStatus;
  // The number of bits required to represent vertex IDs, in the ECL-MIS
  // tiebreak scheme:
  //  ceil(log_2(numVerts + 1))
  int nvBits;
  lno_t minDegree;
  lno_t maxDegree;
  // Bitset representing columns whose status needs to be recomputed
  // These bits are cleared after each refresh.
  bitset_t colUpdateBitset;
  lno_view_t worklist1;
  lno_view_t worklist2;
};

template <typename device_t, typename rowmap_t, typename entries_t, typename labels_t>
struct D2_MIS_Aggregation {
  using exec_space     = typename device_t::execution_space;
  using mem_space      = typename device_t::memory_space;
  using bitset_t       = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type      = typename rowmap_t::non_const_value_type;
  using char_view_t    = Kokkos::View<char*, mem_space>;
  using lno_t          = typename entries_t::non_const_value_type;
  using lno_view_t     = typename entries_t::non_const_type;
  // The type of status/priority values.
  using status_t      = typename std::make_unsigned<lno_t>::type;
  using status_view_t = Kokkos::View<status_t*, mem_space>;
  using range_pol     = Kokkos::RangePolicy<exec_space>;
  using mis2_view     = Kokkos::View<lno_t*, mem_space>;

  D2_MIS_Aggregation(const rowmap_t& rowmap_, const entries_t& entries_)
      : rowmap(rowmap_),
        entries(entries_),
        numVerts(rowmap.extent(0) - 1),
        labels(Kokkos::ViewAllocateWithoutInitializing("AggregateLabels"), numVerts),
        roots("Root Status", numVerts) {
    Kokkos::deep_copy(labels, (lno_t)-1);
  }

  struct Phase1Functor {
    Phase1Functor(lno_t numVerts__, const mis2_view& m1__, const rowmap_t& rowmap__, const entries_t& entries__,
                  const labels_t& labels__, const char_view_t& roots__)
        : numVerts_(numVerts__),
          m1_(m1__),
          rowmap_(rowmap__),
          entries_(entries__),
          labels_(labels__),
          roots_(roots__) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t agg) const {
      lno_t root         = m1_(agg);
      roots_(root)       = 1;
      size_type rowBegin = rowmap_(root);
      size_type rowEnd   = rowmap_(root + 1);
      labels_(root)      = agg;
      for (size_type j = rowBegin; j < rowEnd; j++) {
        lno_t nei = entries_(j);
        if (nei < numVerts_) labels_(nei) = agg;
      }
    }

    lno_t numVerts_;
    mis2_view m1_;
    rowmap_t rowmap_;
    entries_t entries_;
    labels_t labels_;
    char_view_t roots_;
  };

  void createPrimaryAggregates() {
    // Compute an MIS-2
    D2_MIS_RandomPriority<device_t, rowmap_t, entries_t, mis2_view> d2mis(rowmap, entries);
    mis2_view m1 = d2mis.compute();
    // Construct initial aggregates using roots and all direct neighbors
    Kokkos::parallel_for(range_pol(0, m1.extent(0)), Phase1Functor(numVerts, m1, rowmap, entries, labels, roots));
    numAggs = m1.extent(0);
  }

  struct CandAggSizesFunctor {
    CandAggSizesFunctor(lno_t numVerts__, const labels_t& m2__, const rowmap_t& rowmap__, const entries_t& entries__,
                        const labels_t& labels__, const labels_t& candAggSizes__)
        : numVerts_(numVerts__),
          m2_(m2__),
          rowmap_(rowmap__),
          entries_(entries__),
          labels_(labels__),
          candAggSizes_(candAggSizes__) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const {
      lno_t candRoot = m2_(i);
      // Count the number of non-aggregated neighbors, including self
      lno_t aggSize      = 1;
      size_type rowBegin = rowmap_(candRoot);
      size_type rowEnd   = rowmap_(candRoot + 1);
      for (size_type j = rowBegin; j < rowEnd; j++) {
        lno_t nei = entries_(j);
        if (nei == candRoot || nei >= numVerts_) continue;
        if (labels_(nei) == -1) aggSize++;
      }
      candAggSizes_(i) = aggSize;
    }

    lno_t numVerts_;
    labels_t m2_;
    rowmap_t rowmap_;
    entries_t entries_;
    labels_t labels_;
    labels_t candAggSizes_;
  };

  struct ChoosePhase2AggsFunctor {
    ChoosePhase2AggsFunctor(lno_t numVerts__, lno_t numAggs__, const labels_t& m2__, const rowmap_t& rowmap__,
                            const entries_t& entries__, const labels_t& labels__, const labels_t& candAggSizes__,
                            const char_view_t& roots__)
        : numVerts_(numVerts__),
          numAggs_(numAggs__),
          m2_(m2__),
          rowmap_(rowmap__),
          entries_(entries__),
          labels_(labels__),
          candAggSizes_(candAggSizes__),
          roots_(roots__) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lid, bool finalPass) const {
      lno_t aggSize = candAggSizes_(i);
      if (aggSize < 3) return;
      if (finalPass) {
        // Build the aggregate
        lno_t root         = m2_(i);
        lno_t aggID        = numAggs_ + lid;
        labels_(root)      = aggID;
        roots_(root)       = 1;
        size_type rowBegin = rowmap_(root);
        size_type rowEnd   = rowmap_(root + 1);
        for (size_type j = rowBegin; j < rowEnd; j++) {
          lno_t nei = entries_(j);
          if (nei == root || nei >= numVerts_) continue;
          if (labels_(nei) == -1) labels_(nei) = aggID;
        }
      }
      lid++;
    }

    lno_t numVerts_;
    lno_t numAggs_;
    labels_t m2_;
    rowmap_t rowmap_;
    entries_t entries_;
    labels_t labels_;
    labels_t candAggSizes_;
    char_view_t roots_;
  };

  void createSecondaryAggregates() {
    labels_t candAggSizes(Kokkos::ViewAllocateWithoutInitializing("Phase2 Candidate Agg Sizes"), numVerts);
    // Compute a new MIS-2 from only unaggregated nodes
    D2_MIS_RandomPriority<device_t, rowmap_t, entries_t, labels_t> d2mis(rowmap, entries);
    labels_t m2        = d2mis.compute(labels);
    lno_t numCandRoots = m2.extent(0);
    // Compute the sizes of would-be aggregates.
    Kokkos::parallel_for(range_pol(0, numCandRoots),
                         CandAggSizesFunctor(numVerts, m2, rowmap, entries, labels, candAggSizes));
    // Now, filter out the candidate aggs which are big enough, and create those
    // aggregates. Using a scan for this assigns IDs deterministically (unlike
    // an atomic counter).
    lno_t numNewAggs = 0;
    Kokkos::parallel_scan(range_pol(0, numCandRoots),
                          ChoosePhase2AggsFunctor(numVerts, numAggs, m2, rowmap, entries, labels, candAggSizes, roots),
                          numNewAggs);
    numAggs += numNewAggs;
  }

  struct SizeAndConnectivityFunctor {
    SizeAndConnectivityFunctor(lno_t numVerts__, const rowmap_t& rowmap__, const entries_t& entries__,
                               const labels_t& labels__, const labels_t& connectivities__, const labels_t& aggSizes__)
        : numVerts_(numVerts__),
          rowmap_(rowmap__),
          entries_(entries__),
          labels_(labels__),
          connectivities_(connectivities__),
          aggSizes_(aggSizes__) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const {
      lno_t agg = labels_(i);
      if (agg != -1) {
        Kokkos::atomic_increment(&aggSizes_(agg));
        // compute connectivity of i
        size_type rowBegin = rowmap_(i);
        size_type rowEnd   = rowmap_(i + 1);
        int connect        = 0;
        for (size_type j = rowBegin; j < rowEnd; j++) {
          lno_t nei = entries_(j);
          if (nei == i || nei >= numVerts_) continue;
          lno_t neiAgg = labels_(nei);
          if (neiAgg == agg) connect++;
        }
        connectivities_(i) = connect;
      } else
        connectivities_(i) = 0;
    }

    lno_t numVerts_;
    rowmap_t rowmap_;
    entries_t entries_;
    labels_t labels_;
    labels_t connectivities_;
    labels_t aggSizes_;
  };

  struct AssignLeftoverFunctor {
    AssignLeftoverFunctor(lno_t numVerts__, const rowmap_t& rowmap__, const entries_t& entries__,
                          const labels_t& labels__, const labels_t& labelsOld__, const labels_t& connectivities__,
                          const labels_t& aggSizes__, const char_view_t& roots__)
        : numVerts_(numVerts__),
          rowmap_(rowmap__),
          entries_(entries__),
          labels_(labels__),
          labelsOld_(labelsOld__),
          connectivities_(connectivities__),
          aggSizes_(aggSizes__),
          roots_(roots__) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const {
      lno_t agg               = labelsOld_(i);
      lno_t trackedAggs[8]    = {0};
      lno_t trackedConnect[8] = {0};
      char trackedRootAdj[8]  = {0};
      if (roots_(i)) {
        // Never reassign roots
        labels_(i) = agg;
        return;
      }
      lno_t numTracked   = 0;
      size_type rowBegin = rowmap_(i);
      size_type rowEnd   = rowmap_(i + 1);
      for (size_type j = rowBegin; j < rowEnd; j++) {
        lno_t nei = entries_(j);
        if (nei == i || nei >= numVerts_) continue;
        lno_t neiAgg = labelsOld_(nei);
        if (neiAgg == -1 || neiAgg == agg) continue;
        // Try to find neiAgg in tracked
        lno_t trackedID = -1;
        for (lno_t k = 0; k < numTracked; k++) {
          if (trackedAggs[k] == neiAgg) {
            trackedID = k;
            break;
          }
        }
        if (trackedID == -1) {
          if (numTracked < 8) {
            trackedID              = numTracked++;
            trackedAggs[trackedID] = neiAgg;
          } else {
            // Ran out of space, just ignore this neighboring agg
            break;
          }
        }
        // Record the connectivity to this neighboring agg
        if (roots_(nei)) trackedRootAdj[trackedID] = 1;
        trackedConnect[trackedID]++;
      }
      // Now that we know connectivity of this node to (hopefully all)
      // neighboring aggs, decide if it's better to join that agg instead
      char bestRootAdj  = (agg >= 0) ? 1 : 0;
      lno_t bestAgg     = agg;
      lno_t bestConnect = connectivities_(i);
      // If not in an agg, initial bestSize doesn't matter because any
      // neighboring agg has a better connectivity than 0
      lno_t bestSize = (agg >= 0) ? aggSizes_(agg) : 0;
      for (int k = 0; k < numTracked; k++) {
        lno_t s = aggSizes_(trackedAggs[k]);
        // Priorities: adjacent to root > connect > size
        if (trackedRootAdj[k] > bestRootAdj ||
            (trackedRootAdj[k] == bestRootAdj &&
             ((trackedConnect[k] > bestConnect) || (trackedConnect[k] == bestConnect && s < bestSize)))) {
          bestRootAdj = trackedRootAdj[k];
          bestConnect = trackedConnect[k];
          bestSize    = s;
          bestAgg     = trackedAggs[k];
        }
      }
      labels_(i) = bestAgg;
    }

    lno_t numVerts_;
    rowmap_t rowmap_;
    entries_t entries_;
    labels_t labels_;
    labels_t labelsOld_;
    labels_t connectivities_;
    labels_t aggSizes_;
    char_view_t roots_;
  };

  void aggregateLeftovers() {
    // Phase3 is cleanup. All aggregates have already been created, but some
    // vertices might be unaggregated. Compute the current size of each
    // aggregate, and then join each unaggregated node to the smallest
    // neighboring aggregate.
    labels_t labelsOld("old", numVerts);
    Kokkos::deep_copy(labelsOld, labels);
    labels_t connectivities(Kokkos::ViewAllocateWithoutInitializing("connect"), numVerts);
    labels_t aggSizes("Phase3 Agg Sizes", numAggs);
    Kokkos::parallel_for(range_pol(0, numVerts),
                         SizeAndConnectivityFunctor(numVerts, rowmap, entries, labels, connectivities, aggSizes));
    // Now, join vertices to aggregates
    Kokkos::parallel_for(range_pol(0, numVerts), AssignLeftoverFunctor(numVerts, rowmap, entries, labels, labelsOld,
                                                                       connectivities, aggSizes, roots));
  }

  // phase 2 creates new aggregates in between the initial MIS-2 neighborhoods.
  // Effectively slows coarsening rate by adding new aggregates.
  void compute(bool enableSecondaryAggregates) {
    //  * Phase 1: compute MIS-2, construct a 'primary' aggregate from each
    //  in-set point and its neighbors
    createPrimaryAggregates();
    //  * Phase 2:
    //    - Compute an MIS-2 on subgraph of unaggregated nodes
    //    - For each in-set point:
    //      - Count unaggregated neighbors.
    //      - If total agg size would be >= 3, create the new aggregate.
    //    - This is optional: enabling this phase slows coarsening rate (i.e.
    //    coarse graph is larger)
    if (enableSecondaryAggregates) createSecondaryAggregates();
    //  * Phase 3: join still unaggregated (leftover) vertices to a neighboring
    //  aggregate
    //    - Ideally, the smallest neighboring aggregate.
    //    - To remain deterministic, we use the agg sizes from end of
    //    phase 2 and hold those constant during phase 3.
    aggregateLeftovers();
  }

  rowmap_t rowmap;
  entries_t entries;
  lno_t numVerts;
  lno_t numAggs;
  labels_t labels;
  char_view_t roots;
};

}  // namespace Impl
}  // namespace KokkosGraph

#endif
