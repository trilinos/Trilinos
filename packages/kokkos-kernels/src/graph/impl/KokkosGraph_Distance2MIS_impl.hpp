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
// Questions? Contact Brian Kelley (bmkelle@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSGRAPH_DISTANCE2_MIS_IMPL_HPP
#define _KOKKOSGRAPH_DISTANCE2_MIS_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Bitset.hpp"
#include "KokkosKernels_Utils.hpp"
#include <cstdint>

namespace KokkosGraph {
namespace Experimental {
namespace Impl {

template<typename device_t, typename rowmap_t, typename entries_t, typename lno_view_t>
struct D2_MIS_RandomPriority
{
  using exec_space = typename device_t::execution_space;
  using mem_space = typename device_t::memory_space;
  using bitset_t = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t = typename entries_t::non_const_value_type;
  //The type of status/priority values.
  using status_t = typename std::make_unsigned<lno_t>::type;
  using status_view_t = Kokkos::View<status_t*, mem_space>;
  using range_pol = Kokkos::RangePolicy<exec_space>;
  using team_pol = Kokkos::TeamPolicy<exec_space>;
  using team_mem = typename team_pol::member_type;
  using all_worklists_t = Kokkos::View<lno_t**, Kokkos::LayoutLeft, mem_space>;
  using worklist_t = Kokkos::View<lno_t*, Kokkos::LayoutLeft, mem_space>;

  KOKKOS_INLINE_FUNCTION static uint32_t xorshiftHash(uint32_t in)
  {
    uint32_t x = in;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
  }

  // Priority values 0 and max are special, they mean the vertex is
  // in the independent set or eliminated from consideration, respectively.
  // Values in between represent a priority for being added to the set,
  // based on degree and vertex ID as a tiebreak
  //   (higher priority = less preferred to being in the independent set)

  static constexpr status_t IN_SET = 0;
  static constexpr status_t OUT_SET = ~IN_SET;

  D2_MIS_RandomPriority(const rowmap_t& rowmap_, const entries_t& entries_)
    : rowmap(rowmap_), entries(entries_), numVerts(rowmap.extent(0) - 1)
  {
    status_t i = numVerts + 1;
    nvBits = 0;
    while(i)
    {
      i >>= 1;
      nvBits++;
    }
    //Each value in rowStatus represents the status and priority of each row.
    //Each value in colStatus represents the lowest nonzero priority of any row adjacent to the column.
    //  This counts up monotonically as vertices are eliminated (given status OUT_SET)
    rowStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "RowStatus"), numVerts);
    colStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ColStatus"), numVerts);
    allWorklists = Kokkos::View<lno_t**, Kokkos::LayoutLeft, mem_space>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "AllWorklists"), numVerts, 3);
  }

  struct RefreshRowStatus
  {
    RefreshRowStatus(const status_view_t& rowStatus_, const worklist_t& worklist_, lno_t nvBits_, int round)
      : rowStatus(rowStatus_), worklist(worklist_), nvBits(nvBits_)
    {
      hashedRound = xorshiftHash(round);
    }

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
    {
      lno_t i = worklist(w);
      //Combine vertex and round to get some pseudorandom priority bits that change each round
      status_t priority = xorshiftHash(i + hashedRound);
      //Generate unique status per row, with IN_SET < status < OUT_SET,
      int priorityBits = sizeof(status_t) * 8 - nvBits;
      status_t priorityMask = 1;
      priorityMask <<= priorityBits;
      priorityMask--;
      status_t newStatus = (status_t) (i + 1) + ((priority & priorityMask) << nvBits);
      if(newStatus == OUT_SET)
        newStatus--;
      rowStatus(i) = newStatus;
    }

    status_view_t rowStatus;
    worklist_t worklist;
    int nvBits;
    uint32_t hashedRound;
  };

  struct RefreshColStatus
  {
    RefreshColStatus(const status_view_t& colStatus_, const worklist_t& worklist_, const status_view_t& rowStatus_, const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_, lno_t worklistLen_)
      : colStatus(colStatus_), worklist(worklist_), rowStatus(rowStatus_), rowmap(rowmap_), entries(entries_), nv(nv_), worklistLen(worklistLen_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
    {
      lno_t i = worklist(w);
      //iterate over {i} union the neighbors of i, to find
      //minimum status.
      status_t s = OUT_SET;
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      for(size_type j = rowBegin; j <= rowEnd; j++)
      {
        lno_t nei = (j == rowEnd) ? i : entries(j);
        if(nei < nv)
        {
          status_t neiStat = rowStatus(nei);
          if(neiStat < s)
            s = neiStat;
        }
      }
      if(s == IN_SET)
        s = OUT_SET;
      colStatus(i) = s;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const team_mem& t) const
    {
      using MinReducer = Kokkos::Min<status_t>;
      lno_t w = t.league_rank() * t.team_size() + t.team_rank();
      if(w >= worklistLen)
        return;
      lno_t i = worklist(w);
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      lno_t rowLen = rowEnd - rowBegin;
      //iterate over {i} union the neighbors of i, to find
      //minimum status.
      status_t s;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(t, rowLen + 1),
      [&](lno_t j, status_t& ls)
      {
        lno_t nei = (j == rowLen) ? i : entries(rowBegin + j);
        if(nei < nv)
        {
          status_t neiStat = rowStatus(nei);
          if(neiStat < ls)
            ls = neiStat;
        }
      }, MinReducer(s));
      Kokkos::single(Kokkos::PerThread(t),
      [&]()
      {
        if(s == IN_SET)
          s = OUT_SET;
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

  struct DecideSetFunctor
  {
    DecideSetFunctor(const status_view_t& rowStatus_, const status_view_t& colStatus_, const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_, const worklist_t& worklist_, lno_t worklistLen_)
      : rowStatus(rowStatus_), colStatus(colStatus_), rowmap(rowmap_), entries(entries_), nv(nv_), worklist(worklist_), worklistLen(worklistLen_)
    {}

    //Enum values to be used as flags, so that the team policy version can
    //express the neighbor checking as an OR-reduction
    enum
    {
      NEI_OUT_SET = 1,
      NEI_DIFFERENT_STATUS = 2
    };

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
    {
      lno_t i = worklist(w);
      //Processing row i.
      status_t s = rowStatus(i);
      if(s == IN_SET || s == OUT_SET)
        return;
      //s is the status which must be the minimum among all neighbors
      //to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      bool neiOut = false;
      bool neiMismatchS = false;
      for(size_type j = rowBegin; j <= rowEnd; j++)
      {
        lno_t nei = (j == rowEnd) ? i : entries(j);
        if(nei >= nv)
          continue;
        status_t neiStat = colStatus(nei);
        if(neiStat == OUT_SET)
        {
          neiOut = true;
          break;
        }
        else if(neiStat != s)
        {
          neiMismatchS = true;
        }
      }
      if(neiOut)
      {
        //In order to make future progress, need to update the
        //col statuses for all neighbors of i.
        rowStatus(i) = OUT_SET;
      }
      else if(!neiMismatchS)
      {
        //all neighboring col statuses match s, therefore s is the minimum status among all d2 neighbors
        rowStatus(i) = IN_SET;
      }
    }

    KOKKOS_INLINE_FUNCTION void operator()(const team_mem& t) const
    {
      using OrReducer = Kokkos::BOr<int>;
      lno_t w = t.league_rank() * t.team_size() + t.team_rank();
      if(w >= worklistLen)
        return;
      lno_t i = worklist(w);
      //Processing row i.
      status_t s = rowStatus(i);
      if(s == IN_SET || s == OUT_SET)
        return;
      //s is the status which must be the minimum among all neighbors
      //to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      lno_t rowLen = rowEnd - rowBegin;
      int flags = 0;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(t, rowLen + 1),
      [&](lno_t j, int& lflags)
      {
        lno_t nei = (j == rowLen) ? i : entries(rowBegin + j);
        if(nei >= nv)
          return;
        status_t neiStat = colStatus(nei);
        if(neiStat == OUT_SET)
          lflags |= NEI_OUT_SET;
        else if(neiStat != s)
          lflags |= NEI_DIFFERENT_STATUS;
      }, OrReducer(flags));
      Kokkos::single(Kokkos::PerThread(t),
      [&]()
      {
        if(flags & NEI_OUT_SET)
        {
          //In order to make future progress, need to update the
          //col statuses for all neighbors of i.
          rowStatus(i) = OUT_SET;
        }
        else if(!(flags & NEI_DIFFERENT_STATUS))
        {
          //all neighboring col statuses match s, therefore s is the minimum status among all d2 neighbors
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

  struct CountInSet
  {
    CountInSet(const status_view_t& rowStatus_)
      : rowStatus(rowStatus_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet) const
    {
      if(rowStatus(i) == IN_SET)
        lNumInSet++;
    }
    status_view_t rowStatus;
  };

  struct CompactInSet
  {
    CompactInSet(const status_view_t& rowStatus_, const lno_view_t& setList_)
      : rowStatus(rowStatus_), setList(setList_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet, bool finalPass) const
    {
      if(rowStatus(i) == IN_SET)
      {
        if(finalPass)
          setList(lNumInSet) = i;
        lNumInSet++;
      }
    }
    status_view_t rowStatus;
    lno_view_t setList;
  };

  struct InitWorklistFunctor
  {
    InitWorklistFunctor(const worklist_t& worklist_)
      : worklist(worklist_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      worklist(i) = i;
    }
    worklist_t worklist;
  };

  struct CompactWorklistFunctor
  {
    CompactWorklistFunctor(const worklist_t& src_, const worklist_t& dst_, const status_view_t& status_)
      : src(src_), dst(dst_), status(status_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w, lno_t& lNumInSet, bool finalPass) const
    {
      lno_t i = src(w);
      status_t s = status(i);
      if(s != IN_SET && s != OUT_SET)
      {
        //next worklist needs to contain i
        if(finalPass)
          dst(lNumInSet) = i;
        lNumInSet++;
      }
    }

    worklist_t src;
    worklist_t dst;
    status_view_t status;
  };

  lno_view_t compute()
  {
    //Initialize first worklist to 0...numVerts
    worklist_t rowWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 0);
    Kokkos::parallel_for(range_pol(0, numVerts), InitWorklistFunctor(rowWorklist));
    worklist_t colWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 1);
    Kokkos::parallel_for(range_pol(0, numVerts), InitWorklistFunctor(colWorklist));
    worklist_t thirdWorklist = Kokkos::subview(allWorklists, Kokkos::ALL(), 2);
    auto execSpaceEnum = KokkosKernels::Impl::kk_get_exec_space_type<exec_space>();
    bool useTeams = KokkosKernels::Impl::kk_is_gpu_exec_space<exec_space>() && (entries.extent(0) / numVerts >= 16);
    int vectorLength = KokkosKernels::Impl::kk_get_suggested_vector_size(numVerts, entries.extent(0), execSpaceEnum);
    int round = 0;
    lno_t rowWorkLen = numVerts;
    lno_t colWorkLen = numVerts;
    int refreshColTeamSize = 0;
    int decideSetTeamSize = 0;
    if(useTeams)
    {
      team_pol dummyPolicy(1, 1, vectorLength);
      //Compute the recommended team size for RefreshColStatus and DecideSetFunctor (will be constant)
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        refreshColTeamSize = dummyPolicy.team_size_max(refreshCol, Kokkos::ParallelForTag());
      }
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        decideSetTeamSize = dummyPolicy.team_size_max(decideSet, Kokkos::ParallelForTag());
      }
    }
    while(true)
    {
      //Compute new row statuses
      Kokkos::parallel_for(range_pol(0, rowWorkLen), RefreshRowStatus(rowStatus, rowWorklist, nvBits, round));
      //Compute new col statuses
      {
        RefreshColStatus refreshCol(colStatus, colWorklist, rowStatus, rowmap, entries, numVerts, colWorkLen);
        if(useTeams)
          Kokkos::parallel_for(team_pol((colWorkLen + refreshColTeamSize - 1) / refreshColTeamSize, refreshColTeamSize, vectorLength), refreshCol);
        else
          Kokkos::parallel_for(range_pol(0, colWorkLen), refreshCol);
      }
      //Decide row statuses where enough information is available
      {
        DecideSetFunctor decideSet(rowStatus, colStatus, rowmap, entries, numVerts, rowWorklist, rowWorkLen);
        if(useTeams)
          Kokkos::parallel_for(team_pol((rowWorkLen + decideSetTeamSize - 1) / decideSetTeamSize, decideSetTeamSize, vectorLength), decideSet);
        else
          Kokkos::parallel_for(range_pol(0, rowWorkLen), decideSet);
      }
      //Compact row worklist
      Kokkos::parallel_scan(range_pol(0, rowWorkLen), CompactWorklistFunctor(rowWorklist, thirdWorklist, rowStatus), rowWorkLen);
      if(rowWorkLen == 0)
        break;
      std::swap(rowWorklist, thirdWorklist);
      //Compact col worklist
      Kokkos::parallel_scan(range_pol(0, colWorkLen), CompactWorklistFunctor(colWorklist, thirdWorklist, colStatus), colWorkLen);
      std::swap(colWorklist, thirdWorklist);
      round++;
    }
    //now that every vertex has been decided IN_SET/OUT_SET,
    //build a compact list of the vertices which are IN_SET.
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
  all_worklists_t allWorklists;
  //The number of bits required to represent vertex IDs, in the ECL-MIS tiebreak scheme:
  //  ceil(log_2(numVerts + 1))
  int nvBits;
};

//    UNUSED CODE
//    Version of RefreshRowStatus, which does linear interpolation between a degree-based score and a random score.
//    By gradually increasing the interpolation coefficient in favor of random, the MIS can converge much faster than
//    constant priorities.
//
//    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
//    {
//      lno_t i = worklist(w);
//      int degBits = sizeof(status_t) * 8 - nvBits;
//      if(degBits == 0)
//      {
//        //no space to store degree information. Algorithm will still work but will
//        //probably produce a lower quality MIS.
//        rowStatus(i) = i + 1;
//        return;
//      }
//      //Combine vertex and round to get some pseudorandom priority bits that change each round
//      status_t maxDegRange = (((status_t) 1) << degBits) - 2;
//      lno_t deg = rowmap(i + 1) - rowmap(i);
//      //Compute degree-based score and random score
//      float degScore = (float) (deg - minDeg) * invDegRange;
//      float randScore = (xorshiftHash(i + hashedRound) & 0xFFFF) / 65536.f;
//      //Then linearly interpolate using k
//      float finalScore = k * randScore + (1.f - k) * degScore;
//      rowStatus(i) = (status_t) (i + 1) + (((status_t) (finalScore * maxDegRange)) << nvBits);
//    }
//    */

template<typename device_t, typename rowmap_t, typename entries_t, typename lno_view_t>
struct D2_MIS_FixedPriority
{
  using exec_space = typename device_t::execution_space;
  using mem_space = typename device_t::memory_space;
  using bitset_t = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t = typename entries_t::non_const_value_type;
  //The type of status/priority values.
  using status_t = typename std::make_unsigned<lno_t>::type;
  using status_view_t = Kokkos::View<status_t*, mem_space>;
  using range_pol = Kokkos::RangePolicy<exec_space>;

  // Priority values 0 and max are special, they mean the vertex is
  // in the independent set or eliminated from consideration, respectively.
  // Values in between represent a priority for being added to the set,
  // based on degree and vertex ID as a tiebreak
  //   (higher priority = less preferred to being in the independent set)

  static constexpr status_t IN_SET = 0;
  static constexpr status_t OUT_SET = ~IN_SET;

  D2_MIS_FixedPriority(const rowmap_t& rowmap_, const entries_t& entries_)
    : rowmap(rowmap_), entries(entries_), numVerts(rowmap.extent(0) - 1), colUpdateBitset(numVerts),
    worklist1(Kokkos::view_alloc(Kokkos::WithoutInitializing, "WL1"), numVerts),
    worklist2(Kokkos::view_alloc(Kokkos::WithoutInitializing, "WL2"), numVerts)
  {
    status_t i = numVerts + 1;
    nvBits = 0;
    while(i)
    {
      i >>= 1;
      nvBits++;
    }
    //Each value in rowStatus represents the status and priority of each row.
    //Each value in colStatus represents the lowest nonzero priority of any row adjacent to the column.
    //  This counts up monotonically as vertices are eliminated (given status OUT_SET)
    rowStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "RowStatus"), numVerts);
    colStatus = status_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ColStatus"), numVerts);
    KokkosKernels::Impl::graph_min_max_degree<device_t, lno_t, rowmap_t>(rowmap, minDegree, maxDegree);
    //Compute row statuses 
    Kokkos::parallel_for(range_pol(0, numVerts), InitRowStatus(rowStatus, rowmap, numVerts, nvBits, minDegree, maxDegree));
    //Compute col statuses
    Kokkos::parallel_for(range_pol(0, numVerts), InitColStatus(colStatus, rowStatus, rowmap, entries, numVerts));
  }

  struct InitRowStatus
  {
    InitRowStatus(const status_view_t& rowStatus_, const rowmap_t& rowmap_, lno_t nv_, lno_t nvBits_, lno_t minDeg_, lno_t maxDeg_)
      : rowStatus(rowStatus_), rowmap(rowmap_), nv(nv_), nvBits(nvBits_), minDeg(minDeg_), maxDeg(maxDeg_), invDegRange(1.f / (maxDeg - minDeg)) {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      //Generate unique status per row, with IN_SET < status < OUT_SET,
      int degBits = sizeof(status_t) * 8 - nvBits;
      if(degBits == 0)
      {
        //no space to store degree information. Algorithm will still work but will
        //probably produce a lower quality MIS.
        rowStatus(i) = i + 1;
        return;
      }
      status_t maxDegRange = (((status_t) 1) << degBits) - 2;
      lno_t deg = rowmap(i + 1) - rowmap(i);
      float degScore = (float) (deg - minDeg) * invDegRange;
      rowStatus(i) = (status_t) (i + 1) + (((status_t) (degScore * maxDegRange)) << nvBits);
    }

    status_view_t rowStatus;
    rowmap_t rowmap;
    lno_t nv;
    int nvBits;
    lno_t minDeg;
    lno_t maxDeg;
    float invDegRange;
  };

  struct InitColStatus
  {
    InitColStatus(const status_view_t& colStatus_, const status_view_t& rowStatus_, const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_)
      : colStatus(colStatus_), rowStatus(rowStatus_), rowmap(rowmap_), entries(entries_), nv(nv_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      //iterate over {i} union the neighbors of i, to find
      //minimum status.
      status_t s = rowStatus(i);
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      for(size_type j = rowBegin; j < rowEnd; j++)
      {
        lno_t nei = entries(j);
        if(nei != i && nei < nv)
        {
          status_t neiStat = rowStatus(nei);
          if(neiStat < s)
            s = neiStat;
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

  struct IterateStatusFunctor
  {
    IterateStatusFunctor(const status_view_t& rowStatus_, const status_view_t& colStatus_, const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_, const lno_view_t& worklist_, const bitset_t& colUpdateBitset_)
      : rowStatus(rowStatus_), colStatus(colStatus_), rowmap(rowmap_), entries(entries_), nv(nv_), worklist(worklist_), colUpdateBitset(colUpdateBitset_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
    {
      lno_t i = worklist(w);
      //Processing row i.
      status_t s = rowStatus(i);
      //s is the status which must be the minimum among all neighbors
      //to decide that i is IN_SET.
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      bool neiOut = false;
      bool neiMismatchS = false;
      for(size_type j = rowBegin; j <= rowEnd; j++)
      {
        lno_t nei = (j == rowEnd) ? i : entries(j);
        if(nei >= nv)
          continue;
        status_t neiStat = colStatus(nei);
        if(neiStat == OUT_SET)
        {
          neiOut = true;
          break;
        }
        else if(neiStat != s)
        {
          neiMismatchS = true;
        }
      }
      bool statusChanged = neiOut || !neiMismatchS;
      if(neiOut)
      {
        //In order to make future progress, need to update the
        //col statuses for all neighbors of i which have status s.
        //This will increase the minimum to the next smallest row,
        //so that another nearby vertex can be added to the set.
        rowStatus(i) = OUT_SET;
      }
      else if(!neiMismatchS)
      {
        rowStatus(i) = IN_SET;
      }
      if(statusChanged)
      {
        for(size_type j = rowBegin; j <= rowEnd; j++)
        {
          lno_t nei = (j == rowEnd) ? i : entries(j);
          if(nei < nv && colStatus(nei) == s)
            colUpdateBitset.set(nei);
        }
      }
      //else: still undecided
    }

    status_view_t rowStatus;
    status_view_t colStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
    lno_view_t worklist;
    bitset_t colUpdateBitset;
  };

  struct UpdateWorklistFunctor
  {
    UpdateWorklistFunctor(const status_view_t& rowStatus_, const lno_view_t& oldWorklist_, const lno_view_t& newWorklist_)
      : rowStatus(rowStatus_), oldWorklist(oldWorklist_), newWorklist(newWorklist_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w, lno_t& lcount, bool finalPass) const
    {
      //processing row i
      lno_t i = oldWorklist(w);
      //Bit i will be set when it's decided IN_SET/OUT_SET.
      //If clear, vertex i needs to be processed still.
      status_t s = rowStatus(i);
      if(s != IN_SET && s != OUT_SET)
      {
        if(finalPass)
          newWorklist(lcount) = i;
        lcount++;
      }
    }

    status_view_t rowStatus;
    lno_view_t oldWorklist;
    lno_view_t newWorklist;
  };

  struct ColRefreshWorklist
  {
    ColRefreshWorklist(const bitset_t& colUpdateBitset_, const lno_view_t& refreshList_)
      : colUpdateBitset(colUpdateBitset_), refreshList(refreshList_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lindex, bool finalPass) const
    {
      if(colUpdateBitset.test(i))
      {
        if(finalPass)
        {
          refreshList(lindex) = i;
          colUpdateBitset.reset(i);
        }
        lindex++;
      }
    }

    bitset_t colUpdateBitset;
    lno_view_t refreshList;
  };

  struct RefreshColStatus
  {
    RefreshColStatus(const lno_view_t& worklist_, const status_view_t& rowStatus_, const status_view_t& colStatus_, const rowmap_t& rowmap_, const entries_t& entries_, lno_t nv_)
      : worklist(worklist_), rowStatus(rowStatus_), colStatus(colStatus_), rowmap(rowmap_), entries(entries_), nv(nv_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t w) const
    {
      lno_t col = worklist(w);
      status_t minNeiStat = OUT_SET;
      size_type rowBegin = rowmap(col);
      size_type rowEnd = rowmap(col + 1);
      for(size_type j = rowBegin; j <= rowEnd; j++)
      {
        lno_t nei = (j == rowEnd) ? col : entries(j);
        if(nei >= nv)
          continue;
        status_t neiStat = rowStatus(nei);
        if(neiStat < minNeiStat)
          minNeiStat = neiStat;
      }
      if(minNeiStat == IN_SET)
        minNeiStat = OUT_SET;
      colStatus(col) = minNeiStat;
    }

    lno_view_t worklist;
    status_view_t rowStatus;
    status_view_t colStatus;
    rowmap_t rowmap;
    entries_t entries;
    lno_t nv;
  };

  struct InitWorklistFunctor
  {
    InitWorklistFunctor(const lno_view_t& worklist_)
      : worklist(worklist_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      worklist(i) = i;
    }
    lno_view_t worklist;
  };

  struct CountInSet
  {
    CountInSet(const status_view_t& rowStatus_)
      : rowStatus(rowStatus_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet) const
    {
      if(rowStatus(i) == IN_SET)
        lNumInSet++;
    }
    status_view_t rowStatus;
  };

  struct CompactInSet
  {
    CompactInSet(const status_view_t& rowStatus_, const lno_view_t& setList_)
      : rowStatus(rowStatus_), setList(setList_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(lno_t i, lno_t& lNumInSet, bool finalPass) const
    {
      if(rowStatus(i) == IN_SET)
      {
        if(finalPass)
          setList(lNumInSet) = i;
        lNumInSet++;
      }
    }
    status_view_t rowStatus;
    lno_view_t setList;
  };

  lno_view_t compute()
  {
    //Initialize first worklist to 0...numVerts
    Kokkos::parallel_for(range_pol(0, numVerts), InitWorklistFunctor(worklist1));
    lno_t workRemain = numVerts;
    int numIter = 0;
    while(workRemain)
    {
      //do another iteration
      Kokkos::parallel_for(range_pol(0, workRemain),
          IterateStatusFunctor(rowStatus, colStatus, rowmap, entries, numVerts, worklist1, colUpdateBitset));
      //And refresh the column statuses using the other worklist.
      lno_t colsToRefresh;
      Kokkos::parallel_scan(range_pol(0, numVerts),
          ColRefreshWorklist(colUpdateBitset, worklist2), colsToRefresh);
      Kokkos::parallel_for(range_pol(0, colsToRefresh),
          RefreshColStatus(worklist2, rowStatus, colStatus, rowmap, entries, numVerts));
      //then build the next worklist with a scan. Also get the length of the next worklist.
      lno_t newWorkRemain = 0;
      Kokkos::parallel_scan(range_pol(0, workRemain),
          UpdateWorklistFunctor(rowStatus, worklist1, worklist2),
          newWorkRemain);
      //Finally, flip the worklists
      std::swap(worklist1, worklist2);
      workRemain = newWorkRemain;
      numIter++;
    }
    //now that every vertex has been decided IN_SET/OUT_SET,
    //build a compact list of the vertices which are IN_SET.
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
  //The number of bits required to represent vertex IDs, in the ECL-MIS tiebreak scheme:
  //  ceil(log_2(numVerts + 1))
  int nvBits;
  lno_t minDegree;
  lno_t maxDegree;
  //Bitset representing columns whose status needs to be recomputed
  //These bits are cleared after each refresh.
  bitset_t colUpdateBitset;
  lno_view_t worklist1;
  lno_view_t worklist2;
};

template<typename device_t, typename rowmap_t, typename entries_t, typename labels_t>
struct D2_MIS_Coarsening
{
  using exec_space = typename device_t::execution_space;
  using mem_space = typename device_t::memory_space;
  using bitset_t = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;
  using size_type = typename rowmap_t::non_const_value_type;
  using lno_t = typename entries_t::non_const_value_type;
  using lno_view_t = typename entries_t::non_const_type;
  //The type of status/priority values.
  using status_t = typename std::make_unsigned<lno_t>::type;
  using status_view_t = Kokkos::View<status_t*, mem_space>;
  using range_pol = Kokkos::RangePolicy<exec_space>;

  D2_MIS_Coarsening(const rowmap_t& rowmap_, const entries_t& entries_, const labels_t& mis2_)
    : rowmap(rowmap_), entries(entries_), mis2(mis2_),
      numVerts(rowmap.extent(0) - 1),
      labels(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Cluster Labels"), numVerts)
  {
    Kokkos::deep_copy(labels, (lno_t) -1);
  }

  //Phase 1 (over 0...numClusters) labels roots and immediate neighbors of roots.
  struct Phase1Functor
  {
    Phase1Functor(const rowmap_t& rowmap_, const entries_t& entries_, const labels_t& mis2_, lno_t numVerts_, const labels_t& labels_)
      : rowmap(rowmap_), entries(entries_), mis2(mis2_), numVerts(numVerts_), labels(labels_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      lno_t root = mis2(i);
      size_type rowBegin = rowmap(root);
      size_type rowEnd = rowmap(root + 1);
      labels(root) = i;
      for(size_type j = rowBegin; j < rowEnd; j++)
      {
        lno_t nei = entries(j);
        if(nei != root && nei < numVerts)
        {
          labels(nei) = i;
        }
      }
    }

    rowmap_t rowmap;
    entries_t entries;
    labels_t mis2;
    lno_t numVerts;
    labels_t labels;
  };

  KOKKOS_INLINE_FUNCTION static uint32_t xorshiftHash(uint32_t in)
  {
    uint32_t x = in;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
  }

  //Phase 2 (over 0...numVerts) joins unlabeled vertices to the smallest adjacent cluster
  struct Phase2Functor
  {
    Phase2Functor(const rowmap_t& rowmap_, const entries_t& entries_, const labels_t& mis2_, lno_t numVerts_, const labels_t& labels_)
      : rowmap(rowmap_), entries(entries_), mis2(mis2_), numVerts(numVerts_), labels(labels_)
    {}

    KOKKOS_INLINE_FUNCTION void operator()(lno_t i) const
    {
      if(labels(i) != (lno_t) -1)
        return;
      size_type rowBegin = rowmap(i);
      size_type rowEnd = rowmap(i + 1);
      lno_t cluster = -1;
      uint32_t minScore = ~(uint32_t) 0;
      for(size_type j = rowBegin; j < rowEnd; j++)
      {
        lno_t nei = entries(j);
        if(nei == i || nei >= numVerts)
          continue;
        lno_t neiCluster = labels(nei);
        if(neiCluster != -1 && neiCluster != cluster)
        {
          //check if this cluster is smaller
          uint32_t score = xorshiftHash(i + xorshiftHash(neiCluster));
          if(score < minScore)
          {
            cluster = neiCluster;
            minScore = score;
          }
        }
      }
      labels(i) = cluster;
    }

    rowmap_t rowmap;
    entries_t entries;
    labels_t mis2;
    lno_t numVerts;
    labels_t labels;
  };

  labels_t compute()
  {
    lno_t numClusters = mis2.extent(0);
    Kokkos::parallel_for(range_pol(0, numClusters), Phase1Functor(rowmap, entries, mis2, numVerts, labels));
    Kokkos::parallel_for(range_pol(0, numVerts), Phase2Functor(rowmap, entries, mis2, numVerts, labels));
    return labels;
  }

  rowmap_t rowmap;
  entries_t entries;
  labels_t mis2;
  lno_t numVerts;
  labels_t labels;
};

}}}

#endif
