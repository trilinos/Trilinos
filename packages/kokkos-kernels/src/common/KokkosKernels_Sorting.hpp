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
#include "Kokkos_Atomic.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "impl/Kokkos_Timer.hpp"
#include <type_traits>

namespace KokkosKernels {
namespace Impl {

//Radix sort for integers, on a single thread within a team.
//Pros: few diverging branches, so OK for sorting on a single GPU vector lane. Better on CPU cores.
//Con: requires auxiliary storage, and this version only works for integers
template<typename Ordinal, typename ValueType>
KOKKOS_INLINE_FUNCTION void
SerialRadixSort(ValueType* values, ValueType* valuesAux, Ordinal n)
{
  static_assert(std::is_integral<ValueType>::value && std::is_unsigned<ValueType>::value,
      "radixSort can only be run on unsigned integers.");
  if(n <= 1)
    return;
  ValueType maxVal = 0;
  for(Ordinal i = 0; i < n; i++)
  {
    if(maxVal < values[i])
      maxVal = values[i];
  }
  //determine how many significant bits the data has
  int passes = 0;
  while(maxVal)
  {
    maxVal >>= 4;
    passes++;
  }
  //Is the data currently held in values (false) or valuesAux (true)?
  bool inAux = false;
  //sort 4 bits at a time, into 16 buckets
  ValueType mask = 0xF;
  //maskPos counts the low bit index of mask (0, 4, 8, ...)
  Ordinal maskPos = 0;
  for(int p = 0; p < passes; p++)
  {
    //Count the number of elements in each bucket
    Ordinal count[16] = {0};
    Ordinal offset[17];
    if(!inAux)
    {
      for(Ordinal i = 0; i < n; i++)
      {
        count[(values[i] & mask) >> maskPos]++;
      }
    }
    else
    {
      for(Ordinal i = 0; i < n; i++)
      {
        count[(valuesAux[i] & mask) >> maskPos]++;
      }
    }
    offset[0] = 0;
    //get offset as the prefix sum for count
    for(Ordinal i = 0; i < 16; i++)
    {
      offset[i + 1] = offset[i] + count[i];
    }
    //now for each element in [lo, hi), move it to its offset in the other buffer
    //this branch should be ok because whichBuf is the same on all threads
    if(!inAux)
    {
      for(Ordinal i = 0; i < n; i++)
      {
        Ordinal bucket = (values[i] & mask) >> maskPos;
        valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
        count[bucket]--;
      }
    }
    else
    {
      for(Ordinal i = 0; i < n; i++)
      {
        Ordinal bucket = (valuesAux[i] & mask) >> maskPos;
        values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
        count[bucket]--;
      }
    }
    inAux = !inAux;
    mask = mask << 4;
    maskPos += 4;
  }
  //Move values back into main array if they are currently in aux.
  //This is the case if an odd number of rounds were done.
  if(inAux)
  {
    for(Ordinal i = 0; i < n; i++)
    {
      values[i] = valuesAux[i];
    }
  }
}

//Radix sort for integers (no internal parallelism).
//While sorting, also permute "perm" array along with the values.
//Pros: few diverging branches, so good for sorting on a single GPU vector lane.
//Con: requires auxiliary storage, this version only works for integers (although float/double is possible)
template<typename Ordinal, typename ValueType, typename PermType>
KOKKOS_INLINE_FUNCTION void
SerialRadixSort2(ValueType* values, ValueType* valuesAux, PermType* perm, PermType* permAux, Ordinal n)
{
  static_assert(std::is_integral<ValueType>::value && std::is_unsigned<ValueType>::value,
      "radixSort can only be run on unsigned integers.");
  if(n <= 1)
    return;
  ValueType maxVal = 0;
  for(Ordinal i = 0; i < n; i++)
  {
    if(maxVal < values[i])
      maxVal = values[i];
  }
  int passes = 0;
  while(maxVal)
  {
    maxVal >>= 4;
    passes++;
  }
  //Is the data currently held in values (false) or valuesAux (true)?
  bool inAux = false;
  //sort 4 bits at a time, into 16 buckets
  ValueType mask = 0xF;
  //maskPos counts the low bit index of mask (0, 4, 8, ...)
  Ordinal maskPos = 0;
  for(int p = 0; p < passes; p++)
  {
    //Count the number of elements in each bucket
    Ordinal count[16] = {0};
    Ordinal offset[17];
    if(!inAux)
    {
      for(Ordinal i = 0; i < n; i++)
      {
        count[(values[i] & mask) >> maskPos]++;
      }
    }
    else
    {
      for(Ordinal i = 0; i < n; i++)
      {
        count[(valuesAux[i] & mask) >> maskPos]++;
      }
    }
    offset[0] = 0;
    //get offset as the prefix sum for count
    for(Ordinal i = 0; i < 16; i++)
    {
      offset[i + 1] = offset[i] + count[i];
    }
    //now for each element in [lo, hi), move it to its offset in the other buffer
    //this branch should be ok because whichBuf is the same on all threads
    if(!inAux)
    {
      for(Ordinal i = 0; i < n; i++)
      {
        Ordinal bucket = (values[i] & mask) >> maskPos;
        valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
        permAux[offset[bucket + 1] - count[bucket]] = perm[i];
        count[bucket]--;
      }
    }
    else
    {
      for(Ordinal i = 0; i < n; i++)
      {
        Ordinal bucket = (valuesAux[i] & mask) >> maskPos;
        values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
        perm[offset[bucket + 1] - count[bucket]] = permAux[i];
        count[bucket]--;
      }
    }
    inAux = !inAux;
    mask = mask << 4;
    maskPos += 4;
  }
  //Move values back into main array if they are currently in aux.
  //This is the case if an odd number of rounds were done.
  if(inAux)
  {
    for(Ordinal i = 0; i < n; i++)
    {
      values[i] = valuesAux[i];
      perm[i] = permAux[i];
    }
  }
}

template<typename Value>
struct DefaultComparator
{
  KOKKOS_INLINE_FUNCTION bool operator()(const Value lhs, const Value rhs) const
  {
    return lhs < rhs;
  }
};

//Bitonic merge sort (requires only comparison operators and trivially-copyable)
//Pros: In-place, plenty of parallelism for GPUs, and memory references are coalesced
//Con: O(n log^2(n)) serial time is bad on CPUs
//Good diagram of the algorithm at https://en.wikipedia.org/wiki/Bitonic_sorter
template<typename Ordinal, typename ValueType, typename TeamMember, typename Comparator = DefaultComparator<ValueType>>
KOKKOS_INLINE_FUNCTION void
TeamBitonicSort(ValueType* values, Ordinal n, const TeamMember mem, const Comparator& comp = Comparator())
{
  //Algorithm only works on power-of-two input size only.
  //If n is not a power-of-two, will implicitly pretend
  //that values[i] for i >= n is just the max for ValueType, so it never gets swapped
  Ordinal npot = 1;
  Ordinal levels = 0;
  while(npot < n)
  {
    levels++;
    npot <<= 1;
  }
  for(Ordinal i = 0; i < levels; i++)
  {
    for(Ordinal j = 0; j <= i; j++)
    {
      // n/2 pairs of items are compared in parallel
      Kokkos::parallel_for(Kokkos::TeamVectorRange(mem, npot / 2),
        [=](const Ordinal t)
        {
          //How big are the brown/pink boxes?
          Ordinal boxSize = Ordinal(2) << (i - j);
          //Which box contains this thread?
          Ordinal boxID = t >> (i - j);         //t * 2 / boxSize;
          Ordinal boxStart = boxID << (1 + i - j);  //boxID * boxSize
          Ordinal boxOffset = t - (boxStart >> 1);  //t - boxID * boxSize / 2;
          Ordinal elem1 = boxStart + boxOffset;
          if(j == 0)
          {
            //first phase (brown box): within a block, compare with the opposite value in the box
            Ordinal elem2 = boxStart + boxSize - 1 - boxOffset;
            if(elem2 < n)
            {
              //both elements in bounds, so compare them and swap if out of order
              if(comp(values[elem2], values[elem1]))
              {
                ValueType temp = values[elem1];
                values[elem1] = values[elem2];
                values[elem2] = temp;
              }
            }
          }
          else
          {
            //later phases (pink box): within a block, compare with fixed distance (boxSize / 2) apart
            Ordinal elem2 = elem1 + boxSize / 2;
            if(elem2 < n)
            {
              if(comp(values[elem2], values[elem1]))
              {
                ValueType temp = values[elem1];
                values[elem1] = values[elem2];
                values[elem2] = temp;
              }
            }
          }
        });
      mem.team_barrier();
    }
  }
}

//Sort "values", while applying the same swaps to "perm" 
template<typename Ordinal, typename ValueType, typename PermType, typename TeamMember, typename Comparator = DefaultComparator<ValueType>>
KOKKOS_INLINE_FUNCTION void
TeamBitonicSort2(ValueType* values, PermType* perm, Ordinal n, const TeamMember mem, const Comparator& comp = Comparator())
{
  //Algorithm only works on power-of-two input size only.
  //If n is not a power-of-two, will implicitly pretend
  //that values[i] for i >= n is just the max for ValueType, so it never gets swapped
  Ordinal npot = 1;
  Ordinal levels = 0;
  while(npot < n)
  {
    levels++;
    npot <<= 1;
  }
  for(Ordinal i = 0; i < levels; i++)
  {
    for(Ordinal j = 0; j <= i; j++)
    {
      // n/2 pairs of items are compared in parallel
      Kokkos::parallel_for(Kokkos::TeamVectorRange(mem, npot / 2),
        [=](const Ordinal t)
        {
          //How big are the brown/pink boxes?
          Ordinal boxSize = Ordinal(2) << (i - j);
          //Which box contains this thread?
          Ordinal boxID = t >> (i - j);         //t * 2 / boxSize;
          Ordinal boxStart = boxID << (1 + i - j);  //boxID * boxSize
          Ordinal boxOffset = t - (boxStart >> 1);  //t - boxID * boxSize / 2;
          Ordinal elem1 = boxStart + boxOffset;
          if(j == 0)
          {
            //first phase (brown box): within a block, compare with the opposite value in the box
            Ordinal elem2 = boxStart + boxSize - 1 - boxOffset;
            if(elem2 < n)
            {
              //both elements in bounds, so compare them and swap if out of order
              if(comp(values[elem2], values[elem1]))
              {
                ValueType temp1 = values[elem1];
                values[elem1] = values[elem2];
                values[elem2] = temp1;
                PermType temp2 = perm[elem1];
                perm[elem1] = perm[elem2];
                perm[elem2] = temp2;
              }
            }
          }
          else
          {
            //later phases (pink box): within a block, compare with fixed distance (boxSize / 2) apart
            Ordinal elem2 = elem1 + boxSize / 2;
            if(elem2 < n)
            {
              if(comp(values[elem2], values[elem1]))
              {
                ValueType temp1 = values[elem1];
                values[elem1] = values[elem2];
                values[elem2] = temp1;
                PermType temp2 = perm[elem1];
                perm[elem1] = perm[elem2];
                perm[elem2] = temp2;
              }
            }
          }
        });
      mem.team_barrier();
    }
  }
}

//Functor that sorts a view on one team
template<typename View, typename Ordinal, typename TeamMember, typename Comparator>
struct BitonicSingleTeamFunctor
{
  BitonicSingleTeamFunctor(View& v_, const Comparator& comp_) : v(v_), comp(comp_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
  {
    TeamBitonicSort<Ordinal, typename View::value_type, TeamMember, Comparator>(v.data(), v.extent(0), t, comp);
  };
  View v;
  Comparator comp;
};

//Functor that sorts equally sized chunks on each team
template<typename View, typename Ordinal, typename TeamMember, typename Comparator>
struct BitonicChunkFunctor
{
  BitonicChunkFunctor(View& v_, const Comparator& comp_, Ordinal chunkSize_) : v(v_), comp(comp_), chunkSize(chunkSize_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
  {
    Ordinal chunk = t.league_rank();
    Ordinal chunkStart = chunk * chunkSize;
    Ordinal n = chunkSize;
    if(chunkStart + n > Ordinal(v.extent(0)))
      n = v.extent(0) - chunkStart;
    TeamBitonicSort<Ordinal, typename View::value_type, TeamMember, Comparator>(v.data() + chunkStart, n, t, comp);
  };
  View v;
  Comparator comp;
  Ordinal chunkSize;
};

//Functor that does just the first phase (brown) of bitonic sort on equally-sized chunks
template<typename View, typename Ordinal, typename TeamMember, typename Comparator>
struct BitonicPhase1Functor
{
  typedef typename View::value_type Value;
  BitonicPhase1Functor(View& v_, const Comparator& comp_, Ordinal boxSize_, Ordinal teamsPerBox_)
    : v(v_), comp(comp_), boxSize(boxSize_), teamsPerBox(teamsPerBox_)
  {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
  {
    Ordinal box = t.league_rank() / teamsPerBox;
    Ordinal boxStart = boxSize * box;
    Ordinal work = boxSize / teamsPerBox / 2;
    Ordinal workStart = work * (t.league_rank() % teamsPerBox);
    Ordinal workReflect = boxSize - workStart - 1;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(t, work),
      [&](const Ordinal i)
      {
        Ordinal elem1 = boxStart + workStart + i;
        Ordinal elem2 = boxStart + workReflect - i;
        if(elem2 < Ordinal(v.extent(0)))
        {
          if(comp(v(elem2), v(elem1)))
          {
            Value temp = v(elem1);
            v(elem1) = v(elem2);
            v(elem2) = temp;
          }
        }
      });
  };
  View v;
  Comparator comp;
  Ordinal boxSize;
  Ordinal teamsPerBox;
};

//Functor that does the second phase (red) of bitonic sort
template<typename View, typename Ordinal, typename TeamMember, typename Comparator>
struct BitonicPhase2Functor
{
  typedef typename View::value_type Value;
  BitonicPhase2Functor(View& v_, const Comparator& comp_, Ordinal boxSize_, Ordinal teamsPerBox_)
    : v(v_), comp(comp_), boxSize(boxSize_), teamsPerBox(teamsPerBox_)
  {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
  {
    Ordinal logBoxSize = 1;
    while((Ordinal(1) << logBoxSize) < boxSize)
      logBoxSize++;
    Ordinal box = t.league_rank() / teamsPerBox;
    Ordinal boxStart = boxSize * box;
    Ordinal work = boxSize / teamsPerBox / 2;
    Ordinal workStart = boxStart + work * (t.league_rank() % teamsPerBox);
    Ordinal jump = boxSize / 2;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(t, work),
      [&](const Ordinal i)
      {
        Ordinal elem1 = workStart + i;
        Ordinal elem2 = workStart + jump + i;
        if(elem2 < Ordinal(v.extent(0)))
        {
          if(comp(v(elem2), v(elem1)))
          {
            Value temp = v(elem1);
            v(elem1) = v(elem2);
            v(elem2) = temp;
          }
        }
      });
    if(teamsPerBox == 1)
    {
      //This team can finish phase 2 for all the smaller red boxes that follow,
      //since there are no longer cross-team data dependencies
      for(Ordinal subLevel = 1; subLevel < logBoxSize; subLevel++)
      {
        t.team_barrier();
        Ordinal logSubBoxSize = logBoxSize - subLevel;
        Ordinal subBoxSize = Ordinal(1) << logSubBoxSize;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(t, work),
          [&](const Ordinal i)
          {
            Ordinal globalThread = i + t.league_rank() * work;
            Ordinal subBox = globalThread >> (logSubBoxSize - 1);
            Ordinal subBoxStart = subBox << logSubBoxSize;
            Ordinal subBoxOffset = globalThread & ((Ordinal(1) << (logSubBoxSize - 1)) - 1); //i % (subBoxSize / 2)
            Ordinal elem1 = subBoxStart + subBoxOffset;
            //later phases (pink box): within a block, compare with fixed distance (boxSize / 2) apart
            Ordinal elem2 = elem1 + subBoxSize / 2;
            if(elem2 < Ordinal(v.extent(0)))
            {
              if(comp(v(elem2), v(elem1)))
              {
                Value temp = v(elem1);
                v(elem1) = v(elem2);
                v(elem2) = temp;
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

//Version to be called from host on a single array
//Generally ~2x slower than Kokkos::sort() for large arrays (> 50 M elements),
//but faster for smaller arrays.
//
//This is more general than BinSort: bitonic supports any trivially copyable type
//and an arbitrary device-compatible comparison operator (provided through operator() of Comparator)
//If comparator is void, use operator< (which should only be used for primitives)
template<typename View, typename ExecSpace, typename Ordinal, typename Comparator = DefaultComparator<typename View::value_type>>
void bitonicSort(View v, const Comparator& comp = Comparator())
{
  typedef Kokkos::TeamPolicy<ExecSpace> team_policy;
  typedef typename team_policy::member_type team_member;
  Ordinal n = v.extent(0);
  //If n is small, just sort on a single team
  if(n <= Ordinal(1) << 12)
  {
    Kokkos::parallel_for(team_policy(1, Kokkos::AUTO()),
        BitonicSingleTeamFunctor<View, Ordinal, team_member, Comparator>(v, comp));
  }
  else
  {
    Ordinal npot = 1;
    while(npot < n)
      npot <<= 1;
    //Partition the data equally among fixed number of teams
    Ordinal chunkSize = 512;
    Ordinal numTeams = npot / chunkSize;
    //First, sort within teams
    Kokkos::parallel_for(team_policy(numTeams, Kokkos::AUTO()),
        BitonicChunkFunctor<View, Ordinal, team_member, Comparator>(v, comp, chunkSize));
    for(int teamsPerBox = 2; teamsPerBox <= npot / chunkSize; teamsPerBox *= 2)
    {
      Ordinal boxSize = teamsPerBox * chunkSize;
      Kokkos::parallel_for(team_policy(numTeams, Kokkos::AUTO()),
          BitonicPhase1Functor<View, Ordinal, team_member, Comparator>(v, comp, boxSize, teamsPerBox));
      for(int boxDiv = 1; teamsPerBox >> boxDiv; boxDiv++)
      {
        Kokkos::parallel_for(team_policy(numTeams, Kokkos::AUTO()),
            BitonicPhase2Functor<View, Ordinal, team_member, Comparator>(v, comp, boxSize >> boxDiv, teamsPerBox >> boxDiv));
      }
    }
  }
}

}}

#endif

