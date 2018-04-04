/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifndef _KOKKOS_SPADD_HPP
#define _KOKKOS_SPADD_HPP

#include "KokkosKernels_Handle.hpp"

namespace KokkosSparse {
namespace Experimental {

  /*
  Unsorted symbolic algorithm notes:
  -Only needs to sort and merge indices once, in symbolic (sorting is expensive)
  -Can't afford to allocate dense Views for indices/values (assume number of columns is very large)
  -Want numeric() to know exactly where each A/B entry belongs in Ccolinds/Cvalues
  -To accomplish all of these, symbolic() computes arrays Apos and Bpos (both are type clno_nnz_view_t_,
    and have same length as a_entries and b_entries respectively)
    -Apos/Bpos are saved in the handle
  -Apos and Bpos each contain the final index within C row where the A/B entry belongs
  -See UnsortedNumericSumFunctor below for the usage of Apos/Bpos
  */

//Helper macro to check that two types are the same (ignoring const)
#define SAME_TYPE(A, B) std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

  //get C rowmap for sorted input
  template<typename size_type, typename ordinal_type, typename ARowPtrsT, typename BRowPtrsT, typename AColIndsT, typename BColIndsT, typename CRowPtrsT>
  struct SortedCountEntries
  {
    SortedCountEntries(
        const typename ARowPtrsT::const_type Arowptrs_,
        const AColIndsT Acolinds_,
        const typename BRowPtrsT::const_type Browptrs_,
        const BColIndsT Bcolinds_,
        CRowPtrsT Crowcounts_) :
      Arowptrs(Arowptrs_), Acolinds(Acolinds_),
      Browptrs(Browptrs_), Bcolinds(Bcolinds_),
      Crowcounts(Crowcounts_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      //count the union of nonzeros in Arow and Brow
      size_type numEntries = 0;
      size_type ai = 0;
      size_type bi = 0;
      size_type Arowstart = Arowptrs(i);
      size_type Arowlen = Arowptrs(i + 1) - Arowstart;
      size_type Browstart = Browptrs(i);
      size_type Browlen = Browptrs(i + 1) - Browstart;
      while (ai < Arowlen && bi < Browlen)
      {
        //have an entry in C's row
        numEntries++;
        ordinal_type Acol = Acolinds(Arowstart + ai);
        ordinal_type Bcol = Bcolinds(Browstart + bi);
        if(Acol <= Bcol)
          ai++;
        if(Acol >= Bcol)
          bi++;
      }
      //if a and b have different numbers of entries in row i,
      //the next two lines will account for remaining entries in the longer one
      numEntries += Arowlen - ai;
      numEntries += Browlen - bi;
      Crowcounts(i) = numEntries;
    }
    const typename ARowPtrsT::const_type Arowptrs;
    const AColIndsT Acolinds;
    const typename BRowPtrsT::const_type Browptrs;
    const BColIndsT Bcolinds;
    CRowPtrsT Crowcounts;
  };

  //get upper bound for C entries per row (assumes worst case, that entries in A and B on each row are disjoint)
  template<typename size_type, typename ARowPtrsT, typename BRowPtrsT, typename CRowPtrsT>
  struct UnsortedEntriesUpperBound
  {
      UnsortedEntriesUpperBound(
        const typename ARowPtrsT::const_type Arowptrs_,
        const typename BRowPtrsT::const_type Browptrs_,
        CRowPtrsT Crowcounts_) :
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowcounts(Crowcounts_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      Crowcounts(i) = (Arowptrs(i + 1) - Arowptrs(i)) + (Browptrs(i + 1) - Browptrs(i));
    }
    const typename ARowPtrsT::const_type Arowptrs;
    const typename BRowPtrsT::const_type Browptrs;
    CRowPtrsT Crowcounts;
  };

  //Unsorted symbolic: new functors:
  //  -compute uncompressed C (entries only, no values)
  //  -sort uncompressed C entries within row, while permuting A union B permutation array
  //  -compress sorted C entries and A,B perm arrays at the same time, which produces Crowcounts value
  //Inputs: A, B rowptrs/colinds, C uncompressed rowptrs (and allocated C entries)
  //Output: C uncompressed colinds
  template<typename size_type, typename ordinal_type,
                             typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
                             typename AcolindsT, typename BcolindsT, typename CcolindsT>
  struct UnmergedSumFunctor
  {
    UnmergedSumFunctor(const ArowptrsT Arowptrs_, const AcolindsT Acolinds_,
                       const BrowptrsT Browptrs_, const BcolindsT Bcolinds_,
                       CrowptrsT Crowptrs_, CcolindsT Ccolinds_, CcolindsT ABperm_) :
      Arowptrs(Arowptrs_), Acolinds(Acolinds_),
      Browptrs(Browptrs_), Bcolinds(Bcolinds_),
      Crowptrs(Crowptrs_), Ccolinds(Ccolinds_), ABperm(ABperm_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      size_type inserted = 0;
      size_type crowstart = Crowptrs(i);
      size_type arowstart = Arowptrs(i);
      size_type arowlen = Arowptrs(i + 1) - arowstart;
      size_type browstart = Browptrs(i);
      size_type browlen = Browptrs(i + 1) - browstart;
      //Insert all A entries, then all B entries
      for(size_type j = 0; j < arowlen; j++)
      {
        Ccolinds(crowstart + inserted) = Acolinds(arowstart + j);
        ABperm(crowstart + inserted) = j;
        inserted++;
      }
      for(size_type j = 0; j < browlen; j++)
      {
        Ccolinds(crowstart + inserted) = Bcolinds(browstart + j);
        //tell A and B permutation values apart by adding arowlen as a bias to B values
        ABperm(crowstart + inserted) = j + arowlen;
        inserted++;
      }
    }
    const ArowptrsT Arowptrs;
    const AcolindsT Acolinds;
    const BrowptrsT Browptrs;
    const BcolindsT Bcolinds;
    const CrowptrsT Crowptrs;
    CcolindsT Ccolinds;
    CcolindsT ABperm;
  };

  template<typename size_type, typename KeyType, typename ValueType, typename IndexType>
  KOKKOS_INLINE_FUNCTION void
  radixSortKeysAndValues(KeyType* keys, KeyType* keysAux, ValueType* values, ValueType* valuesAux, IndexType n)
  {
    if(n <= 1)
      return;
    //sort 4 bits at a time
    KeyType mask = 0xF;
    bool inAux = false;
    //maskPos counts the low bit index of mask (0, 4, 8, ...)
    IndexType maskPos = 0;
    IndexType sortBits = 0;
    KeyType minKey = keys[0];
    KeyType maxKey = keys[0];
    for(size_type i = 0; i < n; i++)
    {
      if(keys[i] < minKey)
        minKey = keys[i];
      if(keys[i] > maxKey)
        maxKey = keys[i];
    }
    //subtract a bias of minKey so that key range starts at 0
    for(size_type i = 0; i < n; i++)
    {
      keys[i] -= minKey;
    }
    KeyType upperBound = maxKey - minKey;
    while(upperBound)
    {
      upperBound >>= 1;
      sortBits++;
    }
    for(size_type s = 0; s < (sortBits + 3) / 4; s++)
    {
      //Count the number of elements in each bucket
      IndexType count[16] = {0};
      IndexType offset[17] = {0};
      if(!inAux)
      {
        for(IndexType i = 0; i < n; i++)
        {
          count[(keys[i] & mask) >> maskPos]++;
        }
      }
      else
      {
        for(IndexType i = 0; i < n; i++)
        {
          count[(keysAux[i] & mask) >> maskPos]++;
        }
      }
      //get offset as the prefix sum for count
      for(IndexType i = 0; i < 16; i++)
      {
        offset[i + 1] = offset[i] + count[i];
      }
      //now for each element in [lo, hi), move it to its offset in the other buffer
      //this branch should be ok because whichBuf is the same on all threads
      if(!inAux)
      {
        //copy from *Over to *Aux
        for(IndexType i = 0; i < n; i++)
        {
          IndexType bucket = (keys[i] & mask) >> maskPos;
          keysAux[offset[bucket + 1] - count[bucket]] = keys[i];
          valuesAux[offset[bucket + 1] - count[bucket]] = values[i];
          count[bucket]--;
        }
      }
      else
      {
        //copy from *Aux to *Over
        for(IndexType i = 0; i < n; i++)
        {
          IndexType bucket = (keysAux[i] & mask) >> maskPos;
          keys[offset[bucket + 1] - count[bucket]] = keysAux[i];
          values[offset[bucket + 1] - count[bucket]] = valuesAux[i];
          count[bucket]--;
        }
      }
      inAux = !inAux;
      mask = mask << 4;
      maskPos += 4;
    }
    if(inAux)
    {
      //need to deep copy from aux arrays to main
      for(IndexType i = 0; i < n; i++)
      {
        keys[i] = keysAux[i];
        values[i] = valuesAux[i];
      }
    }
    //remove the bias
    for(size_type i = 0; i < n; i++)
    {
      keys[i] += minKey;
    }
  }

  template<typename size_type, typename CrowptrsT, typename CcolindsT>
  struct SortEntriesFunctor
  {
    SortEntriesFunctor(const CrowptrsT Crowptrs_, CcolindsT Ccolinds_, CcolindsT ABperm_) :
      Crowptrs(Crowptrs_),
      Ccolinds(Ccolinds_),
      CcolindsAux("C colind aux", Ccolinds_.dimension_0()),
      ABperm(ABperm_),
      ABpermAux("AB perm aux", ABperm_.dimension_0())
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      //3: Sort each row's colinds (permuting values at same time), then count unique colinds (write that to Crowptr(i))
      //CrowptrTemp tells how many entries in each oversized row
      size_type rowStart = Crowptrs(i);
      size_type rowEnd = Crowptrs(i + 1);
      size_type rowNum = rowEnd - rowStart;
      radixSortKeysAndValues<size_type, typename CcolindsT::non_const_value_type, typename CcolindsT::non_const_value_type>
        (Ccolinds.ptr_on_device() + rowStart, CcolindsAux.ptr_on_device() + rowStart,
         ABperm.ptr_on_device() + rowStart, ABpermAux.ptr_on_device() + rowStart, rowNum);
    }
    CrowptrsT Crowptrs;
    CcolindsT Ccolinds;
    CcolindsT CcolindsAux;
    CcolindsT ABperm;
    CcolindsT ABpermAux;
  };

  template<typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT, typename CcolindsT>
  struct MergeEntriesFunctor
  {
    MergeEntriesFunctor(const ArowptrsT Arowptrs_, const BrowptrsT Browptrs_, const CrowptrsT Crowptrs_, CrowptrsT Crowcounts_,
        const CcolindsT Ccolinds_, const CcolindsT ABperm_, CcolindsT Apos_, CcolindsT Bpos_) :
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowptrs(Crowptrs_),
      Crowcounts(Crowcounts_),
      Ccolinds(Ccolinds_),
      ABperm(ABperm_),
      Apos(Apos_),
      Bpos(Bpos_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      size_type CrowStart = Crowptrs(i);
      size_type CrowEnd = Crowptrs(i + 1);
      size_type ArowStart = Arowptrs(i);
      size_type ArowNum = Arowptrs(i + 1) - ArowStart;
      size_type BrowStart = Browptrs(i);
      ordinal_type CFit = 0; //counting through merged C indices (within row)
      for(size_type Cit = CrowStart; Cit < CrowEnd; Cit++)
      {
        size_type permVal = ABperm(Cit);
        if(permVal < ArowNum)
        {
          //Entry belongs to A
          ordinal_type Aindex = permVal;
          //The Aindex'th entry in row i of A will be added into the CFit'th entry in C
          Apos(ArowStart + Aindex) = CFit;
        }
        else
        {
          //Entry belongs to B
          ordinal_type Bindex = permVal - ArowNum;
          //The Bindex'th entry in row i of B will be added into the CFit'th entry in C
          Bpos(BrowStart + Bindex) = CFit;
        }
        //if NOT merging uncompressed entries Cit and Cit + 1, increment compressed index CFit
        bool mergingWithNext = Cit < CrowEnd - 1 && Ccolinds(Cit) == Ccolinds(Cit + 1);
        if(!mergingWithNext)
          CFit++;
      }
      //at end of the row, know how many entries are in merged C
      Crowcounts(i) = CFit;
    }
    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    CrowptrsT Crowcounts;
    CcolindsT Ccolinds;
    const CcolindsT ABperm;
    CcolindsT Apos;
    CcolindsT Bpos;
  };

  //from tpetra
  template <typename size_type, typename view_type>
  struct parallel_prefix_sum{
    view_type input;
    view_type output;
    typedef typename view_type::non_const_value_type value_type;
    parallel_prefix_sum(view_type in, view_type out): input(in), output(out) {}
    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type& i, value_type& update, const bool fin) const
    {
      size_type n = input.dimension_0();
      value_type curVal = (i < n) ? input(i) : 0;
      if(fin)
      {
        output(i) = update;
      }
      update += (i < n) ? curVal : 0;
    }
  };

  //Symbolic: count entries in each row in C to produce rowmap
  //kernel handle has information about whether it is sorted add or not.
  template <typename KernelHandle,
            typename alno_row_view_t_,
            typename alno_nnz_view_t_,
            typename blno_row_view_t_,
            typename blno_nnz_view_t_,
            typename clno_row_view_t_,
            typename clno_nnz_view_t_>
  void spadd_symbolic(
      KernelHandle* handle, 
      const alno_row_view_t_ a_rowmap,
      const alno_nnz_view_t_ a_entries,
      const blno_row_view_t_ b_rowmap,
      const blno_nnz_view_t_ b_entries,
      clno_row_view_t_ c_rowmap)    //c_rowmap must already be allocated (doesn't need to be initialized)
  {
    typedef typename KernelHandle::SPADDHandleType::execution_space execution_space;
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;
    //Check that A/B/C data types match KernelHandle types, and that C data types are nonconst (doesn't matter if A/B types are const)
    static_assert(SAME_TYPE(typename alno_row_view_t_::non_const_value_type, size_type),
        "add_symbolic: A size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(SAME_TYPE(typename blno_row_view_t_::non_const_value_type, size_type),
        "add_symbolic: B size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(SAME_TYPE(typename clno_row_view_t_::non_const_value_type, size_type),
        "add_symbolic: C size_type must match KernelHandle size_type)");
    static_assert(std::is_same<typename clno_row_view_t_::non_const_value_type, typename clno_row_view_t_::value_type>::value,
        "add_symbolic: C size_type must not be const");
    static_assert(SAME_TYPE(typename alno_nnz_view_t_::non_const_value_type, ordinal_type),
        "add_symbolic: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename blno_nnz_view_t_::non_const_value_type, ordinal_type),
        "add_symbolic: B entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename clno_nnz_view_t_::non_const_value_type, ordinal_type),
        "add_symbolic: C entry type must match KernelHandle entry type (aka nnz_lno_t)");
    static_assert(std::is_same<typename clno_row_view_t_::non_const_value_type, typename clno_row_view_t_::value_type>::value,
        "add_symbolic: C entry type must not be const");
    //symbolic just needs to compute c_rowmap
    //easy for sorted, but for unsorted is easiest to just compute the whole sum
    auto addHandle = handle->get_spadd_handle();
    auto nrows = a_rowmap.dimension_0() - 1;
    typedef Kokkos::RangePolicy<execution_space, ordinal_type> range_type;
    if(addHandle->is_input_sorted())
    {
      clno_row_view_t_ c_rowcounts("C row counts", nrows);
      //call entry count functor to get entry counts per row
      SortedCountEntries<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, alno_nnz_view_t_, blno_nnz_view_t_, clno_row_view_t_>
        countEntries(a_rowmap, a_entries, b_rowmap, b_entries, c_rowcounts);
      Kokkos::parallel_for(range_type(0, nrows), countEntries);
      execution_space::fence();
      //get c_rowmap as cumulative sum
      parallel_prefix_sum<size_type, clno_row_view_t_> prefix(c_rowcounts, c_rowmap);
      Kokkos::parallel_scan(range_type(0, nrows + 1), prefix);
      execution_space::fence();
    }
    else
    {
      //note: scoping individual parts of the process to free views sooner, minimizing peak memory usage
      //run the unsorted c_rowmap upper bound functor (just adds together A and B entry counts row by row)
      clno_row_view_t_ c_rowmap_upperbound("C row counts upper bound", nrows + 1);
      size_t c_nnz_upperbound = 0;
      {
        clno_row_view_t_ c_rowcounts_upperbound("C row counts upper bound", nrows);
        UnsortedEntriesUpperBound<size_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_>
          countEntries(a_rowmap, b_rowmap, c_rowcounts_upperbound);
        Kokkos::parallel_for(range_type(0, nrows), countEntries);
        execution_space::fence();
        //get (temporary) c_rowmap as cumulative sum
        parallel_prefix_sum<size_type, clno_row_view_t_> prefix(c_rowcounts_upperbound, c_rowmap_upperbound);
        Kokkos::parallel_scan(range_type(0, nrows + 1), prefix);
        //compute uncompressed entries of C (just indices, no scalars)
        execution_space::fence();

        auto d_c_nnz_size = Kokkos::subview(c_rowmap_upperbound, nrows);
        auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
        Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
        execution_space::fence();
        c_nnz_upperbound = h_c_nnz_size();
      }
      clno_nnz_view_t_ c_entries_uncompressed("C entries uncompressed", c_nnz_upperbound);
      clno_nnz_view_t_ ab_perm("A and B permuted entry indices", c_nnz_upperbound);
      //compute the unmerged sum
      UnmergedSumFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                         alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_> unmergedSum(
                         a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap_upperbound, c_entries_uncompressed, ab_perm);
      Kokkos::parallel_for(range_type(0, nrows), unmergedSum);
      execution_space::fence();
      //sort the unmerged sum
      SortEntriesFunctor<size_type, clno_row_view_t_, clno_nnz_view_t_>
        sortEntries(c_rowmap_upperbound, c_entries_uncompressed, ab_perm);
      Kokkos::parallel_for(range_type(0, nrows), sortEntries);
      execution_space::fence();
      clno_nnz_view_t_ a_pos("A entry positions", a_entries.dimension_0());
      clno_nnz_view_t_ b_pos("B entry positions", b_entries.dimension_0());
      //merge the entries and compute Apos/Bpos, as well as Crowcounts
      {
        clno_row_view_t_ c_rowcounts("C row counts", nrows);
        MergeEntriesFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_, clno_nnz_view_t_>
          mergeEntries(a_rowmap, b_rowmap, c_rowmap_upperbound, c_rowcounts, c_entries_uncompressed, ab_perm, a_pos, b_pos);
        Kokkos::parallel_for(range_type(0, nrows), mergeEntries);
        execution_space::fence();
        //compute actual c_rowmap
        parallel_prefix_sum<size_type, clno_row_view_t_> prefix(c_rowcounts, c_rowmap);
        Kokkos::parallel_scan(range_type(0, nrows + 1), prefix);
        execution_space::fence();
      }
      addHandle->set_a_b_pos(a_pos, b_pos);
    }
    //provide the number of NNZ in C to user through handle
    //addHandle->set_max_result_nnz(c_rowmap(nrows));

    auto d_c_nnz_size = Kokkos::subview(c_rowmap, nrows);
    auto h_c_nnz_size = Kokkos::create_mirror_view (d_c_nnz_size);
    Kokkos::deep_copy (h_c_nnz_size, d_c_nnz_size);
    execution_space::fence();
    size_type cmax = h_c_nnz_size();
    addHandle->set_max_result_nnz(cmax);


    addHandle->set_call_symbolic();
    addHandle->set_call_numeric(false);
  }

  template<typename size_type, typename ordinal_type,
           typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
           typename AcolindsT, typename BcolindsT, typename CcolindsT,
           typename AvaluesT, typename BvaluesT, typename CvaluesT,
           typename AscalarT, typename BscalarT>
  struct SortedNumericSumFunctor
  {
    SortedNumericSumFunctor(const ArowptrsT Arowptrs_, const BrowptrsT Browptrs_, const CrowptrsT Crowptrs_,
    const AcolindsT Acolinds_, const BcolindsT Bcolinds_, CcolindsT Ccolinds_,
    const AvaluesT Avalues_, const BvaluesT Bvalues_, CvaluesT Cvalues_,
    const AscalarT alpha_, const BscalarT beta_) :
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowptrs(Crowptrs_),
      Acolinds(Acolinds_),
      Bcolinds(Bcolinds_),
      Ccolinds(Ccolinds_),
      Avalues(Avalues_),
      Bvalues(Bvalues_),
      Cvalues(Cvalues_),
      alpha(alpha_),
      beta(beta_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      size_type CrowStart = Crowptrs(i);
      size_type ArowStart = Arowptrs(i);
      size_type ArowEnd = Arowptrs(i + 1);
      size_type Arowlen = ArowEnd - ArowStart;
      size_type BrowStart = Browptrs(i);
      size_type BrowEnd = Browptrs(i + 1);
      size_type Browlen = BrowEnd - BrowStart;
      size_type ai = 0;
      size_type bi = 0;
      size_type ci = 0;
      //add in A entries, while setting C colinds
      while(ai < Arowlen && bi < Browlen)
      {
        ordinal_type Acol = Acolinds(ArowStart + ai);
        ordinal_type Bcol = Bcolinds(BrowStart + bi);
        if(Acol <= Bcol)
        {
          Ccolinds(CrowStart + ci) = Acol;
          Cvalues(CrowStart + ci) += alpha * Avalues(ArowStart + ai);
          ai++;
        }
        if(Acol >= Bcol)
        {
          Ccolinds(CrowStart + ci) = Bcol;
          Cvalues(CrowStart + ci) += beta * Bvalues(BrowStart + bi);
          bi++;
        }
        ci++;
      }
      //append remaining A entries (if any)
      while(ai < Arowlen)
      {
        Ccolinds(CrowStart + ci) = Acolinds(ArowStart + ai);
        Cvalues(CrowStart + ci) = alpha * Avalues(ArowStart + ai);
        ai++;
        ci++;
      }
      //append remaining B entries (if any)
      while(bi < Browlen)
      {
        Ccolinds(CrowStart + ci) = Bcolinds(BrowStart + bi);
        Cvalues(CrowStart + ci) = beta * Bvalues(BrowStart + bi);
        bi++;
        ci++;
      }
    }
    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    const AcolindsT Acolinds;
    const BcolindsT Bcolinds;
    CcolindsT Ccolinds;
    const AvaluesT Avalues;
    const BvaluesT Bvalues;
    CvaluesT Cvalues;
    const AscalarT alpha;
    const BscalarT beta;
  };

  template<typename size_type,
           typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
           typename AcolindsT, typename BcolindsT, typename CcolindsT,
           typename AvaluesT, typename BvaluesT, typename CvaluesT,
           typename AscalarT, typename BscalarT>
  struct UnsortedNumericSumFunctor
  {
    UnsortedNumericSumFunctor(const ArowptrsT Arowptrs_, const BrowptrsT Browptrs_, const CrowptrsT Crowptrs_,
    const AcolindsT Acolinds_, const BcolindsT Bcolinds_, CcolindsT Ccolinds_,
    const AvaluesT Avalues_, const BvaluesT Bvalues_, CvaluesT Cvalues_,
    const AscalarT alpha_, const BscalarT beta_,
    const CcolindsT Apos_, const CcolindsT Bpos_) :
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowptrs(Crowptrs_),
      Acolinds(Acolinds_),
      Bcolinds(Bcolinds_),
      Ccolinds(Ccolinds_),
      Avalues(Avalues_),
      Bvalues(Bvalues_),
      Cvalues(Cvalues_),
      alpha(alpha_),
      beta(beta_),
      Apos(Apos_),
      Bpos(Bpos_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const
    {
      size_type CrowStart = Crowptrs(i);
      size_type ArowStart = Arowptrs(i);
      size_type ArowEnd = Arowptrs(i + 1);
      size_type BrowStart = Browptrs(i);
      size_type BrowEnd = Browptrs(i + 1);
      //add in A entries, while setting C colinds
      for(size_type j = ArowStart; j < ArowEnd; j++)
      {
        Cvalues(CrowStart + Apos(j)) += alpha * Avalues(j);
        Ccolinds(CrowStart + Apos(j)) = Acolinds(j);
      }
      //add in B entries, while setting C colinds
      for(size_type j = BrowStart; j < BrowEnd; j++)
      {
        Cvalues(CrowStart + Bpos(j)) += beta * Bvalues(j);
        Ccolinds(CrowStart + Bpos(j)) = Bcolinds(j);
      }
    }
    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    const AcolindsT Acolinds;
    const BcolindsT Bcolinds;
    CcolindsT Ccolinds;
    const AvaluesT Avalues;
    const BvaluesT Bvalues;
    CvaluesT Cvalues;
    const AscalarT alpha;
    const BscalarT beta;
    const CcolindsT Apos;
    const CcolindsT Bpos;
  };

  template <typename KernelHandle,
            typename alno_row_view_t_,
            typename alno_nnz_view_t_,
            typename ascalar_t_,
            typename ascalar_nnz_view_t_,
            typename blno_row_view_t_,
            typename blno_nnz_view_t_,
            typename bscalar_t_,
            typename bscalar_nnz_view_t_,
            typename clno_row_view_t_,
            typename clno_nnz_view_t_,
            typename cscalar_nnz_view_t_>
  void spadd_numeric(
      KernelHandle* kernel_handle,
      const alno_row_view_t_ a_rowmap,
      const alno_nnz_view_t_ a_entries,
      const ascalar_nnz_view_t_ a_values,
      const ascalar_t_ alpha,
      const blno_row_view_t_ b_rowmap,
      const blno_nnz_view_t_ b_entries,
      const bscalar_nnz_view_t_ b_values,
      const bscalar_t_ beta,
      const clno_row_view_t_ c_rowmap,
      clno_nnz_view_t_ c_entries,
      cscalar_nnz_view_t_ c_values)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;
    typedef typename KernelHandle::nnz_scalar_t scalar_type;
    typedef typename KernelHandle::SPADDHandleType::execution_space execution_space;
    //Check that A/B/C data types match KernelHandle types, and that C data types are nonconst (doesn't matter if A/B types are const)
    static_assert(SAME_TYPE(ascalar_t_, scalar_type), "A scalar type must match handle scalar type");
    static_assert(SAME_TYPE(bscalar_t_, scalar_type), "B scalar type must match handle scalar type");
    static_assert(SAME_TYPE(typename alno_row_view_t_::value_type, size_type),
        "add_symbolic: A size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(SAME_TYPE(typename blno_row_view_t_::value_type, size_type),
        "add_symbolic: B size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(SAME_TYPE(typename clno_row_view_t_::value_type, size_type),
        "add_symbolic: C size_type must match KernelHandle size_type)");
    static_assert(std::is_same<typename clno_row_view_t_::value_type, typename clno_row_view_t_::value_type>::value,
        "add_symbolic: C size_type must not be const");
    static_assert(SAME_TYPE(typename alno_nnz_view_t_::value_type, ordinal_type),
        "add_symbolic: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename blno_nnz_view_t_::value_type, ordinal_type),
        "add_symbolic: B entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename clno_nnz_view_t_::value_type, ordinal_type),
        "add_symbolic: C entry type must match KernelHandle entry type (aka nnz_lno_t)");
    static_assert(std::is_same<typename clno_row_view_t_::non_const_value_type, typename clno_row_view_t_::value_type>::value,
        "add_symbolic: C entry type must not be const");
    static_assert(SAME_TYPE(typename ascalar_nnz_view_t_::value_type, scalar_type),
        "add_symbolic: A scalar type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename bscalar_nnz_view_t_::value_type, scalar_type),
        "add_symbolic: B scalar type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(SAME_TYPE(typename cscalar_nnz_view_t_::value_type, scalar_type),
        "add_symbolic: C scalar type must match KernelHandle entry type (aka nnz_lno_t)");
    static_assert(std::is_same<typename cscalar_nnz_view_t_::non_const_value_type, typename cscalar_nnz_view_t_::value_type>::value,
        "add_symbolic: C scalar type must not be const");
    typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
    auto addHandle = kernel_handle->get_spadd_handle();
    auto nrows = a_rowmap.dimension_0() - 1;
    if(addHandle->is_input_sorted())
    {
      SortedNumericSumFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                                           alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_,
                                           ascalar_nnz_view_t_, bscalar_nnz_view_t_, cscalar_nnz_view_t_,
                                           ascalar_t_, bscalar_t_>
        sortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values, alpha, beta);
      Kokkos::parallel_for(range_type(0, nrows), sortedNumeric);
      execution_space::fence();
    }
    else
    {
      //use a_pos and b_pos (set in the handle by symbolic) to quickly compute C entries and values
      UnsortedNumericSumFunctor<size_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                                           alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_,
                                           ascalar_nnz_view_t_, bscalar_nnz_view_t_, cscalar_nnz_view_t_,
                                           ascalar_t_, bscalar_t_>
        unsortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values, alpha, beta, addHandle->get_a_pos(), addHandle->get_b_pos());
      Kokkos::parallel_for(range_type(0, nrows), unsortedNumeric);
      execution_space::fence();
    }
    addHandle->set_call_numeric();
  }
}
}

#undef SAME_TYPE

#endif

