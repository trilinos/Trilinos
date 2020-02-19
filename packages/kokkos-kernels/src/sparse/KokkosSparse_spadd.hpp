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
#include "KokkosKernels_Sorting.hpp"

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
        ordinal_type nrows_,
        const typename ARowPtrsT::const_type& Arowptrs_,
        const AColIndsT& Acolinds_,
        const typename BRowPtrsT::const_type& Browptrs_,
        const BColIndsT& Bcolinds_,
        const CRowPtrsT& Crowcounts_) :
      nrows(nrows_),
      Arowptrs(Arowptrs_), Acolinds(Acolinds_),
      Browptrs(Browptrs_), Bcolinds(Bcolinds_),
      Crowcounts(Crowcounts_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
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
      if(i == nrows - 1)
      {
        //last workitem also zeros the one-past-end entry of row counts, so
        //that prefix sum is correct
        Crowcounts(nrows) = 0;
      }
    }
    ordinal_type nrows;
    const typename ARowPtrsT::const_type Arowptrs;
    const AColIndsT Acolinds;
    const typename BRowPtrsT::const_type Browptrs;
    const BColIndsT Bcolinds;
    CRowPtrsT Crowcounts;
  };

  //get upper bound for C entries per row (assumes worst case, that entries in A and B on each row are disjoint)
  template<typename size_type, typename ordinal_type, typename ARowPtrsT, typename BRowPtrsT, typename CRowPtrsT>
  struct UnsortedEntriesUpperBound
  {
      UnsortedEntriesUpperBound(
        ordinal_type nrows_,
        const typename ARowPtrsT::const_type& Arowptrs_,
        const typename BRowPtrsT::const_type& Browptrs_,
        const CRowPtrsT& Crowcounts_) :
      nrows(nrows_),
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowcounts(Crowcounts_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
    {
      Crowcounts(i) = (Arowptrs(i + 1) - Arowptrs(i)) + (Browptrs(i + 1) - Browptrs(i));
      if(i == nrows - 1)
      {
        //last workitem also zeros the one-past-end entry of row counts, so
        //that prefix sum is correct
        Crowcounts(nrows) = 0;
      }
    }
    ordinal_type nrows;
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
    UnmergedSumFunctor(ordinal_type nrows_, const ArowptrsT& Arowptrs_, const AcolindsT& Acolinds_,
                       const BrowptrsT& Browptrs_, const BcolindsT& Bcolinds_,
                       const CrowptrsT& Crowptrs_, const CcolindsT& Ccolinds_, const CcolindsT& ABperm_) :
      nrows(nrows_),
      Arowptrs(Arowptrs_), Acolinds(Acolinds_),
      Browptrs(Browptrs_), Bcolinds(Bcolinds_),
      Crowptrs(Crowptrs_), Ccolinds(Ccolinds_), ABperm(ABperm_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
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
    ordinal_type nrows;
    const ArowptrsT Arowptrs;
    const AcolindsT Acolinds;
    const BrowptrsT Browptrs;
    const BcolindsT Bcolinds;
    const CrowptrsT Crowptrs;
    CcolindsT Ccolinds;
    CcolindsT ABperm;
  };

  template<typename ExecSpace, typename size_type, typename ordinal_type, typename CrowptrsT, typename CcolindsT>
  struct SortEntriesFunctor
  {
    SortEntriesFunctor(const CrowptrsT& Crowptrs_, const CcolindsT& Ccolinds_, const CcolindsT& ABperm_) :
      Crowptrs(Crowptrs_),
      Ccolinds(Ccolinds_),
      CcolindsAux("C colind aux", Ccolinds_.extent(0)),
      ABperm(ABperm_),
      ABpermAux("AB perm aux", ABperm_.extent(0))
    {}
    typedef typename Kokkos::TeamPolicy<ExecSpace>::member_type TeamMember;
    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
    {
      //3: Sort each row's colinds (permuting values at same time), then count unique colinds (write that to Crowptr(i))
      //CrowptrTemp tells how many entries in each oversized row
      ordinal_type i = t.league_rank();
      size_type rowStart = Crowptrs(i);
      size_type rowEnd = Crowptrs(i + 1);
      size_type rowNum = rowEnd - rowStart;
      using lno_t = typename CcolindsT::non_const_value_type;
      using unsigned_lno_t = typename std::make_unsigned<lno_t>::type;
      KokkosKernels::Impl::SerialRadixSort2<size_type, unsigned_lno_t, lno_t>
        ((unsigned_lno_t*) Ccolinds.data() + rowStart, (unsigned_lno_t*) CcolindsAux.data() + rowStart,
         ABperm.data() + rowStart, ABpermAux.data() + rowStart, rowNum);
    }
    CrowptrsT Crowptrs;
    CcolindsT Ccolinds;
    CcolindsT CcolindsAux;
    CcolindsT ABperm;
    CcolindsT ABpermAux;
  };

  #ifdef KOKKOS_ENABLE_CUDA
  template<typename size_type, typename ordinal_type, typename CrowptrsT, typename CcolindsT>
  struct SortEntriesFunctor<Kokkos::Cuda, size_type, ordinal_type, CrowptrsT, CcolindsT>
  {
    SortEntriesFunctor(const CrowptrsT& Crowptrs_, CcolindsT& Ccolinds_, CcolindsT& ABperm_) :
      Crowptrs(Crowptrs_),
      Ccolinds(Ccolinds_),
      ABperm(ABperm_)
    {}
    typedef typename Kokkos::TeamPolicy<Kokkos::Cuda>::member_type TeamMember;
    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember t) const
    {
      //3: Sort each row's colinds (permuting values at same time), then count unique colinds (write that to Crowptr(i))
      //CrowptrTemp tells how many entries in each oversized row
      size_type i = t.league_rank();
      size_type rowStart = Crowptrs(i);
      size_type rowEnd = Crowptrs(i + 1);
      size_type rowNum = rowEnd - rowStart;
      KokkosKernels::Impl::TeamBitonicSort2<size_type, typename CcolindsT::non_const_value_type, typename CcolindsT::non_const_value_type, TeamMember>
        (Ccolinds.data() + rowStart, ABperm.data() + rowStart, rowNum, t);
    }
    CrowptrsT Crowptrs;
    CcolindsT Ccolinds;
    CcolindsT ABperm;
  };
  #endif

  template<typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT, typename CcolindsT>
  struct MergeEntriesFunctor
  {
    MergeEntriesFunctor(ordinal_type nrows_, const ArowptrsT& Arowptrs_, const BrowptrsT& Browptrs_, const CrowptrsT& Crowptrs_, const CrowptrsT& Crowcounts_,
        const CcolindsT& Ccolinds_, const CcolindsT& ABperm_, const CcolindsT& Apos_, const CcolindsT& Bpos_) :
      nrows(nrows_),
      Arowptrs(Arowptrs_),
      Browptrs(Browptrs_),
      Crowptrs(Crowptrs_),
      Crowcounts(Crowcounts_),
      Ccolinds(Ccolinds_),
      ABperm(ABperm_),
      Apos(Apos_),
      Bpos(Bpos_)
    {}
    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
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
      if(i == nrows - 1)
        Crowcounts(nrows) = 0;
    }
    ordinal_type nrows;
    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    CrowptrsT Crowcounts;
    CcolindsT Ccolinds;
    const CcolindsT ABperm;
    CcolindsT Apos;
    CcolindsT Bpos;
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
    if(a_rowmap.extent(0) == 0 || a_rowmap.extent(0) == 1)
    {
      //Have 0 rows, so nothing to do except set #nnz to 0
      addHandle->set_max_result_nnz(0);
      //If c_rowmap has a single entry, it must be 0
      if(c_rowmap.extent(0))
        Kokkos::deep_copy(c_rowmap, (size_type) 0);
      addHandle->set_call_symbolic();
      return;
    }
    ordinal_type nrows = a_rowmap.extent(0) - 1;
    typedef Kokkos::RangePolicy<execution_space, ordinal_type> range_type;
    using NoInitialize = Kokkos::ViewAllocateWithoutInitializing;
    if(addHandle->is_input_sorted())
    {
      //call entry count functor to get entry counts per row
      SortedCountEntries<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, alno_nnz_view_t_, blno_nnz_view_t_, clno_row_view_t_>
        countEntries(nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap);
      Kokkos::parallel_for("KokkosSparse::SpAdd::Symbolic::InputSorted::CountEntries", range_type(0, nrows), countEntries);
      KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<clno_row_view_t_, execution_space>(nrows + 1, c_rowmap);
    }
    else
    {
      //note: scoping individual parts of the process to free views sooner, minimizing peak memory usage
      //run the unsorted c_rowmap upper bound functor (just adds together A and B entry counts row by row)
      clno_row_view_t_ c_rowmap_upperbound(NoInitialize("C row counts upper bound"), nrows + 1);
      size_type c_nnz_upperbound = 0;
      {
        UnsortedEntriesUpperBound<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_>
          countEntries(nrows, a_rowmap, b_rowmap, c_rowmap_upperbound);
        Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::CountEntries", range_type(0, nrows), countEntries);
        KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<clno_row_view_t_, execution_space>(nrows + 1, c_rowmap_upperbound);
        Kokkos::deep_copy(c_nnz_upperbound, Kokkos::subview(c_rowmap_upperbound, nrows));
      }
      clno_nnz_view_t_ c_entries_uncompressed(NoInitialize("C entries uncompressed"), c_nnz_upperbound);
      clno_nnz_view_t_ ab_perm(NoInitialize("A and B permuted entry indices"), c_nnz_upperbound);
      //compute the unmerged sum
      UnmergedSumFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                         alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_> unmergedSum(
                         nrows, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap_upperbound, c_entries_uncompressed, ab_perm);
      Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::UnmergedSum", range_type(0, nrows), unmergedSum);
      //sort the unmerged sum
      SortEntriesFunctor<execution_space, size_type, ordinal_type, clno_row_view_t_, clno_nnz_view_t_>
        sortEntries(c_rowmap_upperbound, c_entries_uncompressed, ab_perm);
      Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::SortEntries",
          Kokkos::TeamPolicy<execution_space>(nrows, Kokkos::AUTO()), sortEntries);
      clno_nnz_view_t_ a_pos(NoInitialize("A entry positions"), a_entries.extent(0));
      clno_nnz_view_t_ b_pos(NoInitialize("B entry positions"), b_entries.extent(0));
      //merge the entries and compute Apos/Bpos, as well as Crowcounts
      {
        MergeEntriesFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_, clno_nnz_view_t_>
          mergeEntries(nrows, a_rowmap, b_rowmap, c_rowmap_upperbound, c_rowmap, c_entries_uncompressed, ab_perm, a_pos, b_pos);
        Kokkos::parallel_for("KokkosSparse::SpAdd:Symbolic::InputNotSorted::MergeEntries", range_type(0, nrows), mergeEntries);
        //compute actual c_rowmap
        KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<clno_row_view_t_, execution_space>(nrows + 1, c_rowmap);
      }
      addHandle->set_a_b_pos(a_pos, b_pos);
    }
    //provide the number of NNZ in C to user through handle
    size_type cmax;
    Kokkos::deep_copy(cmax, Kokkos::subview(c_rowmap, nrows));
    addHandle->set_max_result_nnz(cmax);
    addHandle->set_call_symbolic();
    addHandle->set_call_numeric(false);
    //this fence is for accurate timing from host
    execution_space().fence();
  }

  template<typename size_type, typename ordinal_type,
           typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
           typename AcolindsT, typename BcolindsT, typename CcolindsT,
           typename AvaluesT, typename BvaluesT, typename CvaluesT,
           typename AscalarT, typename BscalarT>
  struct SortedNumericSumFunctor
  {
    SortedNumericSumFunctor(const ArowptrsT& Arowptrs_, const BrowptrsT& Browptrs_, const CrowptrsT& Crowptrs_,
    const AcolindsT& Acolinds_, const BcolindsT& Bcolinds_, const CcolindsT& Ccolinds_,
    const AvaluesT& Avalues_, const BvaluesT& Bvalues_, const CvaluesT& Cvalues_,
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
    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
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

  template<typename size_type, typename ordinal_type,
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

    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
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
    //rowmap length can be 0 or 1 if #rows is 0.
    //Otherwise, it's always #rows+1.
    if(a_rowmap.extent(0) == 0 || a_rowmap.extent(0) == 1)
    {
      addHandle->set_call_numeric();
      return;
    }
    ordinal_type nrows = a_rowmap.extent(0) - 1;
    if(addHandle->is_input_sorted())
    {
      SortedNumericSumFunctor<size_type, ordinal_type, alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                                           alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_,
                                           ascalar_nnz_view_t_, bscalar_nnz_view_t_, cscalar_nnz_view_t_,
                                           ascalar_t_, bscalar_t_>
        sortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values, alpha, beta);
      Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputSorted", range_type(0, nrows), sortedNumeric);
    }
    else
    {
      //use a_pos and b_pos (set in the handle by symbolic) to quickly compute C entries and values
      UnsortedNumericSumFunctor<size_type, ordinal_type,
                                           alno_row_view_t_, blno_row_view_t_, clno_row_view_t_,
                                           alno_nnz_view_t_, blno_nnz_view_t_, clno_nnz_view_t_,
                                           ascalar_nnz_view_t_, bscalar_nnz_view_t_, cscalar_nnz_view_t_,
                                           ascalar_t_, bscalar_t_>
        unsortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values, alpha, beta, addHandle->get_a_pos(), addHandle->get_b_pos());
      Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputNotSorted", range_type(0, nrows), unsortedNumeric);
    }
    addHandle->set_call_numeric();
    //this fence is for accurate timing from host
    execution_space().fence();
  }
}
}

#undef SAME_TYPE

#endif

