// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSGRAPH_DEF_HPP
#define TPETRA_CRSGRAPH_DEF_HPP

/// \file Tpetra_CrsGraph_def.hpp
/// \brief Definition of the Tpetra::CrsGraph class
///
/// If you want to use Tpetra::CrsGraph, include "Tpetra_CrsGraph.hpp"
/// (a file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::CrsGraph, include
/// "Tpetra_CrsGraph_decl.hpp".

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include <algorithm>
#include <string>
#include <utility>

namespace Tpetra {

  namespace Details {

    // Implementation of Tpetra::CrsGraph::getLocalDiagOffsets, for
    // the fillComplete case.
    //
    // FIXME (mfh 12 Mar 2016) There's currently no way to make a
    // MemoryUnmanaged Kokkos::StaticCrsGraph.  Thus, we have to do
    // this separately for its column indices.  We want the column
    // indices to be unmanaged because we need to take subviews in
    // this kernel.  Taking a subview of a managed View updates the
    // reference count, which is a thread scalability bottleneck.
    template<class LO, class GO, class Node>
    class GetLocalDiagOffsets {
    public:
      typedef typename Node::device_type device_type;
      // mfh 12 Mar 2016: getLocalDiagOffsets returns offsets as
      // size_t.  However, see Github Issue #213.
      typedef size_t diag_offset_type;
      typedef Kokkos::View<diag_offset_type*, device_type,
                           Kokkos::MemoryUnmanaged> diag_offsets_type;
      typedef typename ::Tpetra::CrsGraph<LO, GO, Node> global_graph_type;
      typedef typename global_graph_type::local_graph_type local_graph_type;
      typedef typename global_graph_type::map_type::local_map_type local_map_type;
      typedef Kokkos::View<const typename local_graph_type::size_type*,
                           Kokkos::LayoutLeft, device_type> row_offsets_type;
      // This is unmanaged for performance, because we need to take
      // subviews inside the functor.
      typedef Kokkos::View<const LO*, Kokkos::LayoutLeft, device_type,
                           Kokkos::MemoryUnmanaged> lcl_col_inds_type;

      GetLocalDiagOffsets (const diag_offsets_type& diagOffsets,
                           const local_map_type& lclRowMap,
                           const local_map_type& lclColMap,
                           const row_offsets_type& ptr,
                           const lcl_col_inds_type& ind,
                           const bool isSorted) :
        diagOffsets_ (diagOffsets),
        lclRowMap_ (lclRowMap),
        lclColMap_ (lclColMap),
        ptr_ (ptr),
        ind_ (ind),
        isSorted_ (isSorted)
      {
        typedef typename device_type::execution_space execution_space;
        typedef Kokkos::RangePolicy<execution_space, LO> policy_type;

        const LO lclNumRows = lclRowMap.getNodeNumElements ();
        policy_type range (0, lclNumRows);
        Kokkos::parallel_for (range, *this);
      }

      KOKKOS_INLINE_FUNCTION void
      operator() (const LO& lclRowInd) const
      {
        const size_t STINV =
          Tpetra::Details::OrdinalTraits<diag_offset_type>::invalid ();
        const GO gblRowInd = lclRowMap_.getGlobalElement (lclRowInd);
        const GO gblColInd = gblRowInd;
        const LO lclColInd = lclColMap_.getLocalElement (gblColInd);

        if (lclColInd == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          diagOffsets_[lclRowInd] = STINV;
        }
        else {
          // Could be empty, but that's OK.
          const LO numEnt = ptr_[lclRowInd+1] - ptr_[lclRowInd];
          // std::pair doesn't have its methods marked as device
          // functions, so we have to use Kokkos::pair.
          auto lclColInds =
            Kokkos::subview (ind_, Kokkos::make_pair (ptr_[lclRowInd],
                                                      ptr_[lclRowInd+1]));
          using ::KokkosSparse::findRelOffset;
          const LO diagOffset =
            findRelOffset<LO, lcl_col_inds_type> (lclColInds, numEnt,
                                                  lclColInd, 0, isSorted_);
          diagOffsets_[lclRowInd] = (diagOffset == numEnt) ? STINV :
            static_cast<diag_offset_type> (diagOffset);
        }
      }

    private:
      diag_offsets_type diagOffsets_;
      local_map_type lclRowMap_;
      local_map_type lclColMap_;
      row_offsets_type ptr_;
      lcl_col_inds_type ind_;
      bool isSorted_;
    };

  } // namespace Details

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            size_t maxNumEntriesPerRow,
            ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , numAllocForAllRows_ (maxNumEntriesPerRow)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,maxNumEntriesPerRow,"
      "pftype,params): ";
    staticAssertions ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      maxNumEntriesPerRow == Teuchos::OrdinalTraits<size_t>::invalid (),
      std::invalid_argument, "The allocation hint maxNumEntriesPerRow must be "
      "a valid size_t value, which in this case means it must not be "
      "Teuchos::OrdinalTraits<size_t>::invalid().");
    resumeFill (params);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const size_t maxNumEntriesPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , numAllocForAllRows_ (maxNumEntriesPerRow)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,colMap,maxNumEntriesPerRow,"
      "pftype,params): ";
    staticAssertions ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      maxNumEntriesPerRow == Teuchos::OrdinalTraits<size_t>::invalid (),
      std::invalid_argument, "The allocation hint maxNumEntriesPerRow must be "
      "a valid size_t value, which in this case means it must not be "
      "Teuchos::OrdinalTraits<size_t>::invalid().");
    resumeFill (params);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , numAllocForAllRows_ (0)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.size ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " << numEntPerRow.size ()
      << " != the local number of rows " << lclNumRows << " as specified by "
      "the input row Map.");

#ifdef HAVE_TPETRA_DEBUG
    for (size_t r = 0; r < lclNumRows; ++r) {
      const size_t curRowCount = numEntPerRow[r];
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
        std::invalid_argument, "numEntPerRow(" << r << ") specifies an invalid "
        "number of entries (Teuchos::OrdinalTraits<size_t>::invalid()).");
    }
#endif // HAVE_TPETRA_DEBUG

    // Deep-copy the input (ArrayRCP, therefore host accessible) into
    // k_numAllocPerRow_.  The latter is a const View, so we have to
    // copy into a nonconst View first, then assign.
    typedef decltype (k_numAllocPerRow_) out_view_type;
    typedef typename out_view_type::non_const_type nc_view_type;
    typedef Kokkos::View<const size_t*,
                         typename nc_view_type::array_layout,
                         Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    nc_view_type numAllocPerRowOut ("Tpetra::CrsGraph::numAllocPerRow",
                                    lclNumRows);
    Kokkos::deep_copy (numAllocPerRowOut, numAllocPerRowIn);
    k_numAllocPerRow_ = numAllocPerRowOut;

    resumeFill (params);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , k_numAllocPerRow_ (numEntPerRow.h_view)
    , numAllocForAllRows_ (0)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.dimension_0 ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " <<
      numEntPerRow.dimension_0 () << " != the local number of rows " <<
      lclNumRows << " as specified by " "the input row Map.");

#ifdef HAVE_TPETRA_DEBUG
    for (size_t r = 0; r < lclNumRows; ++r) {
      const size_t curRowCount = numEntPerRow.h_view(r);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
        std::invalid_argument, "numEntPerRow(" << r << ") specifies an invalid "
        "number of entries (Teuchos::OrdinalTraits<size_t>::invalid()).");
    }
#endif // HAVE_TPETRA_DEBUG

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , k_numAllocPerRow_ (numEntPerRow.h_view)
    , numAllocForAllRows_ (0)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,colMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.dimension_0 ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " <<
      numEntPerRow.dimension_0 () << " != the local number of rows " <<
      lclNumRows << " as specified by " "the input row Map.");

#ifdef HAVE_TPETRA_DEBUG
    for (size_t r = 0; r < lclNumRows; ++r) {
      const size_t curRowCount = numEntPerRow.h_view(r);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
        std::invalid_argument, "numEntPerRow(" << r << ") specifies an invalid "
        "number of entries (Teuchos::OrdinalTraits<size_t>::invalid()).");
    }
#endif // HAVE_TPETRA_DEBUG

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
            ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , numAllocForAllRows_ (0)
    , storageStatus_ (pftype == StaticProfile ?
                      Details::STORAGE_1D_UNPACKED :
                      Details::STORAGE_2D)
    , indicesAreAllocated_ (false)
    , indicesAreLocal_ (false)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    const char tfecfFuncName[] = "CrsGraph(rowMap,colMap,numEntPerRow,pftype,"
      "params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.size ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " << numEntPerRow.size ()
      << " != the local number of rows " << lclNumRows << " as specified by "
      "the input row Map.");

#ifdef HAVE_TPETRA_DEBUG
    for (size_t r = 0; r < lclNumRows; ++r) {
      const size_t curRowCount = numEntPerRow[r];
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
        std::invalid_argument, "numEntPerRow(" << r << ") specifies an invalid "
        "number of entries (Teuchos::OrdinalTraits<size_t>::invalid()).");
    }
#endif // HAVE_TPETRA_DEBUG

    // Deep-copy the input (ArrayRCP, therefore host accessible) into
    // k_numAllocPerRow_.  The latter is a const View, so we have to
    // copy into a nonconst View first, then assign.
    typedef decltype (k_numAllocPerRow_) out_view_type;
    typedef typename out_view_type::non_const_type nc_view_type;
    typedef Kokkos::View<const size_t*,
                         typename nc_view_type::array_layout,
                         Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    nc_view_type numAllocPerRowOut ("Tpetra::CrsGraph::numAllocPerRow",
                                    lclNumRows);
    Kokkos::deep_copy (numAllocPerRowOut, numAllocPerRowIn);
    k_numAllocPerRow_ = numAllocPerRowOut;

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const typename local_graph_type::row_map_type& rowPointers,
            const typename local_graph_type::entries_type::non_const_type& columnIndices,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_(rowMap)
    , colMap_(colMap)
    , globalNumEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalNumDiags_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalMaxNumRowEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , nodeNumEntries_(0)
    , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
    , pftype_(StaticProfile)
    , numAllocForAllRows_(0)
    , storageStatus_ (Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_(true)
    , indicesAreLocal_(true)
    , indicesAreGlobal_(false)
    , fillComplete_(false)
    , indicesAreSorted_(true)
    , noRedundancies_(true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_(true)
  {
    staticAssertions ();
    setAllIndices (rowPointers, columnIndices);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::ArrayRCP<size_t>& rowPointers,
            const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , globalNumEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalNumDiags_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalMaxNumRowEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (StaticProfile)
    , numAllocForAllRows_ (0)
    , storageStatus_ (Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_ (true)
    , indicesAreLocal_ (true)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_ (true)
  {
    staticAssertions ();
    setAllIndices (rowPointers, columnIndices);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const local_graph_type& k_local_graph_,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
    : DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, node_type> (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , lclGraph_ (k_local_graph_)
    , globalNumEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalNumDiags_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalMaxNumRowEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , nodeNumEntries_ (0) // FIXME (mfh 17 Mar 2014) should get from lclGraph_ right now
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (StaticProfile)
    , numAllocForAllRows_ (0)
    , storageStatus_ (Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_ (true)
    , indicesAreLocal_ (true)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , sortGhostsAssociatedWithEachProcessor_(true)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::rcp;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;

    staticAssertions();
    const char tfecfFuncName[] = "CrsGraph(Map,Map,Kokkos::LocalStaticCrsGraph)";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      colMap.is_null (), std::runtime_error,
      ": The input column Map must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_local_graph_.numRows () != rowMap->getNodeNumElements (),
      std::runtime_error,
      ": The input row Map and the input local graph need to have the same "
      "number of rows.  The row Map claims " << rowMap->getNodeNumElements ()
      << " row(s), but the local graph claims " << k_local_graph_.numRows ()
      << " row(s).");
    // NOTE (mfh 17 Mar 2014) getNodeNumRows() returns
    // rowMap_->getNodeNumElements(), but it doesn't have to.
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   k_local_graph_.numRows () != getNodeNumRows (), std::runtime_error,
    //   ": The input row Map and the input local graph need to have the same "
    //   "number of rows.  The row Map claims " << getNodeNumRows () << " row(s), "
    //   "but the local graph claims " << k_local_graph_.numRows () << " row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_lclInds1D_.dimension_0 () != 0 || k_gblInds1D_.dimension_0 () != 0, std::logic_error,
      ": cannot have 1D data structures allocated.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! lclInds2D_.is_null () || ! gblInds2D_.is_null (), std::logic_error,
      ": cannot have 2D data structures allocated.");

    nodeNumAllocated_ = k_local_graph_.row_map (getNodeNumRows ());
    nodeNumEntries_ = k_local_graph_.row_map (getNodeNumRows ());

    // NOTE (mfh 17 Mar 2014) We also need a version of this CrsGraph
    // constructor that takes a domain and range Map, as well as a row
    // and column Map.  In that case, we must pass the domain and
    // range Map into the following method.
    setDomainRangeMaps (rowMap_, rowMap_);
    makeImportExport ();

    k_lclInds1D_ = lclGraph_.entries;
    k_rowPtrs_ = lclGraph_.row_map;

    typename local_graph_type::row_map_type d_ptrs = lclGraph_.row_map;
    typename local_graph_type::entries_type d_inds = lclGraph_.entries;

    // Reset local properties
    upperTriangular_ = true;
    lowerTriangular_ = true;
    nodeMaxNumRowEntries_ = 0;
    nodeNumDiags_         = 0;

    // Compute triangular properties
    const size_t numLocalRows = getNodeNumRows ();
    for (size_t localRow = 0; localRow < numLocalRows; ++localRow) {
      const GO globalRow = rowMap_->getGlobalElement (localRow);
      const LO rlcid = colMap_->getLocalElement (globalRow);

      // It's entirely possible that the local matrix has no entries
      // in the column corresponding to the current row.  In that
      // case, the column Map may not necessarily contain that GID.
      // This is why we check whether rlcid is "invalid" (which means
      // that globalRow is not a GID in the column Map).
      if (rlcid != Teuchos::OrdinalTraits<LO>::invalid ()) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          rlcid + 1 >= static_cast<LO> (d_ptrs.dimension_0 ()),
          std::runtime_error, ": The given row Map and/or column Map is/are "
          "not compatible with the provided local graph.");
        if (d_ptrs(rlcid) != d_ptrs(rlcid + 1)) {
          const size_t smallestCol =
            static_cast<size_t> (d_inds(d_ptrs(rlcid)));
          const size_t largestCol =
            static_cast<size_t> (d_inds(d_ptrs(rlcid + 1)-1));
          if (smallestCol < localRow) {
            upperTriangular_ = false;
          }
          if (localRow < largestCol) {
            lowerTriangular_ = false;
          }
          for (size_t i = d_ptrs(rlcid); i < d_ptrs(rlcid + 1); ++i) {
            if (d_inds(i) == rlcid) {
              ++nodeNumDiags_;
            }
          }
        }
        nodeMaxNumRowEntries_ =
          std::max (static_cast<size_t> (d_ptrs(rlcid + 1) - d_ptrs(rlcid)),
                    nodeMaxNumRowEntries_);
      }
    }

    haveLocalConstants_ = true;
    computeGlobalConstants ();

    fillComplete_ = true;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  ~CrsGraph ()
  {}


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Teuchos::ParameterList>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getValidParameters () const
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;

    RCP<ParameterList> params = parameterList ("Tpetra::CrsGraph");

    // Make a sublist for the Import.
    RCP<ParameterList> importSublist = parameterList ("Import");

    // FIXME (mfh 02 Apr 2012) We should really have the Import and
    // Export objects fill in these lists.  However, we don't want to
    // create an Import or Export unless we need them.  For now, we
    // know that the Import and Export just pass the list directly to
    // their Distributor, so we can create a Distributor here
    // (Distributor's constructor is a lightweight operation) and have
    // it fill in the list.

    // Fill in Distributor default parameters by creating a
    // Distributor and asking it to do the work.
    Distributor distributor (rowMap_->getComm (), importSublist);
    params->set ("Import", *importSublist, "How the Import performs communication.");

    // Make a sublist for the Export.  For now, it's a clone of the
    // Import sublist.  It's not a shallow copy, though, since we
    // might like the Import to do communication differently than the
    // Export.
    params->set ("Export", *importSublist, "How the Export performs communication.");

    return params;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<const Teuchos::ParameterList> validParams =
      getValidParameters ();
    params->validateParametersAndSetDefaults (*validParams);
    this->setMyParamList (params);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalNumRows () const
  {
    return rowMap_->getGlobalNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalNumCols () const
  {
    const char tfecfFuncName[] = "getGlobalNumCols: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete () || getDomainMap ().is_null (), std::runtime_error,
      "The graph does not have a domain Map.  You may not call this method in "
      "that case.");
    return getDomainMap ()->getGlobalNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeNumRows () const
  {
    return rowMap_.is_null () ? static_cast<size_t> (0) :
      rowMap_->getNodeNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeNumCols () const
  {
    const char tfecfFuncName[] = "getNodeNumCols: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error,
      "The graph does not have a column Map.  You may not call this method "
      "unless the graph has a column Map.  This requires either that a custom "
      "column Map was given to the constructor, or that fillComplete() has "
      "been called.");
    return colMap_.is_null () ? static_cast<size_t> (0) :
      colMap_->getNodeNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeNumDiags () const
  {
    return nodeNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalNumDiags () const
  {
    return globalNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNode () const
  {
    return rowMap_.is_null () ? Teuchos::null : rowMap_->getNode ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getRowMap () const
  {
    return rowMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getColMap () const
  {
    return colMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getDomainMap () const
  {
    return domainMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getRangeMap () const
  {
    return rangeMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getImporter () const
  {
    return importer_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getExporter () const
  {
    return exporter_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  hasColMap () const
  {
    return ! colMap_.is_null ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isStorageOptimized () const
  {
    // FIXME (mfh 07 Aug 2014) Why wouldn't storage be optimized if
    // getNodeNumRows() is zero?

    const bool isOpt = indicesAreAllocated_ &&
      k_numRowEntries_.dimension_0 () == 0 &&
      getNodeNumRows () > 0;

#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "isStorageOptimized";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isOpt && getProfileType () == DynamicProfile, std::logic_error,
      ": The matrix claims to have optimized storage, but getProfileType() "
      "returns DynamicProfile.  This should never happen.  Please report this "
      "bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    return isOpt;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  ProfileType
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getProfileType () const
  {
    return pftype_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalNumEntries () const
  {
    return globalNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeNumEntries () const
  {
    return nodeNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalMaxNumRowEntries () const
  {
    return globalMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeMaxNumRowEntries () const
  {
    return nodeMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isFillComplete () const
  {
    return fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isFillActive () const
  {
    return ! fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isUpperTriangular () const
  {
    return upperTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isLowerTriangular () const
  {
    return lowerTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isLocallyIndexed () const
  {
    return indicesAreLocal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isGloballyIndexed () const
  {
    return indicesAreGlobal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeAllocationSize () const
  {
    return nodeNumAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getComm () const
  {
    return rowMap_.is_null () ? Teuchos::null : rowMap_->getComm ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  GlobalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getIndexBase () const
  {
    return rowMap_->getIndexBase ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  indicesAreAllocated () const
  {
    return indicesAreAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isSorted () const
  {
    return indicesAreSorted_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  isMerged () const
  {
    return noRedundancies_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  setLocallyModified ()
  {
    // FIXME (mfh 07 May 2013) How do we know that the change
    // introduced a redundancy, or even that it invalidated the sorted
    // order of indices?  CrsGraph has always made this conservative
    // guess.  It could be a bit costly to check at insertion time,
    // though.
    indicesAreSorted_ = false;
    noRedundancies_ = false;

    // We've modified the graph, so we'll have to recompute local
    // constants like the number of diagonal entries on this process.
    haveLocalConstants_ = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  allocateIndices (const ELocalGlobal lg)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ArrayRCP<size_t>::size_type size_type;
    typedef typename local_graph_type::row_map_type::non_const_type
      non_const_row_map_type;
    typedef typename local_graph_type::entries_type::non_const_type
      lcl_col_inds_type;
    typedef Kokkos::View<GlobalOrdinal*,
      typename lcl_col_inds_type::array_layout,
      device_type> gbl_col_inds_type;
    const char tfecfFuncName[] = "allocateIndices: ";

    // This is a protected function, only callable by us.  If it was
    // called incorrectly, it is our fault.  That's why the tests
    // below throw std::logic_error instead of std::invalid_argument.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed () && lg == GlobalIndices, std::logic_error,
      "The graph is locally indexed, but Tpetra code is calling this method "
      "with lg=GlobalIndices.  Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed () && lg == LocalIndices, std::logic_error,
      "The graph is globally indexed, but Tpetra code is calling this method "
      "with lg=LocalIndices.  Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated (), std::logic_error, "The graph's indices are "
      "already allocated, but Tpetra code is calling allocateIndices) again.  "
      "Please report this bug to the Tpetra developers.");

    const size_t numRows = getNodeNumRows ();

    if (getProfileType () == StaticProfile) {
      //
      //  STATIC ALLOCATION PROFILE
      //
      non_const_row_map_type k_rowPtrs ("Tpetra::CrsGraph::ptr", numRows + 1);

      if (k_numAllocPerRow_.dimension_0 () != 0) {
        // It's OK to throw std::invalid_argument here, because we
        // haven't incurred any side effects yet.  Throwing that
        // exception (and not, say, std::logic_error) implies that the
        // instance can recover.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          k_numAllocPerRow_.dimension_0 () != numRows, std::invalid_argument,
          "k_numAllocPerRow_ is allocated (has length != 0), but its length = "
          << k_numAllocPerRow_.dimension_0 () << " != numRows = " << numRows
          << ".");

        // k_numAllocPerRow_ is a host View, but k_rowPtrs (the thing
        // we want to compute here) lives on device.  That's OK;
        // computeOffsetsFromCounts can handle this case.
        using ::Tpetra::Details::computeOffsetsFromCounts;

        // FIXME (mfh 27 Jun 2016) Currently, computeOffsetsFromCounts
        // doesn't attempt to check its input for "invalid" flag
        // values.  For now, we omit that feature of the sequential
        // code disabled below.
        computeOffsetsFromCounts (k_rowPtrs, k_numAllocPerRow_);
      }
      else {
        // It's OK to throw std::invalid_argument here, because we
        // haven't incurred any side effects yet.  Throwing that
        // exception (and not, say, std::logic_error) implies that the
        // instance can recover.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numAllocForAllRows_ == Teuchos::OrdinalTraits<size_t>::invalid (),
           std::invalid_argument, "numAllocForAllRows_ has an invalid value, "
           "namely Teuchos::OrdinalTraits<size_t>::invalid() = " <<
           Teuchos::OrdinalTraits<size_t>::invalid () << ".");

        using ::Tpetra::Details::computeOffsetsFromConstantCount;
        computeOffsetsFromConstantCount (k_rowPtrs, numAllocForAllRows_);
      }

      // "Commit" the resulting row offsets.
      k_rowPtrs_ = k_rowPtrs;

      // FIXME (mfh 05,11 Aug 2014) This assumes UVM, since k_rowPtrs_
      // is currently a device View.  Should instead use a DualView.
      const size_type numInds = static_cast<size_type> (k_rowPtrs_(numRows));
      if (lg == LocalIndices) {
        k_lclInds1D_ = lcl_col_inds_type ("Tpetra::CrsGraph::ind", numInds);
      }
      else {
        k_gblInds1D_ = gbl_col_inds_type ("Tpetra::CrsGraph::ind", numInds);
      }
      nodeNumAllocated_ = numInds;
      storageStatus_ = Details::STORAGE_1D_UNPACKED;
    }
    else {
      //
      //  DYNAMIC ALLOCATION PROFILE
      //

      const bool useNumAllocPerRow = (k_numAllocPerRow_.dimension_0 () != 0);

      if (lg == LocalIndices) {
        lclInds2D_ = arcp<Array<LocalOrdinal> > (numRows);
        nodeNumAllocated_ = 0;
        for (size_t i = 0; i < numRows; ++i) {
          const size_t howMany = useNumAllocPerRow ?
            k_numAllocPerRow_(i) : numAllocForAllRows_;
          nodeNumAllocated_ += howMany;
          if (howMany > 0) {
            lclInds2D_[i].resize (howMany);
          }
        }
      }
      else { // allocate global indices
        gblInds2D_ = arcp<Array<GlobalOrdinal> > (numRows);
        nodeNumAllocated_ = 0;
        for (size_t i = 0; i < numRows; ++i) {
          const size_t howMany = useNumAllocPerRow ?
            k_numAllocPerRow_(i) : numAllocForAllRows_;
          nodeNumAllocated_ += howMany;
          if (howMany > 0) {
            gblInds2D_[i].resize (howMany);
          }
        }
      }
      storageStatus_ = Details::STORAGE_2D;
    }

    indicesAreLocal_  = (lg == LocalIndices);
    indicesAreGlobal_ = (lg == GlobalIndices);

    if (numRows > 0) { // reallocate k_numRowEntries_ & fill w/ 0s
      using Kokkos::ViewAllocateWithoutInitializing;
      typedef decltype (k_numRowEntries_) row_ent_type;
      const char label[] = "Tpetra::CrsGraph::numRowEntries";

      row_ent_type numRowEnt (ViewAllocateWithoutInitializing (label), numRows);
      Kokkos::deep_copy (numRowEnt, static_cast<size_t> (0)); // fill w/ 0s
      k_numRowEntries_ = numRowEnt; // "commit" our allocation
    }

    // done with these
    numAllocForAllRows_ = 0;
    k_numAllocPerRow_ = decltype (k_numAllocPerRow_) ();
    indicesAreAllocated_ = true;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayView<const LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalView (const RowInfo rowinfo) const
  {
    using Kokkos::subview;
    typedef LocalOrdinal LO;
    typedef Kokkos::View<const LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return Teuchos::ArrayView<const LO> ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        row_view_type rowView = subview (row_view_type (k_lclInds1D_), rng);
        const LO* const rowViewRaw = (len == 0) ? NULL : rowView.ptr_on_device ();
        return Teuchos::ArrayView<const LO> (rowViewRaw, len, Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        return lclInds2D_[rowinfo.localRow] ();
      }
      else {
        return Teuchos::ArrayView<const LO> (); // nothing in the row to view
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  LocalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalViewRawConst (const LocalOrdinal*& lclInds,
                        LocalOrdinal& numEnt,
                        const RowInfo& rowinfo) const
  {
    lclInds = NULL;
    numEnt = 0;

    if (rowinfo.allocSize != 0) {
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
#ifdef HAVE_TPETRA_DEBUG
        if (rowinfo.offset1D + rowinfo.allocSize >
            static_cast<size_t> (k_lclInds1D_.dimension_0 ())) {
          return static_cast<LocalOrdinal> (-1);
        }
#endif // HAVE_TPETRA_DEBUG
        lclInds = &k_lclInds1D_[rowinfo.offset1D];
        numEnt = rowinfo.allocSize;
      }
      else { // 2-D storage
#ifdef HAVE_TPETRA_DEBUG
        if (rowinfo.localRow >= static_cast<size_t> (lclInds2D_.size ())) {
          return static_cast<LocalOrdinal> (-1);
        }
#endif // HAVE_TPETRA_DEBUG
        // Use a const reference so we don't touch the ArrayRCP's ref
        // count, since ArrayRCP's ref count is not thread safe.
        const auto& curRow = lclInds2D_[rowinfo.localRow];
        if (! curRow.empty ()) {
          lclInds = curRow.getRawPtr ();
          numEnt = curRow.size ();
        }
      }
    }
    return static_cast<LocalOrdinal> (0);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayView<LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalViewNonConst (const RowInfo rowinfo)
  {
    using Kokkos::subview;
    typedef LocalOrdinal LO;
    typedef Kokkos::View<LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) { // nothing in the row to view
      return Teuchos::ArrayView<LO> ();
    }
    else {
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        row_view_type rowView = subview (row_view_type (k_lclInds1D_), rng);
        LO* const rowViewRaw = (len == 0) ? NULL : rowView.ptr_on_device ();
        return Teuchos::ArrayView<LO> (rowViewRaw, len, Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        return lclInds2D_[rowinfo.localRow] ();
      }
      else {
        return Teuchos::ArrayView<LO> (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Kokkos::View<const LocalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalKokkosRowView (const RowInfo& rowinfo) const
  {
    typedef LocalOrdinal LO;
    typedef Kokkos::View<const LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (k_lclInds1D_), rng);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        // Use a reference, to avoid doing the array lookup twice.
        //
        // NOTE (mfh 18 Jul 2016) Don't create a Teuchos::ArrayView,
        // because this is not thread safe in a debug build.
        Teuchos::Array<LocalOrdinal>& lclInds = lclInds2D_[rowinfo.localRow];
        const LO* rowPtr = lclInds.getRawPtr ();
        const auto rowSize = lclInds.size ();
        return row_view_type (rowPtr, rowSize);
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Kokkos::View<LocalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalKokkosRowViewNonConst (const RowInfo& rowinfo)
  {
    typedef LocalOrdinal LO;
    typedef Kokkos::View<LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (k_lclInds1D_), rng);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        Teuchos::ArrayView<LO> rowAv = lclInds2D_[rowinfo.localRow] ();
        return row_view_type (rowAv.getRawPtr (), rowAv.size ());
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Kokkos::View<const GlobalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalKokkosRowView (const RowInfo& rowinfo) const
  {
    typedef GlobalOrdinal GO;
    typedef Kokkos::View<const GO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (this->k_gblInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_gblInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (this->k_gblInds1D_), rng);
      }
      else if (! this->gblInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        Teuchos::ArrayView<const GO> rowAv = this->gblInds2D_[rowinfo.localRow] ();
        // FIXME (mfh 26 Nov 2015) This assumes UVM, because it
        // assumes that host code can access device memory through
        // Teuchos::ArrayView.
        return row_view_type (rowAv.getRawPtr (), rowAv.size ());
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayView<const GlobalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalView (const RowInfo rowinfo) const
  {
    Teuchos::ArrayView<const GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (k_gblInds1D_.dimension_0 () != 0) {
        auto rng = std::make_pair (rowinfo.offset1D,
                                   rowinfo.offset1D + rowinfo.allocSize);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_gblInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        Kokkos::View<const GlobalOrdinal*, execution_space,
          Kokkos::MemoryUnmanaged> k_gblInds1D_unmanaged = k_gblInds1D_;
        view = Kokkos::Compat::getConstArrayView (Kokkos::subview (k_gblInds1D_unmanaged, rng));
      }
      else if (! gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  LocalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalViewRawConst (const GlobalOrdinal*& gblInds,
                         LocalOrdinal& numEnt,
                         const RowInfo& rowinfo) const
  {
    gblInds = NULL;
    numEnt = 0;

    if (rowinfo.allocSize != 0) {
      if (k_gblInds1D_.dimension_0 () != 0) { // 1-D storage
#ifdef HAVE_TPETRA_DEBUG
        if (rowinfo.offset1D + rowinfo.allocSize >
            static_cast<size_t> (k_gblInds1D_.dimension_0 ())) {
          return static_cast<LocalOrdinal> (-1);
        }
#endif // HAVE_TPETRA_DEBUG
        gblInds = &k_gblInds1D_[rowinfo.offset1D];
        numEnt = rowinfo.allocSize;
      }
      else {
#ifdef HAVE_TPETRA_DEBUG
        if (rowinfo.localRow >= static_cast<size_t> (gblInds2D_.size ())) {
          return static_cast<LocalOrdinal> (-1);
        }
#endif // HAVE_TPETRA_DEBUG
        const auto& curRow = gblInds2D_[rowinfo.localRow];
        if (! curRow.empty ()) {
          gblInds = curRow.getRawPtr ();
          numEnt = curRow.size ();
        }
      }
    }
    return static_cast<LocalOrdinal> (0);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayView<GlobalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalViewNonConst (const RowInfo rowinfo)
  {
    Teuchos::ArrayView<GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (k_gblInds1D_.dimension_0 () != 0) {
        auto rng = std::make_pair (rowinfo.offset1D,
                                   rowinfo.offset1D + rowinfo.allocSize);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_gblInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        Kokkos::View<GlobalOrdinal*, execution_space,
          Kokkos::MemoryUnmanaged> k_gblInds1D_unmanaged = k_gblInds1D_;
        view = Kokkos::Compat::getArrayView (Kokkos::subview (k_gblInds1D_unmanaged, rng));
      }
      else if (! gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  RowInfo
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getRowInfo (const LocalOrdinal myRow) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getRowInfo: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::logic_error,
      "Late catch! Graph does not have row info anymore.  "
      "Error should have been caught earlier.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    if (! hasRowInfo () || rowMap_.is_null () ||
        ! rowMap_->isNodeLocalElement (myRow)) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }

    ret.localRow = static_cast<size_t> (myRow);
    if (nodeNumAllocated_ != 0 && nodeNumAllocated_ != STINV) {
      // graph data structures have the info that we need
      //
      // if static graph, offsets tell us the allocation size
      if (getProfileType() == StaticProfile) {
        ret.offset1D  = k_rowPtrs_(myRow);
        ret.allocSize = k_rowPtrs_(myRow+1) - k_rowPtrs_(myRow);
        if (k_numRowEntries_.dimension_0 () == 0) {
          ret.numEntries = ret.allocSize;
        } else {
          ret.numEntries = k_numRowEntries_(myRow);
        }
      }
      else {
        ret.offset1D = STINV;
        if (isLocallyIndexed ()) {
          ret.allocSize = lclInds2D_[myRow].size ();
        }
        else {
          ret.allocSize = gblInds2D_[myRow].size ();
        }
        ret.numEntries = k_numRowEntries_(myRow);
      }
    }
    else if (nodeNumAllocated_ == 0) {
      // have performed allocation, but the graph has no allocation or entries
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else if (! indicesAreAllocated ()) {
      // haven't performed allocation yet; probably won't hit this code
      //
      // FIXME (mfh 07 Aug 2014) We want graph's constructors to
      // allocate, rather than doing lazy allocation at first insert.
      // This will make k_numAllocPerRow_ obsolete.
      const bool useNumAllocPerRow = (k_numAllocPerRow_.dimension_0 () != 0);
      if (useNumAllocPerRow) {
        ret.allocSize = k_numAllocPerRow_(myRow); // this is a host View
      } else {
        ret.allocSize = numAllocForAllRows_;
      }
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else {
      // don't know how we ended up here...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  RowInfo
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getRowInfoFromGlobalRowIndex (const GlobalOrdinal gblRow) const
  {
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    if (! this->hasRowInfo () || this->rowMap_.is_null ()) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }
    const LocalOrdinal myRow = this->rowMap_->getLocalElement (gblRow);
    if (myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }

    ret.localRow = static_cast<size_t> (myRow);
    if (nodeNumAllocated_ != 0 && nodeNumAllocated_ != STINV) {
      // graph data structures have the info that we need
      //
      // if static graph, offsets tell us the allocation size
      if (getProfileType() == StaticProfile) {
        ret.offset1D  = k_rowPtrs_(myRow);
        ret.allocSize = k_rowPtrs_(myRow+1) - k_rowPtrs_(myRow);
        if (k_numRowEntries_.dimension_0 () == 0) {
          ret.numEntries = ret.allocSize;
        } else {
          ret.numEntries = k_numRowEntries_(myRow);
        }
      }
      else {
        ret.offset1D = STINV;
        if (isLocallyIndexed ()) {
          ret.allocSize = lclInds2D_[myRow].size ();
        }
        else {
          ret.allocSize = gblInds2D_[myRow].size ();
        }
        ret.numEntries = k_numRowEntries_(myRow);
      }
    }
    else if (nodeNumAllocated_ == 0) {
      // have performed allocation, but the graph has no allocation or entries
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else if (! indicesAreAllocated ()) {
      // haven't performed allocation yet; probably won't hit this code
      //
      // FIXME (mfh 07 Aug 2014) We want graph's constructors to
      // allocate, rather than doing lazy allocation at first insert.
      // This will make k_numAllocPerRow_ obsolete.
      const bool useNumAllocPerRow = (k_numAllocPerRow_.dimension_0 () != 0);
      if (useNumAllocPerRow) {
        ret.allocSize = k_numAllocPerRow_(myRow); // this is a host View
      } else {
        ret.allocSize = numAllocForAllRows_;
      }
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else {
      // don't know how we ended up here...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  staticAssertions () const
  {
    using Teuchos::OrdinalTraits;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;

    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal):
    //     This is so that we can store local indices in the memory
    //     formerly occupied by global indices.
    static_assert (sizeof (GlobalOrdinal) >= sizeof (LocalOrdinal),
                   "Tpetra::CrsGraph: sizeof(GlobalOrdinal) must be >= sizeof(LocalOrdinal).");
    // Assumption: max(size_t) >= max(LocalOrdinal)
    //     This is so that we can represent any LocalOrdinal as a size_t.
    static_assert (sizeof (size_t) >= sizeof (LocalOrdinal),
                   "Tpetra::CrsGraph: sizeof(size_t) must be >= sizeof(LocalOrdinal).");
    static_assert (sizeof(GST) >= sizeof(size_t),
                   "Tpetra::CrsGraph: sizeof(Tpetra::global_size_t) must be >= sizeof(size_t).");

    // FIXME (mfh 30 Sep 2015) We're not using
    // Teuchos::CompileTimeAssert any more.  Can we do these checks
    // with static_assert?

    // can't call max() with CompileTimeAssert, because it isn't a
    // constant expression; will need to make this a runtime check
    const char msg[] = "Tpetra::CrsGraph: Object cannot be created with the "
      "given template arguments: size assumptions are not valid.";
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (OrdinalTraits<LO>::max ()) > OrdinalTraits<size_t>::max (),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<GST> (OrdinalTraits<LO>::max ()) > static_cast<GST> (OrdinalTraits<GO>::max ()),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (OrdinalTraits<GO>::max ()) > OrdinalTraits<GST>::max(),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      OrdinalTraits<size_t>::max () > OrdinalTraits<GST>::max (),
      std::runtime_error, msg);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertIndices (const RowInfo& rowinfo,
                 const SLocalGlobalViews &newInds,
                 const ELocalGlobal lg,
                 const ELocalGlobal I)
  {
    using Teuchos::ArrayView;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      lg != GlobalIndices && lg != LocalIndices, std::invalid_argument,
      "Tpetra::CrsGraph::insertIndices: lg must be either GlobalIndices or "
      "LocalIndices.");
#endif // HAVE_TPETRA_DEBUG
    size_t numNewInds = 0;
    if (lg == GlobalIndices) { // input indices are global
      ArrayView<const GlobalOrdinal> new_ginds = newInds.ginds;
      numNewInds = new_ginds.size();
      if (I == GlobalIndices) { // store global indices
        ArrayView<GlobalOrdinal> gind_view = getGlobalViewNonConst(rowinfo);
        std::copy(new_ginds.begin(), new_ginds.end(), gind_view.begin()+rowinfo.numEntries);
      }
      else if (I == LocalIndices) { // store local indices
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        typename ArrayView<const GlobalOrdinal>::iterator         in = new_ginds.begin();
        const typename ArrayView<const GlobalOrdinal>::iterator stop = new_ginds.end();
        typename ArrayView<LocalOrdinal>::iterator out = lind_view.begin()+rowinfo.numEntries;
        while (in != stop) {
          *out++ = colMap_->getLocalElement (*in++);
        }
      }
    }
    else if (lg == LocalIndices) { // input indices are local
      ArrayView<const LocalOrdinal> new_linds = newInds.linds;
      numNewInds = new_linds.size();
      if (I == LocalIndices) { // store local indices
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        std::copy(new_linds.begin(), new_linds.end(), lind_view.begin()+rowinfo.numEntries);
      }
      else if (I == GlobalIndices) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Tpetra::CrsGraph::"
          "insertIndices: the case where the input indices are local and the "
          "indices to write are global (lg=LocalIndices, I=GlobalIndices) is "
          "not implemented, because it does not make sense." << std::endl <<
          "If you have correct local column indices, that means the graph has "
          "a column Map.  In that case, you should be storing local indices.");
      }
    }

    k_numRowEntries_(rowinfo.localRow) += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();
    return numNewInds;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertGlobalIndicesImpl (const LocalOrdinal myRow,
                           const Teuchos::ArrayView<const GlobalOrdinal>& indices)
  {
    const char tfecfFuncName[] = "insertGlobalIndicesImpl";

    RowInfo rowInfo = getRowInfo(myRow);
    const size_t numNewInds = indices.size();
    const size_t newNumEntries = rowInfo.numEntries + numNewInds;
    if (newNumEntries > rowInfo.allocSize) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() == StaticProfile, std::runtime_error,
        ": new indices exceed statically allocated graph structure.");

      // update allocation, doubling size to reduce # reallocations
      size_t newAllocSize = 2*rowInfo.allocSize;
      if (newAllocSize < newNumEntries) {
        newAllocSize = newNumEntries;
      }
      gblInds2D_[myRow].resize(newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);
    }

    // Copy new indices at end of global index array
    if (k_gblInds1D_.dimension_0 () != 0) {
      const size_t numIndsToCopy = static_cast<size_t> (indices.size ());
      const size_t offset = rowInfo.offset1D + rowInfo.numEntries;
      for (size_t k = 0; k < numIndsToCopy; ++k) {
        k_gblInds1D_[offset + k] = indices[k];
      }
    }
    else {
      std::copy(indices.begin(), indices.end(),
                gblInds2D_[myRow].begin()+rowInfo.numEntries);
    }

    k_numRowEntries_(myRow) += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();

#ifdef HAVE_TPETRA_DEBUG
    {
      const size_t chkNewNumEntries = getNumEntriesInLocalRow (myRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        chkNewNumEntries != newNumEntries, std::logic_error,
        ": Internal logic error. Please contact Tpetra team.");
    }
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertLocalIndicesImpl (const LocalOrdinal myRow,
                          const Teuchos::ArrayView<const LocalOrdinal>& indices)
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::subview;
    using Kokkos::View;
    typedef LocalOrdinal LO;
    const char* tfecfFuncName ("insertLocallIndicesImpl");

    const RowInfo rowInfo = this->getRowInfo(myRow);
    const size_t numNewInds = indices.size();
    const size_t newNumEntries = rowInfo.numEntries + numNewInds;
    if (newNumEntries > rowInfo.allocSize) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() == StaticProfile, std::runtime_error,
        ": new indices exceed statically allocated graph structure.");

      // update allocation, doubling size to reduce number of reallocations
      size_t newAllocSize = 2*rowInfo.allocSize;
      if (newAllocSize < newNumEntries)
        newAllocSize = newNumEntries;
      lclInds2D_[myRow].resize(newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);
    }

    // Store the new indices at the end of row myRow.
    if (k_lclInds1D_.dimension_0 () != 0) {
      typedef View<const LO*, execution_space, MemoryUnmanaged> input_view_type;
      typedef View<LO*, execution_space, MemoryUnmanaged> row_view_type;

      input_view_type inputInds (indices.getRawPtr (), indices.size ());
      const size_t start = rowInfo.offset1D + rowInfo.numEntries; // end of row
      const std::pair<size_t, size_t> rng (start, start + newNumEntries);
      // mfh 23 Nov 2015: Don't just create a subview of k_lclInds1D_
      // directly, because that first creates a _managed_ subview,
      // then returns an unmanaged version of that.  That touches the
      // reference count, which costs performance in a measurable way.
      row_view_type myInds = subview (row_view_type (k_lclInds1D_), rng);
      Kokkos::deep_copy (myInds, inputInds);
    }
    else {
      std::copy (indices.begin (), indices.end (),
                 lclInds2D_[myRow].begin () + rowInfo.numEntries);
    }

    k_numRowEntries_(myRow) += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();
#ifdef HAVE_TPETRA_DEBUG
    {
      const size_t chkNewNumEntries = getNumEntriesInLocalRow (myRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        chkNewNumEntries != newNumEntries, std::logic_error,
        ": Internal logic error. Please contact Tpetra team.");
    }
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  sortRowIndices (const RowInfo rowinfo)
  {
    if (rowinfo.numEntries > 0) {
      Teuchos::ArrayView<LocalOrdinal> inds_view =
        this->getLocalViewNonConst (rowinfo);
      std::sort (inds_view.begin (), inds_view.begin () + rowinfo.numEntries);
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  mergeRowIndices (RowInfo rowinfo)
  {
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "mergeRowIndices: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized (), std::logic_error, "The graph is already storage "
      "optimized, so we shouldn't be merging any indices.  "
      "Please report this bug to the Tpetra developers.");

    ArrayView<LocalOrdinal> inds_view = this->getLocalViewNonConst (rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = std::unique(beg,end);
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the
    // assignment below will destroy the packed structure
    TEUCHOS_TEST_FOR_EXCEPT( isStorageOptimized () && mergedEntries != rowinfo.numEntries );
#endif // HAVE_TPETRA_DEBUG

    k_numRowEntries_(rowinfo.localRow) = mergedEntries;
    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  setDomainRangeMaps (const Teuchos::RCP<const map_type>& domainMap,
                      const Teuchos::RCP<const map_type>& rangeMap)
  {
    // simple pointer comparison for equality
    if (domainMap_ != domainMap) {
      domainMap_ = domainMap;
      importer_ = null;
    }
    if (rangeMap_ != rangeMap) {
      rangeMap_  = rangeMap;
      exporter_ = null;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  clearGlobalConstants ()
  {
    globalNumEntries_       = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    globalNumDiags_         = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    globalMaxNumRowEntries_ = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    haveGlobalConstants_    = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  checkInternalState () const
  {
#ifdef HAVE_TPETRA_DEBUG
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const size_t         STI = Teuchos::OrdinalTraits<size_t>::invalid ();
    const char err[] = "Tpetra::CrsGraph::checkInternalState: Likely internal "
      "logic error.  Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object
    // always remains in a valid state
    // the graph should have been allocated with a row map
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap_ == null,                     std::logic_error, err );
    // am either complete or active
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() == isFillComplete(),  std::logic_error, err );
    // if the graph has been fill completed, then all maps should be present
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == true && (colMap_ == null || rangeMap_ == null || domainMap_ == null), std::logic_error, err );
    // if storage has been optimized, then indices should have been allocated (even if trivially so)
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && indicesAreAllocated() == false, std::logic_error, err );
    // if storage has been optimized, then number of allocated is now the number of entries
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && getNodeAllocationSize() != getNodeNumEntries(), std::logic_error, err );
    // if graph doesn't have the global constants, then they should all be marked as invalid
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == false && ( globalNumEntries_ != GSTI || globalNumDiags_ != GSTI || globalMaxNumRowEntries_ != GSTI ), std::logic_error, err );
    // if the graph has global cosntants, then they should be valid.
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ == GSTI || globalNumDiags_ == GSTI || globalMaxNumRowEntries_ == GSTI ), std::logic_error, err );
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ < nodeNumEntries_ || globalNumDiags_ < nodeNumDiags_ || globalMaxNumRowEntries_ < nodeMaxNumRowEntries_ ),
                        std::logic_error, err );
    // if indices are allocated, then the allocation specifications should have been released
    TEUCHOS_TEST_FOR_EXCEPTION(
      indicesAreAllocated () && (numAllocForAllRows_ != 0 ||
                                 k_numAllocPerRow_.dimension_0 () != 0),
      std::logic_error, err );
    // if indices are not allocated, then information dictating allocation quantities should be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesAreAllocated () && (nodeNumAllocated_ != STI ||
                                   nodeNumEntries_ != 0),
      std::logic_error, err );
    // if storage is optimized, then profile should be static
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() && pftype_ != StaticProfile,                                                               std::logic_error, err );

    // If k_rowPtrs_ exists (has nonzero size), it must have N+1 rows,
    // and k_rowPtrs_(N) must equal k_gblInds1D_.dimension_0() (if
    // globally indexed) or k_lclInds1D_.dimension_0() (if locally
    // indexed).

    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed () && k_rowPtrs_.dimension_0 () != 0 &&
      (static_cast<size_t> (k_rowPtrs_.dimension_0 ()) != getNodeNumRows () + 1 ||
       k_rowPtrs_(getNodeNumRows ()) != static_cast<size_t> (k_gblInds1D_.dimension_0 ())),
      std::logic_error, err );

    TEUCHOS_TEST_FOR_EXCEPTION(
      isLocallyIndexed () && k_rowPtrs_.dimension_0 () != 0 &&
      (static_cast<size_t> (k_rowPtrs_.dimension_0 ()) != getNodeNumRows () + 1 ||
       k_rowPtrs_(getNodeNumRows ()) != static_cast<size_t> (k_lclInds1D_.dimension_0 ())),
      std::logic_error, err );

    // If profile is dynamic and indices are allocated, then 2-D
    // allocations of column index storage (either local or global)
    // must be present.
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == DynamicProfile &&
      indicesAreAllocated () &&
      getNodeNumRows () > 0 &&
      lclInds2D_.is_null () && gblInds2D_.is_null (),
      std::logic_error, err );

    // If profile is dynamic and the calling process owns nonzero
    // rows, then k_numRowEntries_ and 2-D storage of column indices
    // (whether local or global) must be present.
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == DynamicProfile &&
      indicesAreAllocated () &&
      getNodeNumRows () > 0 &&
      (k_numRowEntries_.dimension_0 () == 0 ||
       (lclInds2D_.is_null () && gblInds2D_.is_null ())),
      std::logic_error, err );

    // if profile is dynamic, then 1D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == DynamicProfile &&
      (k_lclInds1D_.dimension_0 () != 0 || k_gblInds1D_.dimension_0 () != 0),
      std::logic_error, err );
    // if profile is dynamic, then row offsets should not be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == DynamicProfile && k_rowPtrs_.dimension_0 () != 0,
      std::logic_error, err );
    // if profile is static and we have allocated non-trivially, then
    // 1D allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == StaticProfile && indicesAreAllocated () &&
      getNodeAllocationSize () > 0 && k_lclInds1D_.dimension_0 () == 0 &&
      k_gblInds1D_.dimension_0 () == 0,
      std::logic_error, err);
    // if profile is static, then 2D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      pftype_ == StaticProfile && (lclInds2D_ != null || gblInds2D_ != null),
      std::logic_error, err );

    // if indices are not allocated, then none of the buffers should be.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! indicesAreAllocated () &&
      ((k_rowPtrs_.dimension_0 () != 0 || k_numRowEntries_.dimension_0 () != 0) ||
       k_lclInds1D_.dimension_0 () != 0 || lclInds2D_ != null ||
       k_gblInds1D_.dimension_0 () != 0 || gblInds2D_ != null),
      std::logic_error, err );

    // indices may be local or global only if they are allocated
    // (numAllocated is redundant; could simply be indicesAreLocal_ ||
    // indicesAreGlobal_)
    TEUCHOS_TEST_FOR_EXCEPTION( (indicesAreLocal_ || indicesAreGlobal_) && ! indicesAreAllocated_, std::logic_error, err );
    // indices may be local or global, but not both
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && indicesAreGlobal_ == true,                                                  std::logic_error, err );
    // if indices are local, then global allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      indicesAreLocal_ && (k_gblInds1D_.dimension_0 () != 0 || gblInds2D_ != null),
      std::logic_error, err );
    // if indices are global, then local allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      indicesAreGlobal_ && (k_lclInds1D_.dimension_0 () != 0 || lclInds2D_ != null),
      std::logic_error, err );
    // if indices are local, then local allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      indicesAreLocal_ && getNodeAllocationSize () > 0 &&
      k_lclInds1D_.dimension_0 () == 0 && getNodeNumRows () > 0 &&
      lclInds2D_.is_null (),
      std::logic_error, err);
    // if indices are global, then global allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION(
      indicesAreGlobal_ && getNodeAllocationSize () > 0 &&
      k_gblInds1D_.dimension_0 () == 0 && getNodeNumRows () > 0 &&
      gblInds2D_.is_null (),
      std::logic_error, err);
    // if indices are allocated, then we should have recorded how many were allocated
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == true  && nodeNumAllocated_ == STI,                                       std::logic_error, err );
    // if indices are not allocated, then the allocation size should be marked invalid
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && nodeNumAllocated_ != STI,                                       std::logic_error, err );
    // check the actual allocations
    if (indicesAreAllocated()) {
      size_t actualNumAllocated = 0;
      if (pftype_ == DynamicProfile) {
        if (isGloballyIndexed() && gblInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += gblInds2D_[r].size();
          }
        }
        else if (isLocallyIndexed() && lclInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += lclInds2D_[r].size();
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
      }
      else if (k_rowPtrs_.dimension_0 () != 0) { // pftype_ == StaticProfile
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_rowPtrs_.dimension_0 ()) != getNodeNumRows () + 1,
          std::logic_error, err);

        actualNumAllocated = k_rowPtrs_(getNodeNumRows ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          isLocallyIndexed () &&
          static_cast<size_t> (k_lclInds1D_.dimension_0 ()) != actualNumAllocated,
          std::logic_error, err );
        TEUCHOS_TEST_FOR_EXCEPTION(
          isGloballyIndexed () &&
          static_cast<size_t> (k_gblInds1D_.dimension_0 ()) != actualNumAllocated,
          std::logic_error, err );
        TEUCHOS_TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
      }
    }
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    const LocalOrdinal lrow = rowMap_->getLocalElement (globalRow);
    if (hasRowInfo () && lrow != OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowinfo = this->getRowInfo (lrow);
      return rowinfo.numEntries;
    } else {
      return OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNumEntriesInLocalRow (LocalOrdinal localRow) const
  {
    if (hasRowInfo () && rowMap_->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = this->getRowInfo (localRow);
      return rowinfo.numEntries;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNumAllocatedEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement (globalRow);
    if (hasRowInfo () && lrow != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowinfo = this->getRowInfo (lrow);
      return rowinfo.allocSize;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNumAllocatedEntriesInLocalRow (LocalOrdinal localRow) const
  {
    if (hasRowInfo () && rowMap_->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = this->getRowInfo (localRow);
      return rowinfo.allocSize;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayRCP<const size_t>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodeRowPtrs () const
  {
    using Kokkos::ViewAllocateWithoutInitializing;
    using Kokkos::create_mirror_view;
    using Teuchos::ArrayRCP;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::non_const_value_type row_offset_type;
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] = "Tpetra::CrsGraph::getNodeRowPtrs: ";
    const char suffix[] = "  Please report this bug to the Tpetra developers.";
#endif // HAVE_TPETRA_DEBUG
    const size_t size = k_rowPtrs_.dimension_0 ();
    const bool same = Kokkos::Impl::is_same<size_t, row_offset_type>::value;

    if (size == 0) {
      return ArrayRCP<const size_t> ();
    }

    ArrayRCP<const row_offset_type> ptr_rot;
    ArrayRCP<const size_t> ptr_st;
    if (same) { // size_t == row_offset_type
      // NOTE (mfh 22 Mar 2015) In a debug build of Kokkos, the result
      // of create_mirror_view might actually be a new allocation.
      // This helps with debugging when there are two memory spaces.
      typename row_map_type::HostMirror ptr_h = create_mirror_view (k_rowPtrs_);
      Kokkos::deep_copy (ptr_h, k_rowPtrs_);
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        ptr_h.dimension_0 () != k_rowPtrs_.dimension_0 (), std::logic_error,
        prefix << "size_t == row_offset_type, but ptr_h.dimension_0() = "
        << ptr_h.dimension_0 () << " != k_rowPtrs_.dimension_0() = "
        << k_rowPtrs_.dimension_0 () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        same && size != 0 && k_rowPtrs_.ptr_on_device () == NULL, std::logic_error,
        prefix << "size_t == row_offset_type and k_rowPtrs_.dimension_0() = "
        << size << " != 0, but k_rowPtrs_.ptr_on_device() == NULL." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION(
        same && size != 0 && ptr_h.ptr_on_device () == NULL, std::logic_error,
        prefix << "size_t == row_offset_type and k_rowPtrs_.dimension_0() = "
        << size << " != 0, but create_mirror_view(k_rowPtrs_).ptr_on_device() "
        "== NULL." << suffix);
#endif // HAVE_TPETRA_DEBUG
      ptr_rot = Kokkos::Compat::persistingView (ptr_h);
    }
    else { // size_t != row_offset_type
      typedef Kokkos::View<size_t*, device_type> ret_view_type;
      ret_view_type ptr_d (ViewAllocateWithoutInitializing ("ptr"), size);
      ::Tpetra::Details::copyOffsets (ptr_d, k_rowPtrs_);
      typename ret_view_type::HostMirror ptr_h = create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h, ptr_d);
      ptr_st = Kokkos::Compat::persistingView (ptr_h);
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      same && size != 0 && ptr_rot.is_null (), std::logic_error,
      prefix << "size_t == row_offset_type and size = " << size
      << " != 0, but ptr_rot is null." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! same && size != 0 && ptr_st.is_null (), std::logic_error,
      prefix << "size_t != row_offset_type and size = " << size
      << " != 0, but ptr_st is null." << suffix);
#endif // HAVE_TPETRA_DEBUG

    // If size_t == row_offset_type, return a persisting host view of
    // k_rowPtrs_.  Otherwise, return a size_t host copy of k_rowPtrs_.
#ifdef HAVE_TPETRA_DEBUG
    ArrayRCP<const size_t> retval =
      Kokkos::Impl::if_c<same,
        ArrayRCP<const row_offset_type>,
        ArrayRCP<const size_t> >::select (ptr_rot, ptr_st);
    TEUCHOS_TEST_FOR_EXCEPTION(
      size != 0 && retval.is_null (), std::logic_error,
      prefix << "size = " << size << " != 0, but retval is null." << suffix);
    return retval;
#else
    return Kokkos::Impl::if_c<same,
      ArrayRCP<const row_offset_type>,
      ArrayRCP<const size_t> >::select (ptr_rot, ptr_st);
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::ArrayRCP<const LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNodePackedIndices () const
  {
    return Kokkos::Compat::persistingView (k_lclInds1D_);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalRowCopy (LocalOrdinal localRow,
                   const Teuchos::ArrayView<LocalOrdinal>&indices,
                   size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "getLocalRowCopy: ";

    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed () && ! hasColMap (), std::runtime_error,
      "Tpetra::CrsGraph::getLocalRowCopy: The graph is globally indexed and "
      "does not have a column Map yet.  That means we don't have local indices "
      "for columns yet, so it doesn't make sense to call this method.  If the "
      "graph doesn't have a column Map yet, you should call fillComplete on "
      "it first.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo(), std::runtime_error,
      "Graph row information was deleted at fillComplete.");
#endif // HAVE_TPETRA_DEBUG

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowinfo = getRowInfo (localRow);
    // No side effects on error.
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < theNumEntries, std::runtime_error,
      "Specified storage (size==" << indices.size () << ") does not suffice "
      "to hold all " << theNumEntries << " entry/ies for this row.");
    numEntries = theNumEntries;

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (isLocallyIndexed ()) {
        ArrayView<const LO> lview = getLocalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = lview[j];
        }
      }
      else if (isGloballyIndexed ()) {
        ArrayView<const GO> gview = getGlobalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = colMap_->getLocalElement (gview[j]);
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalRowCopy (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<GlobalOrdinal>& indices,
                    size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "getGlobalRowCopy: ";
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      "Graph row information was deleted at fillComplete.");
#endif // HAVE_TPETRA_DEBUG

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowinfo = getRowInfoFromGlobalRowIndex (globalRow);
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < theNumEntries, std::runtime_error,
      "Specified storage (size==" << indices.size () << ") does not suffice "
      "to hold all " << theNumEntries << " entry/ies for this row.");
    numEntries = theNumEntries; // first side effect

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (isLocallyIndexed ()) {
        ArrayView<const LocalOrdinal> lview = getLocalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = colMap_->getGlobalElement (lview[j]);
        }
      }
      else if (isGloballyIndexed ()) {
        ArrayView<const GlobalOrdinal> gview = getGlobalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = gview[j];
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalRowView (LocalOrdinal localRow,
                   Teuchos::ArrayView<const LocalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getLocalRowView: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "The graph's indices are "
      "currently stored as global indices, so we cannot return a view with "
      "local column indices, whether or not the graph has a column Map.  If "
      "the graph _does_ have a column Map, use getLocalRowCopy() instead.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error, "Graph row information was "
      "deleted at fillComplete().");
#endif // HAVE_TPETRA_DEBUG

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowInfo = getRowInfo (localRow);
    indices = Teuchos::null;
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = this->getLocalView (rowInfo);
      // getLocalView returns a view of the _entire_ row, including
      // any extra space at the end (which 1-D unpacked storage
      // might have, for example).  That's why we have to take a
      // subview of the returned view.
      indices = indices (0, rowInfo.numEntries);
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) !=
       getNumEntriesInLocalRow (localRow), std::logic_error, "indices.size() "
       "= " << indices.size () << " != getNumEntriesInLocalRow(localRow=" <<
       localRow << ") = " << getNumEntriesInLocalRow (localRow) <<
       ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getGlobalRowView (GlobalOrdinal globalRow,
                    Teuchos::ArrayView<const GlobalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getGlobalRowView: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error, "The graph's indices are "
      "currently stored as local indices, so we cannot return a view with "
      "global column indices.  Use getGlobalRowCopy() instead.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      "Graph row information was deleted at fillComplete().");
#endif // HAVE_TPETRA_DEBUG

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowInfo = getRowInfoFromGlobalRowIndex (globalRow);
    indices = Teuchos::null;
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = (this->getGlobalView (rowInfo)) (0, rowInfo.numEntries);
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) != getNumEntriesInGlobalRow (globalRow),
       std::logic_error, "indices.size() = " << indices.size ()
       << " != getNumEntriesInGlobalRow(globalRow=" << globalRow << ") = "
       << getNumEntriesInGlobalRow (globalRow)
       << ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertLocalIndices (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal>& indices)
  {
    const char tfecfFuncName[] = "insertLocalIndices";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error,
      ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error,
      ": graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error,
      ": cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! rowMap_->isNodeLocalElement (localRow), std::runtime_error,
      ": row does not belong to this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices);
    }

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, if the graph has a column Map, test whether
    // any of the given column indices are not in the column Map.
    // Keep track of the invalid column indices so we can tell the
    // user about them.
    if (hasColMap ()) {
      using Teuchos::Array;
      using Teuchos::toString;
      using std::endl;
      typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type size_type;

      const map_type& colMap = * (getColMap ());
      Array<LocalOrdinal> badColInds;
      bool allInColMap = true;
      for (size_type k = 0; k < indices.size (); ++k) {
        if (! colMap.isNodeLocalElement (indices[k])) {
          allInColMap = false;
          badColInds.push_back (indices[k]);
        }
      }
      if (! allInColMap) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::insertLocalIndices: You attempted to insert "
          "entries in owned row " << localRow << ", at the following column "
          "indices: " << toString (indices) << "." << endl;
        os << "Of those, the following indices are not in the column Map on "
          "this process: " << toString (badColInds) << "." << endl << "Since "
          "the graph has a column Map already, it is invalid to insert entries "
          "at those locations.";
        TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
      }
    }
#endif // HAVE_TPETRA_DEBUG

    insertLocalIndicesImpl (localRow, indices);

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isLocallyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertLocalIndices (const LocalOrdinal localRow,
                      const LocalOrdinal numEnt,
                      const LocalOrdinal inds[])
  {
    Teuchos::ArrayView<const LocalOrdinal> indsT (inds, numEnt);
    this->insertLocalIndices (localRow, indsT);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertLocalIndicesFiltered (const LocalOrdinal localRow,
                              const Teuchos::ArrayView<const LocalOrdinal>& indices)
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "insertLocalIndicesFiltered";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed() == true, std::runtime_error,
      ": graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasColMap() == false, std::runtime_error,
      ": cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
      ": row does not belong to this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices);
    }

     // If we have a column map, use it to filter the entries.
    if (hasColMap ()) {
      Teuchos::Array<LO> filtered_indices (indices);
      SLocalGlobalViews inds_view;
      SLocalGlobalNCViews inds_ncview;
      inds_ncview.linds = filtered_indices();
      const size_t numFilteredEntries =
        filterIndices<LocalIndices>(inds_ncview);
      inds_view.linds = filtered_indices (0, numFilteredEntries);
      insertLocalIndicesImpl(localRow, inds_view.linds);
    }
    else {
      insertLocalIndicesImpl(localRow, indices);
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isLocallyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertGlobalIndices (const GlobalOrdinal grow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& indices)
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;
    const char tfecfFuncName[] = "insertGlobalIndices";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed() == true, std::runtime_error,
      ": graph indices are local; use insertLocalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": You are not allowed to call this method if fill is not active.  "
      "If fillComplete has been called, you must first call resumeFill "
      "before you may insert indices.");
    if (! indicesAreAllocated ()) {
      allocateIndices (GlobalIndices);
    }
    const LO myRow = rowMap_->getLocalElement (grow);
    if (myRow != Teuchos::OrdinalTraits<LO>::invalid ()) {
#ifdef HAVE_TPETRA_DEBUG
      if (hasColMap ()) {
        using std::endl;
        const map_type& colMap = * (getColMap ());
        // In a debug build, keep track of the nonowned ("bad") column
        // indices, so that we can display them in the exception
        // message.  In a release build, just ditch the loop early if
        // we encounter a nonowned column index.
        Array<GO> badColInds;
        bool allInColMap = true;
        for (size_type k = 0; k < indices.size (); ++k) {
          if (! colMap.isNodeGlobalElement (indices[k])) {
            allInColMap = false;
            badColInds.push_back (indices[k]);
          }
        }
        if (! allInColMap) {
          std::ostringstream os;
          os << "Tpetra::CrsGraph::insertGlobalIndices: You attempted to insert "
            "entries in owned row " << grow << ", at the following column "
            "indices: " << toString (indices) << "." << endl;
          os << "Of those, the following indices are not in the column Map on "
            "this process: " << toString (badColInds) << "." << endl << "Since "
            "the matrix has a column Map already, it is invalid to insert "
            "entries at those locations.";
          TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
        }
      }
#endif // HAVE_TPETRA_DEBUG
      insertGlobalIndicesImpl (myRow, indices);
    }
    else { // a nonlocal row
      const size_type numIndices = indices.size ();
      // This creates the Array if it doesn't exist yet.  std::map's
      // operator[] does a lookup each time, so it's better to pull
      // nonlocals_[grow] out of the loop.
      std::vector<GO>& nonlocalRow = nonlocals_[grow];
      for (size_type k = 0; k < numIndices; ++k) {
        nonlocalRow.push_back (indices[k]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isGloballyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertGlobalIndices (const GlobalOrdinal globalRow,
                       const LocalOrdinal numEnt,
                       const GlobalOrdinal inds[])
  {
    Teuchos::ArrayView<const GlobalOrdinal> indsT (inds, numEnt);
    this->insertGlobalIndices (globalRow, indsT);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  insertGlobalIndicesFiltered (const GlobalOrdinal grow,
                               const Teuchos::ArrayView<const GlobalOrdinal>& indices)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "insertGlobalIndicesFiltered";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed() == true, std::runtime_error,
      ": graph indices are local; use insertLocalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": You are not allowed to call this method if fill is not active.  "
      "If fillComplete has been called, you must first call resumeFill "
      "before you may insert indices.");
    if (! indicesAreAllocated ()) {
      allocateIndices (GlobalIndices);
    }
    const LO myRow = rowMap_->getLocalElement (grow);
    if (myRow != Teuchos::OrdinalTraits<LO>::invalid ()) {
      // If we have a column map, use it to filter the entries.
      if (hasColMap ()) {
        Array<GO> filtered_indices(indices);
        SLocalGlobalViews inds_view;
        SLocalGlobalNCViews inds_ncview;
        inds_ncview.ginds = filtered_indices();
        const size_t numFilteredEntries =
          filterIndices<GlobalIndices> (inds_ncview);
        inds_view.ginds = filtered_indices (0, numFilteredEntries);
        insertGlobalIndicesImpl(myRow, inds_view.ginds);
      }
      else {
       insertGlobalIndicesImpl(myRow, indices);
      }
    }
    else { // a nonlocal row
      typedef typename ArrayView<const GO>::size_type size_type;
      const size_type numIndices = indices.size ();
      for (size_type k = 0; k < numIndices; ++k) {
        nonlocals_[grow].push_back (indices[k]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isGloballyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  removeLocalIndices (LocalOrdinal lrow)
  {
    const char tfecfFuncName[] = "removeLocalIndices: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error, "requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized (), std::runtime_error,
      "cannot remove indices after optimizeStorage() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "graph indices are global.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! rowMap_->isNodeLocalElement (lrow), std::runtime_error,
      "Local row " << lrow << " is not in the row Map on the calling process.");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices);
    }

    // FIXME (mfh 13 Aug 2014) What if they haven't been cleared on
    // all processes?
    clearGlobalConstants ();

    if (k_numRowEntries_.dimension_0 () != 0) {
      const size_t oldNumEntries = k_numRowEntries_(lrow);
      nodeNumEntries_ -= oldNumEntries;
      k_numRowEntries_(lrow) = 0;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNumEntriesInLocalRow (lrow) != 0 ||
      ! indicesAreAllocated () ||
      ! isLocallyIndexed (), std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  setAllIndices (const typename local_graph_type::row_map_type& rowPointers,
                 const typename local_graph_type::entries_type::non_const_type& columnIndices)
  {
    const char tfecfFuncName[] = "setAllIndices: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap () || getColMap ().is_null (), std::runtime_error,
      "The graph must have a column Map before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (rowPointers.size ()) != this->getNodeNumRows () + 1,
      std::runtime_error, "rowPointers.size() = " << rowPointers.size () <<
      " != this->getNodeNumRows()+1 = " << (this->getNodeNumRows () + 1) <<
      ".");

    // FIXME (mfh 07 Aug 2014) We need to relax this restriction,
    // since the future model will be allocation at construction, not
    // lazy allocation on first insert.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_lclInds1D_.dimension_0 () != 0 || k_gblInds1D_.dimension_0 () != 0,
      std::runtime_error, "You may not call this method if 1-D data structures "
      "are already allocated.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclInds2D_ != Teuchos::null || gblInds2D_ != Teuchos::null,
      std::runtime_error, "You may not call this method if 2-D data structures "
      "are already allocated.");

    const size_t localNumEntries = rowPointers(getNodeNumRows ());

    indicesAreAllocated_ = true;
    indicesAreLocal_     = true;
    pftype_              = StaticProfile; // if the profile wasn't static before, it sure is now.
    k_lclInds1D_         = columnIndices;
    k_rowPtrs_           = rowPointers;
    nodeNumAllocated_    = localNumEntries;
    nodeNumEntries_      = localNumEntries;

    // Build the local graph.
    lclGraph_ = local_graph_type (k_lclInds1D_, k_rowPtrs_);

    // These normally get cleared out at the end of allocateIndices.
    // It makes sense to clear them out here, because at the end of
    // this method, the graph is allocated on the calling process.
    numAllocForAllRows_ = 0;
    k_numAllocPerRow_ = decltype (k_numAllocPerRow_) ();

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  setAllIndices (const Teuchos::ArrayRCP<size_t>& rowPointers,
                 const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices)
  {
    using Kokkos::View;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::array_layout layout_type;
    typedef typename row_map_type::non_const_value_type row_offset_type;
    typedef View<size_t*, layout_type , Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged> input_view_type;
    typedef typename row_map_type::non_const_type nc_row_map_type;

    const size_t size = static_cast<size_t> (rowPointers.size ());
    const bool same = Kokkos::Impl::is_same<size_t, row_offset_type>::value;
    input_view_type ptr_in (rowPointers.getRawPtr (), size);

    nc_row_map_type ptr_rot ("Tpetra::CrsGraph::ptr", size);

    if (same) { // size_t == row_offset_type
      // This compile-time logic ensures that the compiler never sees
      // an assignment of View<row_offset_type*, ...> to View<size_t*,
      // ...> unless size_t == row_offset_type.
      input_view_type ptr_decoy (rowPointers.getRawPtr (), size); // never used
      Kokkos::deep_copy (Kokkos::Impl::if_c<same,
                           nc_row_map_type,
                           input_view_type>::select (ptr_rot, ptr_decoy),
                         ptr_in);
    }
    else { // size_t != row_offset_type
      // CudaUvmSpace != HostSpace, so this will be false in that case.
      const bool inHostMemory =
        Kokkos::Impl::is_same<typename row_map_type::memory_space,
          Kokkos::HostSpace>::value;
      if (inHostMemory) {
        // Copy (with cast from size_t to row_offset_type, with bounds
        // checking if necessary) to ptr_rot.
        ::Tpetra::Details::copyOffsets (ptr_rot, ptr_in);
      }
      else { // Copy input row offsets to device first.
        //
        // FIXME (mfh 24 Mar 2015) If CUDA UVM, running in the host's
        // execution space would avoid the double copy.
        //
        View<size_t*, layout_type ,execution_space > ptr_st ("Tpetra::CrsGraph::ptr", size);
        Kokkos::deep_copy (ptr_st, ptr_in);
        // Copy on device (casting from size_t to row_offset_type,
        // with bounds checking if necessary) to ptr_rot.  This
        // executes in the output View's execution space, which is the
        // same as execution_space.
        ::Tpetra::Details::copyOffsets (ptr_rot, ptr_st);
      }
    }

    Kokkos::View<LocalOrdinal*, layout_type , execution_space > k_ind =
      Kokkos::Compat::getKokkosViewDeepCopy<device_type> (columnIndices ());
    setAllIndices (ptr_rot, k_ind);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getNumEntriesPerLocalRowUpperBound (Teuchos::ArrayRCP<const size_t>& boundPerLocalRow,
                                      size_t& boundForAllLocalRows,
                                      bool& boundSameForAllLocalRows) const
  {
    // The three output arguments.  We assign them to the actual
    // output arguments at the end, in order to implement
    // transactional semantics.
    Teuchos::ArrayRCP<const size_t> numEntriesPerRow;
    size_t numEntriesForAll = 0;
    bool allRowsSame = true;

    const ptrdiff_t numRows = static_cast<ptrdiff_t> (this->getNodeNumRows ());

    if (! this->indicesAreAllocated ()) {
      if (k_numAllocPerRow_.dimension_0 () != 0) {
        // This is a shallow copy; the ArrayRCP wraps the View in a
        // custom destructor, which ensures correct deallocation if
        // that is the only reference to the View.  Furthermore, this
        // View is a host View, so this doesn't assume UVM.
        numEntriesPerRow = Kokkos::Compat::persistingView (k_numAllocPerRow_);
        allRowsSame = false; // conservatively; we don't check the array
      }
      else {
        numEntriesForAll = numAllocForAllRows_;
        allRowsSame = true;
      }
    }
    else if (k_numRowEntries_.dimension_0 () != 0) {
      // This is a shallow copy; the ArrayRCP wraps the View in a
      // custom destructor, which ensures correct deallocation if that
      // is the only reference to the View.  Furthermore, this View is
      // a host View, so this doesn't assume UVM.
      numEntriesPerRow = Kokkos::Compat::persistingView (k_numRowEntries_);
      allRowsSame = false; // conservatively; we don't check the array
    }
    else if (this->nodeNumAllocated_ == 0) {
      numEntriesForAll = 0;
      allRowsSame = true;
    }
    else {
      // left with the case that we have optimized storage. in this
      // case, we have to construct a list of row sizes.
      TEUCHOS_TEST_FOR_EXCEPTION(
        this->getProfileType () != StaticProfile, std::logic_error,
        "Tpetra::CrsGraph::getNumEntriesPerRowUpperBound: "
        "The graph is not StaticProfile, but storage appears to be optimized.  "
        "Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        numRows != 0 && k_rowPtrs_.dimension_0 () == 0, std::logic_error,
        "Tpetra::CrsGraph::getNumEntriesPerRowUpperBound: "
        "The graph has " << numRows << " (> 0) row" << (numRows != 1 ? "s" : "")
        << " on the calling process, but the k_rowPtrs_ array has zero entries.  "
        "Please report this bug to the Tpetra developers.");

      Teuchos::ArrayRCP<size_t> numEnt;
      if (numRows != 0) {
        numEnt = Teuchos::arcp<size_t> (numRows);
      }

      // We have to iterate through the row offsets anyway, so we
      // might as well check whether all rows' bounds are the same.
      bool allRowsReallySame = false;
      for (ptrdiff_t i = 0; i < numRows; ++i) {
        numEnt[i] = k_rowPtrs_(i+1) - k_rowPtrs_(i);
        if (i != 0 && numEnt[i] != numEnt[i-1]) {
          allRowsReallySame = false;
        }
      }
      if (allRowsReallySame) {
        if (numRows == 0) {
          numEntriesForAll = 0;
        } else {
          numEntriesForAll = numEnt[1] - numEnt[0];
        }
        allRowsSame = true;
      }
      else {
        numEntriesPerRow = numEnt; // Teuchos::arcp_const_cast<const size_t> (numEnt);
        allRowsSame = false; // conservatively; we don't check the array
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      numEntriesForAll != 0 && numEntriesPerRow.size () != 0, std::logic_error,
      "Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound: "
      "numEntriesForAll and numEntriesPerRow are not consistent.  The former "
      "is nonzero (" << numEntriesForAll << "), but the latter has nonzero "
      "size " << numEntriesPerRow.size () << ".  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numEntriesForAll != 0 && ! allRowsSame, std::logic_error,
      "Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound: "
      "numEntriesForAll and allRowsSame are not consistent.  The former "
      "is nonzero (" << numEntriesForAll << "), but the latter is false.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numEntriesPerRow.size () != 0 && allRowsSame, std::logic_error,
      "Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound: "
      "numEntriesPerRow and allRowsSame are not consistent.  The former has "
      "nonzero length " << numEntriesForAll << ", but the latter is true.  "
      "Please report this bug to the Tpetra developers.");

    boundPerLocalRow = numEntriesPerRow;
    boundForAllLocalRows = numEntriesForAll;
    boundSameForAllLocalRows = allRowsSame;
  }


  // TODO: in the future, globalAssemble() should use import/export functionality
  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  globalAssemble ()
  {
    using Teuchos::Array;
    using Teuchos::Comm;
    using Teuchos::gatherAll;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using Teuchos::toString;
    using Teuchos::waitAll;
    using std::endl;
    using std::make_pair;
    using std::pair;
    typedef GlobalOrdinal GO;
    typedef typename std::map<GO, std::vector<GO> >::const_iterator NLITER;
    typedef typename Array<GO>::size_type size_type;

    const char tfecfFuncName[] = "globalAssemble"; // for exception macro
    RCP<const Comm<int> > comm = getComm();

    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive(), std::runtime_error,
      ": requires that fill is active.");
    // Determine if any nodes have global entries to share
    {
      size_t MyNonlocals = nonlocals_.size(), MaxGlobalNonlocals;
      reduceAll (*comm, REDUCE_MAX, MyNonlocals, outArg (MaxGlobalNonlocals));
      if (MaxGlobalNonlocals == 0) {
        return;  // no entries to share
      }
    }

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images:
    //   globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int, GO> > IdsAndRows;
    std::map<GO, int> NLR2Id;
    Teuchos::SerialDenseMatrix<int, char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all
      // nonowned rows.  Compute list of rows for which we have data.
      Array<GO> NLRs;
      std::set<GO> setOfRows;
      for (NLITER iter = nonlocals_.begin (); iter != nonlocals_.end (); ++iter) {
        setOfRows.insert (iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize (setOfRows.size ());
      std::copy (setOfRows.begin (), setOfRows.end (), NLRs.begin ());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        const LookupStatus stat =
          rowMap_->getRemoteIndexList (NLRs (), NLRIds ());
        int lclerror = ( stat == IDNotPresent ? 1 : 0 );
        int gblerror;
        reduceAll<int, int> (*comm, REDUCE_MAX, lclerror, outArg (gblerror));
        if (gblerror != 0) {
          const int myRank = comm->getRank ();
          std::ostringstream os;
          os << "On one or more processes in the communicator, "
             << "there were insertions into rows of the graph that do not "
             << "exist in the row Map on any process in the communicator."
             << endl << "This process " << myRank << " is "
             << (lclerror == 0 ? "not " : "") << "one of those offending "
             << "processes." << endl;
          if (lclerror != 0) {
            // If NLRIds[k] is -1, then NLRs[k] is a row index not in
            // the row Map.  Collect this list of invalid row indices
            // for display in the exception message.
            Array<GO> invalidNonlocalRows;
            for (size_type k = 0; k < NLRs.size (); ++k) {
              if (NLRIds[k] == -1) {
                invalidNonlocalRows.push_back (NLRs[k]);
              }
            }
            const size_type numInvalid = invalidNonlocalRows.size ();
            os << "On this process, " << numInvalid << " nonlocal row"
               << (numInvalid != 1 ? "s " : " ") << " were inserted that are "
               << "not in the row Map on any process." << endl;
            // Don't print _too_ many nonlocal rows.
            if (numInvalid <= 100) {
              os << "Offending row indices: "
                 << toString (invalidNonlocalRows ()) << endl;
            }
          }
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            gblerror != 0, std::runtime_error,
            ": nonlocal entries correspond to invalid rows."
            << endl << os.str ());
        }
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages, 0);
      typename Array<GO>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        // IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j) {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      gatherAll (*getComm(), numImages, localNeighbors.getRawPtr(),
                 numImages * numImages, globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images.
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    //////////////////////////////////////////////////////////////////////////////////////

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j) {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    const size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GO, GO> > IJSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int, GO> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) {
      int id = IdAndRow->first;
      GO row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(sendIDs[numSends] != id,
          std::logic_error, ": internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename std::vector<GO>::const_iterator j = nonlocals_[row].begin ();
           j != nonlocals_[row].end (); ++j) {
        IJSendBuffer.push_back (pair<GlobalOrdinal, GlobalOrdinal> (row, *j));
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<typename Array<int>::size_type> (numSends) != sendIDs.size (),
      std::logic_error, ": internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    //////////////////////////////////////////////////////////////////////////////////////
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    //////////////////////////////////////////////////////////////////////////////////////

    // Array of pending nonblocking communication requests.  It's OK
    // to mix nonblocking send and receive requests in the same
    // waitAll() call.
    Array<RCP<Teuchos::CommRequest<int> > > requests;

    // perform non-blocking sends: send sizes to our recipients
    for (size_t s = 0; s < numSends ; ++s) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (isend<int, size_t> (*comm,
                                              rcp (&sendSizes[s], false),
                                              sendIDs[s]));
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<size_t> recvSizes (numRecvs);
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (ireceive (*comm, rcp (&recvSizes[r], false), recvIDs[r]));
    }
    // Wait on all the nonblocking sends and receives.
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    // This doesn't necessarily deallocate the array.
    requests.resize (0);

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJSendBuffer
    Array<ArrayView<pair<GO,GO> > > sendBuffers(numSends,null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpSendBuf =
        arcp (sendBuffers[s].getRawPtr(), 0, sendBuffers[s].size(), false);
      requests.push_back (isend<int, pair<GO,GO> > (*comm, tmpSendBuf, sendIDs[s]));
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate (recvSizes.begin(), recvSizes.end(), 0);
    Array<pair<GO,GO> > IJRecvBuffer (totalRecvSize);
    // from the size info, build the ArrayViews into IJRecvBuffer
    Array<ArrayView<pair<GO,GO> > > recvBuffers (numRecvs, null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpRecvBuf =
        arcp (recvBuffers[r].getRawPtr(), 0, recvBuffers[r].size(), false);
      requests.push_back (ireceive (*comm, tmpRecvBuf, recvIDs[r]));
    }
    // perform waits
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node,
    //       so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<pair<GO,GO> >::const_iterator ij = IJRecvBuffer.begin();
         ij != IJRecvBuffer.end(); ++ij)
    {
      insertGlobalIndicesFiltered (ij->first, tuple<GO> (ij->second));
    }
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "resumeFill";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! hasRowInfo(), std::runtime_error,
      ": Sorry, you cannot resume fill of the CrsGraph, since the graph's row "
      "information was deleted in fillComplete().");

#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *rowMap_->getComm() );
#endif // HAVE_TPETRA_DEBUG
    clearGlobalConstants();
    if (params != null) this->setParameterList (params);
    lowerTriangular_  = false;
    upperTriangular_  = false;
    // either still sorted/merged or initially sorted/merged
    indicesAreSorted_ = true;
    noRedundancies_ = true;
    fillComplete_ = false;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive() || isFillComplete(), std::logic_error,
      "::resumeFill(): At end of method, either fill is not active or fill is "
      "complete.  This violates stated post-conditions.  Please report this bug "
      "to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    // If the graph already has domain and range Maps, don't clobber
    // them.  If it doesn't, use the current row Map for both the
    // domain and range Maps.
    //
    // NOTE (mfh 28 Sep 2014): If the graph was constructed without a
    // column Map, and column indices are inserted which are not in
    // the row Map on any process, this will cause troubles.  However,
    // that is not a common case for most applications that we
    // encounter, and checking for it might require more
    // communication.
    Teuchos::RCP<const map_type> domMap = this->getDomainMap ();
    if (domMap.is_null ()) {
      domMap = this->getRowMap ();
    }
    Teuchos::RCP<const map_type> ranMap = this->getRangeMap ();
    if (ranMap.is_null ()) {
      ranMap = this->getRowMap ();
    }
    this->fillComplete (domMap, ranMap, params);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                const Teuchos::RCP<const map_type>& rangeMap,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "fillComplete";

#ifdef HAVE_TPETRA_DEBUG
    rowMap_->getComm ()->barrier ();
#endif // HAVE_TPETRA_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, ": Graph fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");

    const int numProcs = getComm ()->getSize ();

    //
    // Read and set parameters
    //

    // Does the caller want to sort remote GIDs (within those owned by
    // the same process) in makeColMap()?
    if (! params.is_null ()) {
      if (params->isParameter ("sort column map ghost gids")) {
        sortGhostsAssociatedWithEachProcessor_ =
          params->get<bool> ("sort column map ghost gids",
                             sortGhostsAssociatedWithEachProcessor_);
      }
      else if (params->isParameter ("Sort column Map ghost GIDs")) {
        sortGhostsAssociatedWithEachProcessor_ =
          params->get<bool> ("Sort column Map ghost GIDs",
                             sortGhostsAssociatedWithEachProcessor_);
      }
    }

    // If true, the caller promises that no process did nonlocal
    // changes since the last call to fillComplete.
    bool assertNoNonlocalInserts = false;
    if (! params.is_null ()) {
      assertNoNonlocalInserts =
        params->get<bool> ("No Nonlocal Changes", assertNoNonlocalInserts);
    }

    //
    // Allocate indices, if they haven't already been allocated
    //
    if (! indicesAreAllocated ()) {
      if (hasColMap ()) {
        // We have a column Map, so use local indices.
        allocateIndices (LocalIndices);
      } else {
        // We don't have a column Map, so use global indices.
        allocateIndices (GlobalIndices);
      }
    }

    //
    // Do global assembly, if requested and if the communicator
    // contains more than one process.
    //
    const bool mayNeedGlobalAssemble = ! assertNoNonlocalInserts && numProcs > 1;
    if (mayNeedGlobalAssemble) {
      // This first checks if we need to do global assembly.
      // The check costs a single all-reduce.
      globalAssemble ();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numProcs > 1 && nonlocals_.size() > 0, std::runtime_error,
        ":" << std::endl << "The graph's communicator contains only one "
        "process, but there are nonlocal entries.  " << std::endl <<
        "This probably means that invalid entries were added to the graph.");
    }

    // Set domain and range Map.  This may clear the Import / Export
    // objects if the new Maps differ from any old ones.
    setDomainRangeMaps (domainMap, rangeMap);

    // If the graph does not already have a column Map (either from
    // the user constructor calling the version of the constructor
    // that takes a column Map, or from a previous fillComplete call),
    // then create it.
    if (! hasColMap ()) {
      makeColMap ();
    }

    // Make indices local, if they aren't already.
    // The method doesn't do any work if the indices are already local.
    makeIndicesLocal ();

    if (! isSorted ()) {
      // If this process has no indices, then CrsGraph considers it
      // already trivially sorted.  Thus, this method need not be
      // called on all processes in the row Map's communicator.
      sortAllIndices ();
    }

    if (! isMerged()) {
      mergeAllIndices ();
    }
    makeImportExport (); // Make Import and Export objects, if necessary
    computeGlobalConstants ();
    fillLocalGraph (params);
    fillComplete_ = true;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == true || isFillComplete() == false, std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                            const Teuchos::RCP<const map_type>& rangeMap,
                            const Teuchos::RCP<const import_type>& importer,
                            const Teuchos::RCP<const export_type>& exporter,
                            const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "expertStaticFillComplete: ";
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string label;
    if(!params.is_null())
      label = params->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-Setup"))));
#endif


    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      domainMap.is_null () || rangeMap.is_null (),
      std::runtime_error, "The input domain Map and range Map must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      pftype_ != StaticProfile, std::runtime_error, "You may not call this "
      "method unless the graph is StaticProfile.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete () || ! hasColMap (), std::runtime_error, "You may not "
      "call this method unless the graph has a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNodeNumRows () > 0 && k_rowPtrs_.dimension_0 () == 0,
      std::runtime_error, "The calling process has getNodeNumRows() = "
      << getNodeNumRows () << " > 0 rows, but the row offsets array has not "
      "been set.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (k_rowPtrs_.dimension_0 ()) != getNodeNumRows () + 1,
      std::runtime_error, "The row offsets array has length " <<
      k_rowPtrs_.dimension_0 () << " != getNodeNumRows()+1 = " <<
      (getNodeNumRows () + 1) << ".");

    // Note: We don't need to do the following things which are normally done in fillComplete:
    // allocateIndices, globalAssemble, makeColMap, makeIndicesLocal, sortAllIndices, mergeAllIndices

    // Note: Need to do this so computeGlobalConstants & fillLocalGraph work
    //
    // The first assignment is always true if the graph has 1-D
    // storage (StaticProfile).  The second assignment is only true if
    // storage is packed.
    nodeNumAllocated_ = k_rowPtrs_(getNodeNumRows ());
    nodeNumEntries_ = nodeNumAllocated_;

    // Constants from allocateIndices
    //
    // mfh 08 Aug 2014: numAllocForAllRows_ and k_numAllocPerRow_ go
    // away once the graph is allocated.  expertStaticFillComplete
    // either presumes that the graph is allocated, or "allocates" it.
    //
    // FIXME (mfh 08 Aug 2014) The goal for the Kokkos refactor
    // version of CrsGraph is to allocate in the constructor, not
    // lazily on first insert.  That will make both
    // numAllocForAllRows_ and k_numAllocPerRow_ obsolete.
    numAllocForAllRows_  = 0;
    k_numAllocPerRow_    = decltype (k_numAllocPerRow_) ();
    indicesAreAllocated_ = true;

    // Constants from makeIndicesLocal
    //
    // The graph has a column Map, so its indices had better be local.
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;

    // set domain/range map: may clear the import/export objects
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-Maps"))));
#endif
    setDomainRangeMaps (domainMap, rangeMap);

    // Presume the user sorted and merged the arrays first
    indicesAreSorted_ = true;
    noRedundancies_ = true;

    // makeImportExport won't create a new importer/exporter if I set one here first.
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXcheckI"))));
#endif

    importer_ = Teuchos::null;
    exporter_ = Teuchos::null;
    if (importer != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! importer->getSourceMap ()->isSameAs (*getDomainMap ()) ||
        ! importer->getTargetMap ()->isSameAs (*getColMap ()),
        std::invalid_argument,": importer does not match matrix maps.");
      importer_ = importer;

    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXcheckE"))));
#endif

    if (exporter != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! exporter->getSourceMap ()->isSameAs (*getRowMap ()) ||
        ! exporter->getTargetMap ()->isSameAs (*getRangeMap ()),
        std::invalid_argument,": exporter does not match matrix maps.");
      exporter_ = exporter;
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXmake"))));
#endif

    makeImportExport ();

    // Compute the constants
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-cGC"))));
#endif
    computeGlobalConstants ();

    // Since we have a StaticProfile, fillLocalGraph will do the right thing...
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-fLG"))));
#endif
    fillLocalGraph (params);
    fillComplete_ = true;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-cIS"))));
#endif
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::computeOffsetsFromCounts;
    using Kokkos::create_mirror_view;
    typedef decltype (k_numRowEntries_) row_entries_type;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::non_const_type non_const_row_map_type;
    typedef typename local_graph_type::entries_type::non_const_type lclinds_1d_type;
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "fillLocalGraph (called from fillComplete or "
      "expertStaticFillComplete): ";
#endif // HAVE_TPETRA_DEBUG
    const size_t lclNumRows = this->getNodeNumRows ();

    // This method's goal is to fill in the two arrays (compressed
    // sparse row format) that define the sparse graph's structure.
    //
    // Use the nonconst version of row_map_type for ptr_d, because
    // the latter is const and we need to modify ptr_d here.
    non_const_row_map_type ptr_d;
    row_map_type ptr_d_const;
    lclinds_1d_type ind_d;

    bool requestOptimizedStorage = true;
    if (! params.is_null () && ! params->get ("Optimize Storage", true)) {
      requestOptimizedStorage = false;
    }
    if (getProfileType () == DynamicProfile) {
      // Pack 2-D storage (DynamicProfile) into 1-D packed storage.
      //
      // DynamicProfile means that the graph's column indices are
      // currently stored in a 2-D "unpacked" format, in the
      // arrays-of-arrays lclInds2D_.  We allocate 1-D storage
      // (ind_d) and then copy from 2-D storage (lclInds2D_) into 1-D
      // storage (ind_d).
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (k_numRowEntries_.dimension_0 ()) != lclNumRows,
         std::logic_error, "(DynamicProfile branch) "
         "k_numRowEntries_.dimension_0() = "
         << k_numRowEntries_.dimension_0 ()
         << " != getNodeNumRows() = " << lclNumRows << "");
#endif // HAVE_TPETRA_DEBUG

      // Pack the row offsets into ptr_d, by doing a sum-scan of the
      // array of valid entry counts per row (k_numRowEntries_).  The
      // pack method can handle its counts input being a host View.
      //
      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      {
        // Allocate the packed row offsets array.
        ptr_d = non_const_row_map_type ("Tpetra::CrsGraph::ptr", lclNumRows+1);
        typename row_entries_type::const_type numRowEnt_h = k_numRowEntries_;
        // This function can handle that numRowEnt_h lives on host.
        lclTotalNumEntries = computeOffsetsFromCounts (ptr_d, numRowEnt_h);
        ptr_d_const = ptr_d;
      }

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (ptr_d.dimension_0 ()) != lclNumRows + 1,
         std::logic_error, "(DynamicProfile branch) After packing ptr_d, "
         "ptr_d.dimension_0() = " << ptr_d.dimension_0 () << " != "
         "(lclNumRows+1) = " << (lclNumRows+1) << ".");
      {
        // mfh 26 Jun 2016: Don't assume UVM.  ptr_d lives in device
        // memory, so if we want to check its last entry on host, copy
        // it back to host explicitly.
        auto ptr_d_ent_d = Kokkos::subview (ptr_d, lclNumRows);
        auto ptr_d_ent_h = create_mirror_view (ptr_d_ent_d);
        Kokkos::deep_copy (ptr_d_ent_h, ptr_d_ent_d);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (ptr_d_ent_h () != lclTotalNumEntries, std::logic_error,
           "(DynamicProfile branch) After packing ptr_d, ptr_d(lclNumRows = "
           << lclNumRows << ") = " << ptr_d_ent_h () << " != total number of "
           "entries on the calling process = " << lclTotalNumEntries << ".");
      }
#endif // HAVE_TPETRA_DEBUG

      // Allocate the array of packed column indices.
      ind_d = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);
      // Pack the column indices.  We have to do this sequentially on
      // host, since lclInds2D_ is an ArrayRCP<Array<LO>>, which
      // doesn't work in parallel kernels (its iterators aren't even
      // thread safe in debug mode).
      {
        auto ptr_h = create_mirror_view (ptr_d);
        Kokkos::deep_copy (ptr_h, ptr_d); // we need the entries on host
        auto ind_h = create_mirror_view (ind_d); // we'll fill this on host

        // k_numRowEntries_ is a host View already, so we can use it here.
        typename row_entries_type::const_type numRowEnt_h = k_numRowEntries_;
        for (size_t row = 0; row < lclNumRows; ++row) {
          const size_t numEnt = numRowEnt_h(row);
          std::copy (lclInds2D_[row].begin (),
                     lclInds2D_[row].begin () + numEnt,
                     ind_h.ptr_on_device () + ptr_h(row));
        }
        Kokkos::deep_copy (ind_d, ind_h);
      }

#ifdef HAVE_TPETRA_DEBUG
      // Sanity check of packed row offsets.
      if (ptr_d.dimension_0 () != 0) {
        const size_t numOffsets = static_cast<size_t> (ptr_d.dimension_0 ());

        // mfh 26 Jun 2016: Don't assume UVM.  ptr_d lives in device
        // memory, so if we want to check its last entry on host, copy
        // it back to host explicitly.
        auto ptr_d_ent_d = Kokkos::subview (ptr_d, numOffsets - 1);
        auto ptr_d_ent_h = create_mirror_view (ptr_d_ent_d);
        Kokkos::deep_copy (ptr_d_ent_h, ptr_d_ent_d);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (ptr_d_ent_h ()) != ind_d.dimension_0 (),
           std::logic_error, "(DynamicProfile branch) After packing column "
           "indices, ptr_d(" << (numOffsets-1) << ") = " << ptr_d_ent_h ()
           << " != ind_d.dimension_0() = " << ind_d.dimension_0 () << ".");
      }
#endif // HAVE_TPETRA_DEBUG
    }
    else if (getProfileType () == StaticProfile) {
      // StaticProfile means that the graph's column indices are
      // currently stored in a 1-D format, with row offsets in
      // k_rowPtrs_ and local column indices in k_lclInds1D_.

#ifdef HAVE_TPETRA_DEBUG
      // StaticProfile also means that the graph's array of row
      // offsets must already be allocated.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (k_rowPtrs_.dimension_0 () == 0, std::logic_error,
         "(StaticProfile branch) k_rowPtrs_ has size zero, but shouldn't");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (k_rowPtrs_.dimension_0 () != lclNumRows + 1, std::logic_error,
         "(StaticProfile branch) k_rowPtrs_.dimension_0() = "
         << k_rowPtrs_.dimension_0 () << " != (lclNumRows + 1) = "
         << (lclNumRows + 1) << ".");
      {
        const size_t numOffsets = k_rowPtrs_.dimension_0 ();

        // mfh 26 Jun 2016: Don't assume UVM.  k_rowPtrs_ lives in
        // device memory, so if we want to check its last entry on
        // host, copy it back to host explicitly.
        auto k_rowPtrs_ent_d = Kokkos::subview (k_rowPtrs_, numOffsets - 1);
        auto k_rowPtrs_ent_h = create_mirror_view (k_rowPtrs_ent_d);
        Kokkos::deep_copy (k_rowPtrs_ent_h, k_rowPtrs_ent_d);

        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numOffsets != 0 &&
           k_lclInds1D_.dimension_0 () != k_rowPtrs_ent_h (),
           std::logic_error, "(StaticProfile branch) numOffsets = " <<
           numOffsets << " != 0 and k_lclInds1D_.dimension_0() = " <<
           k_lclInds1D_.dimension_0 () << " != k_rowPtrs_(" << numOffsets <<
           ") = " << k_rowPtrs_ent_h () << ".");
      }
#endif // HAVE_TPETRA_DEBUG

      if (nodeNumEntries_ != nodeNumAllocated_) {
        // The graph's current 1-D storage is "unpacked."  This means
        // the row offsets may differ from what the final row offsets
        // should be.  This could happen, for example, if the user
        // specified StaticProfile in the constructor and set an upper
        // bound on the number of entries in each row, but didn't fill
        // all those entries.

#ifdef HAVE_TPETRA_DEBUG
        if (k_rowPtrs_.dimension_0 () != 0) {
          const size_t numOffsets =
            static_cast<size_t> (k_rowPtrs_.dimension_0 ());
          // mfh 26 Jun 2016: Don't assume UVM.  k_rowPtrs_ lives in
          // device memory, so if we want to check its last entry on
          // host, copy it back to host explicitly.
          auto k_rowPtrs_ent_d = Kokkos::subview (k_rowPtrs_, numOffsets - 1);
          auto k_rowPtrs_ent_h = create_mirror_view (k_rowPtrs_ent_d);
          Kokkos::deep_copy (k_rowPtrs_ent_h, k_rowPtrs_ent_d);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (k_rowPtrs_ent_h () != static_cast<size_t> (k_lclInds1D_.dimension_0 ()),
             std::logic_error, "(StaticProfile unpacked branch) Before "
             "allocating or packing, k_rowPtrs_(" << (numOffsets-1) << ") = "
             << k_rowPtrs_ent_h () << " != k_lclInds1D_.dimension_0() = "
             << k_lclInds1D_.dimension_0 () << ".");
        }
#endif // HAVE_TPETRA_DEBUG

        // Pack the row offsets into ptr_d, by doing a sum-scan of the
        // array of valid entry counts per row (k_numRowEntries_).

        // Total number of entries in the matrix on the calling
        // process.  We will compute this in the loop below.  It's
        // cheap to compute and useful as a sanity check.
        size_t lclTotalNumEntries = 0;
        {
          // Allocate the packed row offsets array.
          ptr_d = non_const_row_map_type ("Tpetra::CrsGraph::ptr", lclNumRows + 1);
          ptr_d_const = ptr_d;

          // It's ok that k_numRowEntries_ is a host View; the
          // function can handle this.
          typename row_entries_type::const_type numRowEnt_h = k_numRowEntries_;
#ifdef HAVE_TPETRA_DEBUG
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (numRowEnt_h.dimension_0 ()) != lclNumRows,
             std::logic_error, "(StaticProfile unpacked branch) "
             "numRowEnt_h.dimension_0() = " << numRowEnt_h.dimension_0 ()
             << " != getNodeNumRows() = " << lclNumRows << "");
#endif // HAVE_TPETRA_DEBUG

          lclTotalNumEntries = computeOffsetsFromCounts (ptr_d, numRowEnt_h);

#ifdef HAVE_TPETRA_DEBUG
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (ptr_d.dimension_0 ()) != lclNumRows + 1,
             std::logic_error, "(StaticProfile unpacked branch) After "
             "allocating ptr_d, ptr_d.dimension_0() = " << ptr_d.dimension_0 ()
             << " != lclNumRows+1 = " << (lclNumRows+1) << ".");
          {
            // mfh 26 Jun 2016: Don't assume UVM.  ptr_d lives in device
            // memory, so if we want to check its last entry on host, copy
            // it back to host explicitly.
            auto ptr_d_ent_d = Kokkos::subview (ptr_d, lclNumRows);
            auto ptr_d_ent_h = create_mirror_view (ptr_d_ent_d);
            Kokkos::deep_copy (ptr_d_ent_h, ptr_d_ent_d);
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
              (ptr_d_ent_h () != lclTotalNumEntries, std::logic_error,
               "Tpetra::CrsGraph::fillLocalGraph: In StaticProfile unpacked "
               "branch, after filling ptr_d, ptr_d(lclNumRows=" << lclNumRows
               << ") = " << ptr_d_ent_h () << " != total number of entries on "
               "the calling process = " << lclTotalNumEntries << ".");
          }
#endif // HAVE_TPETRA_DEBUG
        }

        // Allocate the array of packed column indices.
        ind_d = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);

        // k_rowPtrs_ and k_lclInds1D_ are currently unpacked.  Pack
        // them, using the packed row offsets array ptr_d that we
        // created above.
        //
        // FIXME (mfh 08 Aug 2014) If "Optimize Storage" is false (in
        // CrsMatrix?), we need to keep around the unpacked row
        // offsets and column indices.

        // Pack the column indices from unpacked k_lclInds1D_ into
        // packed ind_d.  We will replace k_lclInds1D_ below.
        typedef pack_functor<
          typename local_graph_type::entries_type::non_const_type,
          row_map_type> inds_packer_type;
        inds_packer_type f (ind_d, k_lclInds1D_, ptr_d, k_rowPtrs_);
        Kokkos::parallel_for (lclNumRows, f);

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (ptr_d.dimension_0 () == 0, std::logic_error, "(StaticProfile "
           "\"Optimize Storage\"=true branch) After packing, "
           "ptr_d.dimension_0() = 0.  This probably means k_rowPtrs_ was "
           "never allocated.");
        if (ptr_d.dimension_0 () != 0) {
          const size_t numOffsets = static_cast<size_t> (ptr_d.dimension_0 ());
          // mfh 26 Jun 2016: Don't assume UVM.  ptr_d lives in device
          // memory, so if we want to check its last entry on host, copy
          // it back to host explicitly.
          auto ptr_d_ent_d = Kokkos::subview (ptr_d, numOffsets - 1);
          auto ptr_d_ent_h = create_mirror_view (ptr_d_ent_d);
          Kokkos::deep_copy (ptr_d_ent_h, ptr_d_ent_d);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (ptr_d_ent_h ()) != ind_d.dimension_0 (),
             std::logic_error, "(StaticProfile \"Optimize Storage\"=true "
             "branch) After packing, ptr_d(" << (numOffsets-1) << ") = " <<
             ptr_d_ent_h () << " != ind_d.dimension_0() = " <<
             ind_d.dimension_0 () << ".");
        }
#endif // HAVE_TPETRA_DEBUG
      }
      else { // We don't have to pack, so just set the pointers.
        ptr_d_const = k_rowPtrs_;
        ind_d = k_lclInds1D_;

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (ptr_d_const.dimension_0 () == 0, std::logic_error, "(StaticProfile "
           "\"Optimize Storage\"=false branch) ptr_d_const.dimension_0() = 0.  "
           "This probably means that k_rowPtrs_ was never allocated.");
        if (ptr_d_const.dimension_0 () != 0) {
          const size_t numOffsets =
            static_cast<size_t> (ptr_d_const.dimension_0 ());
          // mfh 26 Jun 2016: Don't assume UVM.  ptr_d_const lives in
          // device memory, so if we want to check its last entry on
          // host, copy it back to host explicitly.
          auto ptr_d_const_ent_d = Kokkos::subview (ptr_d_const, numOffsets - 1);
          auto ptr_d_const_ent_h = create_mirror_view (ptr_d_const_ent_d);
          Kokkos::deep_copy (ptr_d_const_ent_h, ptr_d_const_ent_d);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (ptr_d_const_ent_h ()) != ind_d.dimension_0 (),
             std::logic_error, "(StaticProfile \"Optimize Storage\"=false "
             "branch) ptr_d_const(" << (numOffsets-1) << ") = " <<
             ptr_d_const_ent_h () << " != ind_d.dimension_0() = " <<
             ind_d.dimension_0 () << ".");
        }
#endif // HAVE_TPETRA_DEBUG
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    // Extra sanity checks.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (ptr_d_const.dimension_0 ()) != lclNumRows + 1,
       std::logic_error, "After packing, ptr_d_const.dimension_0() = " <<
       ptr_d_const.dimension_0 () << " != lclNumRows+1 = " << (lclNumRows+1)
       << ".");
    if (ptr_d_const.dimension_0 () != 0) {
      const size_t numOffsets = static_cast<size_t> (ptr_d_const.dimension_0 ());
      // mfh 26 Jun 2016: Don't assume UVM.  ptr_d_const lives in
      // device memory, so if we want to check its last entry on
      // host, copy it back to host explicitly.
      auto ptr_d_const_ent_d = Kokkos::subview (ptr_d_const, numOffsets - 1);
      auto ptr_d_const_ent_h = create_mirror_view (ptr_d_const_ent_d);
      Kokkos::deep_copy (ptr_d_const_ent_h, ptr_d_const_ent_d);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (ptr_d_const_ent_h ()) != ind_d.dimension_0 (),
         std::logic_error, "After packing, ptr_d_const(" << (numOffsets-1) <<
         ") = " << ptr_d_const_ent_h () << " != ind_d.dimension_0() = " <<
         ind_d.dimension_0 () << ".");
    }
#endif // HAVE_TPETRA_DEBUG

    if (requestOptimizedStorage) {
      // With optimized storage, we don't need to store the 2-D column
      // indices array-of-arrays, or the array of row entry counts.

      // Free graph data structures that are only needed for 2-D or
      // unpacked 1-D storage.
      lclInds2D_ = Teuchos::null;
      k_numRowEntries_ = row_entries_type ();

      // Keep the new 1-D packed allocations.
      k_rowPtrs_   = ptr_d_const;
      k_lclInds1D_ = ind_d;

      // Storage is packed now, so the number of allocated entries is
      // the same as the actual number of entries.
      nodeNumAllocated_ = nodeNumEntries_;
      // The graph is definitely StaticProfile now, whether or not it
      // was before.
      pftype_ = StaticProfile;
    }

    // FIXME (mfh 28 Aug 2014) "Local Graph" sublist no longer used.

    // Build the local graph.
    lclGraph_ = local_graph_type (ind_d, ptr_d_const);

    // TODO (mfh 13 Mar 2014) getNodeNumDiags(), isUpperTriangular(),
    // and isLowerTriangular() depend on computeGlobalConstants(), in
    // particular the part where it looks at the local matrix.  You
    // have to use global indices to determine which entries are
    // diagonal, or above or below the diagonal.  However, lower or
    // upper triangularness is a local property.
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  replaceColMap (const Teuchos::RCP<const map_type>& newColMap)
  {
    // NOTE: This safety check matches the code, but not the documentation of Crsgraph
    //
    // FIXME (mfh 18 Aug 2014) This will break if the calling process
    // has no entries, because in that case, currently it is neither
    // locally nor globally indexed.  This will change once we get rid
    // of lazy allocation (so that the constructor allocates indices
    // and therefore commits to local vs. global).
    const char tfecfFuncName[] = "replaceColMap: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed () || isGloballyIndexed (), std::runtime_error,
      "Requires matching maps and non-static graph.");
    colMap_ = newColMap;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  reindexColumns (const Teuchos::RCP<const map_type>& newColMap,
                  const Teuchos::RCP<const import_type>& newImport,
                  const bool sortIndicesInEachRow)
  {
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
    typedef typename local_graph_type::entries_type::non_const_type col_inds_type;
    const char tfecfFuncName[] = "reindexColumns: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete (), std::runtime_error, "The graph is fill complete "
      "(isFillComplete() returns true).  You must call resumeFill() before "
      "you may call this method.");

    // mfh 19 Aug 2014: This method does NOT redistribute data; it
    // doesn't claim to do the work of an Import or Export.  This
    // means that for all processes, the calling process MUST own all
    // column indices, in both the old column Map (if it exists) and
    // the new column Map.  We check this via an all-reduce.
    //
    // Some processes may be globally indexed, others may be locally
    // indexed, and others (that have no graph entries) may be
    // neither.  This method will NOT change the graph's current
    // state.  If it's locally indexed, it will stay that way, and
    // vice versa.  It would easy to add an option to convert indices
    // from global to local, so as to save a global-to-local
    // conversion pass.  However, we don't do this here.  The intended
    // typical use case is that the graph already has a column Map and
    // is locally indexed, and this is the case for which we optimize.

    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());

    // Attempt to convert indices to the new column Map's version of
    // local.  This will fail if on the calling process, the graph has
    // indices that are not on that process in the new column Map.
    // After the local conversion attempt, we will do an all-reduce to
    // see if any processes failed.

    // If this is false, then either the graph contains a column index
    // which is invalid in the CURRENT column Map, or the graph is
    // locally indexed but currently has no column Map.  In either
    // case, there is no way to convert the current local indices into
    // global indices, so that we can convert them into the new column
    // Map's local indices.  It's possible for this to be true on some
    // processes but not others, due to replaceColMap.
    bool allCurColIndsValid = true;
    // On the calling process, are all valid current column indices
    // also in the new column Map on the calling process?  In other
    // words, does local reindexing suffice, or should the user have
    // done an Import or Export instead?
    bool localSuffices = true;

    // Final arrays for the local indices.  We will allocate exactly
    // one of these ONLY if the graph is locally indexed on the
    // calling process, and ONLY if the graph has one or more entries
    // (is not empty) on the calling process.  In that case, we
    // allocate the first (1-D storage) if the graph has a static
    // profile, else we allocate the second (2-D storage).
    typename local_graph_type::entries_type::non_const_type newLclInds1D;
    Teuchos::ArrayRCP<Teuchos::Array<LO> > newLclInds2D;

    // If indices aren't allocated, that means the calling process
    // owns no entries in the graph.  Thus, there is nothing to
    // convert, and it trivially succeeds locally.
    if (indicesAreAllocated ()) {
      if (isLocallyIndexed ()) {
        if (hasColMap ()) { // locally indexed, and currently has a column Map
          const map_type& oldColMap = * (getColMap ());
          if (pftype_ == StaticProfile) {
            // Allocate storage for the new local indices.
            RCP<node_type> node = getRowMap ()->getNode ();
            newLclInds1D =
              col_inds_type ("Tpetra::CrsGraph::ind", nodeNumAllocated_);
            // Attempt to convert the new indices locally.
            for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
              const RowInfo rowInfo = this->getRowInfo (lclRow);
              const size_t beg = rowInfo.offset1D;
              const size_t end = beg + rowInfo.numEntries;
              for (size_t k = beg; k < end; ++k) {
                // FIXME (mfh 21 Aug 2014) This assumes UVM.  Should
                // use a DualView instead.
                const LO oldLclCol = k_lclInds1D_(k);
                if (oldLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                  allCurColIndsValid = false;
                  break; // Stop at the first invalid index
                }
                const GO gblCol = oldColMap.getGlobalElement (oldLclCol);

                // The above conversion MUST succeed.  Otherwise, the
                // current local index is invalid, which means that
                // the graph was constructed incorrectly.
                if (gblCol == Teuchos::OrdinalTraits<GO>::invalid ()) {
                  allCurColIndsValid = false;
                  break; // Stop at the first invalid index
                }
                else {
                  const LO newLclCol = newColMap->getLocalElement (gblCol);
                  if (newLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                    localSuffices = false;
                    break; // Stop at the first invalid index
                  }
                  // FIXME (mfh 21 Aug 2014) This assumes UVM.  Should
                  // use a DualView instead.
                  newLclInds1D(k) = newLclCol;
                }
              } // for each entry in the current row
            } // for each locally owned row
          }
          else { // pftype_ == DynamicProfile
            // Allocate storage for the new local indices.  We only
            // allocate the outer array here; we will allocate the
            // inner arrays below.
            newLclInds2D = Teuchos::arcp<Teuchos::Array<LO> > (lclNumRows);

            // Attempt to convert the new indices locally.
            for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
              const RowInfo rowInfo = this->getRowInfo (lclRow);
              newLclInds2D.resize (rowInfo.allocSize);

              Teuchos::ArrayView<const LO> oldLclRowView = getLocalView (rowInfo);
              Teuchos::ArrayView<LO> newLclRowView = (newLclInds2D[lclRow]) ();

              for (size_t k = 0; k < rowInfo.numEntries; ++k) {
                const LO oldLclCol = oldLclRowView[k];
                if (oldLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                  allCurColIndsValid = false;
                  break; // Stop at the first invalid index
                }
                const GO gblCol = oldColMap.getGlobalElement (oldLclCol);

                // The above conversion MUST succeed.  Otherwise, the
                // local index is invalid and the graph is wrong.
                if (gblCol == Teuchos::OrdinalTraits<GO>::invalid ()) {
                  allCurColIndsValid = false;
                  break; // Stop at the first invalid index
                }
                else {
                  const LO newLclCol = newColMap->getLocalElement (gblCol);
                  if (newLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                    localSuffices = false;
                    break; // Stop at the first invalid index.
                  }
                  newLclRowView[k] = newLclCol;
                }
              } // for each entry in the current row
            } // for each locally owned row
          } // pftype_
        }
        else { // locally indexed, but no column Map
          // This case is only possible if replaceColMap() was called
          // with a null argument on the calling process.  It's
          // possible, but it means that this method can't possibly
          // succeed, since we have no way of knowing how to convert
          // the current local indices to global indices.
          allCurColIndsValid = false;
        }
      }
      else { // globally indexed
        // If the graph is globally indexed, we don't need to save
        // local indices, but we _do_ need to know whether the current
        // global indices are valid in the new column Map.  We may
        // need to do a getRemoteIndexList call to find this out.
        //
        // In this case, it doesn't matter whether the graph currently
        // has a column Map.  We don't need the old column Map to
        // convert from global indices to the _new_ column Map's local
        // indices.  Furthermore, we can use the same code, whether
        // the graph is static or dynamic profile.

        // Test whether the current global indices are in the new
        // column Map on the calling process.
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const RowInfo rowInfo = this->getRowInfo (lclRow);
          Teuchos::ArrayView<const GO> oldGblRowView = getGlobalView (rowInfo);
          for (size_t k = 0; k < rowInfo.numEntries; ++k) {
            const GO gblCol = oldGblRowView[k];
            if (! newColMap->isNodeGlobalElement (gblCol)) {
              localSuffices = false;
              break; // Stop at the first invalid index
            }
          } // for each entry in the current row
        } // for each locally owned row
      } // locally or globally indexed
    } // whether indices are allocated

    // Do an all-reduce to check both possible error conditions.
    int lclSuccess[2];
    lclSuccess[0] = allCurColIndsValid ? 1 : 0;
    lclSuccess[1] = localSuffices ? 1 : 0;
    int gblSuccess[2];
    gblSuccess[0] = 0;
    gblSuccess[1] = 0;
    RCP<const Teuchos::Comm<int> > comm =
      getRowMap ().is_null () ? Teuchos::null : getRowMap ()->getComm ();
    if (! comm.is_null ()) {
      reduceAll<int, int> (*comm, REDUCE_MIN, 2, lclSuccess, gblSuccess);
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      gblSuccess[0] == 0, std::runtime_error, "It is not possible to continue."
      "  The most likely reason is that the graph is locally indexed, but the "
      "column Map is missing (null) on some processes, due to a previous call "
      "to replaceColMap().");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      gblSuccess[1] == 0, std::runtime_error, "On some process, the graph "
      "contains column indices that are in the old column Map, but not in the "
      "new column Map (on that process).  This method does NOT redistribute "
      "data; it does not claim to do the work of an Import or Export operation."
      "  This means that for all processess, the calling process MUST own all "
      "column indices, in both the old column Map and the new column Map.  In "
      "this case, you will need to do an Import or Export operation to "
      "redistribute data.");

    // Commit the results.
    if (isLocallyIndexed ()) {
      if (pftype_ == StaticProfile) {
        k_lclInds1D_ = newLclInds1D;
      } else { // dynamic profile
        lclInds2D_ = newLclInds2D;
      }
      // We've reindexed, so we don't know if the indices are sorted.
      //
      // FIXME (mfh 17 Sep 2014) It could make sense to check this,
      // since we're already going through all the indices above.  We
      // could also sort each row in place; that way, we would only
      // have to make one pass over the rows.
      indicesAreSorted_ = false;
      if (sortIndicesInEachRow) {
        // NOTE (mfh 17 Sep 2014) The graph must be locally indexed in
        // order to call this method.
        //
        // FIXME (mfh 17 Sep 2014) This violates the strong exception
        // guarantee.  It would be better to sort the new index arrays
        // before committing them.
        sortAllIndices ();
      }
    }
    colMap_ = newColMap;

    if (newImport.is_null ()) {
      // FIXME (mfh 19 Aug 2014) Should use the above all-reduce to
      // check whether the input Import is null on any process.
      //
      // If the domain Map hasn't been set yet, we can't compute a new
      // Import object.  Leave it what it is; it should be null, but
      // it doesn't matter.  If the domain Map _has_ been set, then
      // compute a new Import object if necessary.
      if (! domainMap_.is_null ()) {
        if (! domainMap_->isSameAs (* newColMap)) {
          importer_ = Teuchos::rcp (new import_type (domainMap_, newColMap));
        } else {
          importer_ = Teuchos::null; // don't need an Import
        }
      }
    } else {
      // The caller gave us an Import object.  Assume that it's valid.
      importer_ = newImport;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                               const Teuchos::RCP<const import_type>& newImporter)
  {
    const char prefix[] = "Tpetra::CrsGraph::replaceDomainMapAndImporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION(
      colMap_.is_null (), std::invalid_argument, prefix << "You may not call "
      "this method unless the graph already has a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      newDomainMap.is_null (), std::invalid_argument,
      prefix << "The new domain Map must be nonnull.");

#ifdef HAVE_TPETRA_DEBUG
    if (newImporter.is_null ()) {
      // It's not a good idea to put expensive operations in a macro
      // clause, even if they are side effect - free, because macros
      // don't promise that they won't evaluate their arguments more
      // than once.  It's polite for them to do so, but not required.
      const bool colSameAsDom = colMap_->isSameAs (*newDomainMap);
      TEUCHOS_TEST_FOR_EXCEPTION(
        colSameAsDom, std::invalid_argument, "If the new Import is null, "
        "then the new domain Map must be the same as the current column Map.");
    }
    else {
      const bool colSameAsTgt =
        colMap_->isSameAs (* (newImporter->getTargetMap ()));
      const bool newDomSameAsSrc =
        newDomainMap->isSameAs (* (newImporter->getSourceMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! colSameAsTgt || ! newDomSameAsSrc, std::invalid_argument, "If the "
        "new Import is nonnull, then the current column Map must be the same "
        "as the new Import's target Map, and the new domain Map must be the "
        "same as the new Import's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG

    domainMap_ = newDomainMap;
    importer_ = Teuchos::rcp_const_cast<import_type> (newImporter);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::local_graph_type
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalGraph () const
  {
    return lclGraph_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  computeGlobalConstants ()
  {
    using Teuchos::outArg;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(! hasColMap(), std::logic_error, "Tpetra::"
      "CrsGraph::computeGlobalConstants: At this point, the graph should have "
      "a column Map, but it does not.  Please report this bug to the Tpetra "
      "developers.");
#endif // HAVE_TPETRA_DEBUG

    // If necessary, (re)compute the local constants: nodeNumDiags_,
    // lowerTriangular_, upperTriangular_, and nodeMaxNumRowEntries_.
    if (! haveLocalConstants_) {
      // We have actually already computed nodeNumEntries_.
      // nodeNumEntries_ gets updated by methods that insert or remove
      // indices (including setAllIndices and
      // expertStaticFillComplete).  Before fillComplete, its count
      // may include duplicate column indices in the same row.
      // However, mergeRowIndices and mergeRowIndicesAndValues both
      // subtract off merged indices in each row from the total count.
      // Thus, nodeNumEntries_ _should_ be accurate at this point,
      // meaning that we don't have to re-count it here.

      // Reset local properties
      upperTriangular_ = true;
      lowerTriangular_ = true;
      nodeMaxNumRowEntries_ = 0;
      nodeNumDiags_         = 0;

      // At this point, we know that we have both a row Map and a column Map.
      const map_type& rowMap = *rowMap_;
      const map_type& colMap = *colMap_;

      // Go through all the entries of the graph.  Count the number of
      // diagonal elements we encounter, and figure out whether the
      // graph is lower or upper triangular.  Diagonal elements are
      // determined using global indices, with respect to the whole
      // graph.  However, lower or upper triangularity is a local
      // property, and is determined using local indices.
      //
      // At this point, indices have already been sorted in each row.
      // That makes finding out whether the graph is lower / upper
      // triangular easier.
      if (indicesAreAllocated () && nodeNumAllocated_ > 0) {
        const LO numLocalRows = static_cast<LO> (this->getNodeNumRows ());
        for (LO localRow = 0; localRow < numLocalRows; ++localRow) {
          const GO globalRow = rowMap.getGlobalElement (localRow);
          // Find the local (column) index for the diagonal entry.
          // This process might not necessarily own _any_ entries in
          // the current row.  If it doesn't, skip this row.  It won't
          // affect any of the attributes (nodeNumDiagons_,
          // upperTriangular_, lowerTriangular_, or
          // nodeMaxNumRowEntries_) which this loop sets.
          const LO rlcid = colMap.getLocalElement (globalRow);
            // This process owns one or more entries in the current row.
            const RowInfo rowInfo = this->getRowInfo (localRow);
            ArrayView<const LO> rview = getLocalView (rowInfo);
            typename ArrayView<const LO>::iterator beg, end, cur;
            beg = rview.begin();
            end = beg + rowInfo.numEntries;
            if (beg != end) {
              for (cur = beg; cur != end; ++cur) {
                // is this the diagonal?
                if (rlcid == *cur) ++nodeNumDiags_;
              }
              // Local column indices are sorted in each row.  That means
              // the smallest column index in this row (on this process)
              // is *beg, and the largest column index in this row (on
              // this process) is *(end - 1).  We know that end - 1 is
              // valid because beg != end.
              const size_t smallestCol = static_cast<size_t> (*beg);
              const size_t largestCol = static_cast<size_t> (*(end - 1));

              if (smallestCol < static_cast<size_t> (localRow)) {
                upperTriangular_ = false;
              }
              if (static_cast<size_t> (localRow) < largestCol) {
                lowerTriangular_ = false;
              }
            }
            // Update the max number of entries over all rows.
            nodeMaxNumRowEntries_ = std::max (nodeMaxNumRowEntries_, rowInfo.numEntries);
        }
      }
      haveLocalConstants_ = true;
    } // if my process doesn't have local constants

    // Compute global constants from local constants.  Processes that
    // already have local constants still participate in the
    // all-reduces, using their previously computed values.
    if (haveGlobalConstants_ == false) {
      // Promote all the nodeNum* and nodeMaxNum* quantities from
      // size_t to global_size_t, when doing the all-reduces for
      // globalNum* / globalMaxNum* results.
      //
      // FIXME (mfh 07 May 2013) Unfortunately, we either have to do
      // this in two all-reduces (one for the sum and the other for
      // the max), or use a custom MPI_Op that combines the sum and
      // the max.  The latter might even be slower than two
      // all-reduces on modern network hardware.  It would also be a
      // good idea to use nonblocking all-reduces (MPI 3), so that we
      // don't have to wait around for the first one to finish before
      // starting the second one.
      GST lcl[2], gbl[2];
      lcl[0] = static_cast<GST> (nodeNumEntries_);
      lcl[1] = static_cast<GST> (nodeNumDiags_);
      reduceAll<int,GST> (*getComm (), Teuchos::REDUCE_SUM,
                          2, lcl, gbl);
      globalNumEntries_ = gbl[0];
      globalNumDiags_   = gbl[1];
      reduceAll<int,GST> (*getComm (), Teuchos::REDUCE_MAX,
                          static_cast<GST> (nodeMaxNumRowEntries_),
                          outArg (globalMaxNumRowEntries_));
      haveGlobalConstants_ = true;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  makeIndicesLocal ()
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename local_graph_type::entries_type::non_const_type
      lcl_col_inds_type;
    typedef Kokkos::View<GO*, typename lcl_col_inds_type::array_layout,
      device_type> gbl_col_inds_type;
    typedef decltype (k_numRowEntries_) row_entries_type;
    const char tfecfFuncName[] = "makeIndicesLocal: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::logic_error, "The graph does not have a "
       "column Map yet.  This method should never be called in that case.  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->getColMap ().is_null (), std::logic_error, "The graph claims "
       "that it has a column Map, because hasColMap() returns true.  However, "
       "the result of getColMap() is null.  This should never happen.  Please "
       "report this bug to the Tpetra developers.");

    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
    const map_type& colMap = * (this->getColMap ());

    if (isGloballyIndexed () && lclNumRows != 0) {
      // This is a host View.
      typename row_entries_type::const_type h_numRowEnt = k_numRowEntries_;

      // allocate data for local indices
      if (getProfileType () == StaticProfile) {
        // If GO and LO are the same size, we can reuse the existing
        // array of 1-D index storage to convert column indices from
        // GO to LO.  Otherwise, we'll just allocate a new buffer.
        constexpr bool LO_GO_same = std::is_same<LO, GO>::value;
        if (nodeNumAllocated_ && LO_GO_same) {
          k_lclInds1D_ = Kokkos::Impl::if_c<LO_GO_same,
            t_GlobalOrdinal_1D,
            lcl_col_inds_type>::select (k_gblInds1D_, k_lclInds1D_);
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (k_rowPtrs_.dimension_0 () == 0, std::logic_error,
             "k_rowPtrs_.dimension_0() == 0.  This should never happen at this "
             "point.  Please report this bug to the Tpetra developers.");

          // FIXME (mfh 26 Jun 2016) This assumes UVM.
          const size_t numEnt = k_rowPtrs_[lclNumRows];
          k_lclInds1D_ = lcl_col_inds_type ("Tpetra::CrsGraph::lclind", numEnt);
        }

        // FIXME (mfh 26 Jun 2016) Convert to a Kokkos parallel_reduce
        // functor (reduce over error code).  Use LocalMap for
        // global-to-local index conversions in the parallel loop.
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          // FIXME (mfh 26 Jun 2016) This assumes UVM (for k_rowPtrs_).
          const size_t offset = k_rowPtrs_(lclRow);
          // NOTE (mfh 26 Jun 2016) It's always legal to cast the
          // number of entries in a row to LO, as long as the row
          // doesn't have too many duplicate entries.
          const LO numEnt = static_cast<LO> (h_numRowEnt(lclRow));
          for (LO j = 0; j < numEnt; ++j) {
            // FIXME (mfh 26 Jun 2016) This assumes UVM (for k_gblInds1D_).
            const GO gid = k_gblInds1D_(offset + j);
            const LO lid = colMap.getLocalElement (gid);
            // FIXME (mfh 26 Jun 2016) This assumes UVM (for k_lclInds1D_).
            k_lclInds1D_(offset + j) = lid;
#ifdef HAVE_TPETRA_DEBUG
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
              (lid == Teuchos::OrdinalTraits<LO>::invalid(), std::logic_error,
               "In local row " << lclRow << ", global column " << gid << " is "
               "not in the column Map.  This should never happen.  Please "
               "report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
          }
        }
        // We've converted column indices from global to local, so we
        // can deallocate the global column indices (which we know are
        // in 1-D storage, because the graph has static profile).
        k_gblInds1D_ = gbl_col_inds_type ();
      }
      else {  // the graph has dynamic profile (2-D index storage)
        lclInds2D_ = Teuchos::arcp<Teuchos::Array<LO> > (lclNumRows);
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          if (! gblInds2D_[lclRow].empty ()) {
            const GO* const ginds = gblInds2D_[lclRow].getRawPtr ();
            // NOTE (mfh 26 Jun 2016) It's always legal to cast the
            // number of entries in a row to LO, as long as the row
            // doesn't have too many duplicate entries.
            const LO rna = static_cast<LO> (gblInds2D_[lclRow].size ());
            const LO numEnt = static_cast<LO> (h_numRowEnt(lclRow));
            lclInds2D_[lclRow].resize (rna);
            LO* const linds = lclInds2D_[lclRow].getRawPtr ();
            for (LO j = 0; j < numEnt; ++j) {
              const GO gid = ginds[j];
              const LO lid = colMap.getLocalElement (gid);
              linds[j] = lid;
#ifdef HAVE_TPETRA_DEBUG
              TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
                (linds[j] == Teuchos::OrdinalTraits<LO>::invalid(),
                 std::logic_error,
                 "Global column ginds[j=" << j << "]=" << ginds[j]
                 << " of local row " << lclRow << " is not in the column Map.  "
                 "This should never happen.  Please report this bug to the "
                 "Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
            }
          }
        }
        gblInds2D_ = Teuchos::null;
      }
    } // globallyIndexed() && lclNumRows > 0

    lclGraph_ = local_graph_type (k_lclInds1D_, k_rowPtrs_);
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  sortAllIndices ()
  {
    typedef LocalOrdinal LO;

    // this should be called only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed () );
    if (isSorted () == false) {
      // FIXME (mfh 06 Mar 2014) This would be a good place for a
      // thread-parallel kernel.
      const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const RowInfo rowInfo = this->getRowInfo (lclRow);
        this->sortRowIndices (rowInfo);
      }
    }
    indicesAreSorted_ = true; // we just sorted every row
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  makeColMap ()
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "makeColMap";

    if (hasColMap ()) { // The graph already has a column Map.
      // FIXME (mfh 26 Feb 2013): This currently prevents structure
      // changes that affect the column Map.
      return;
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      ": The graph is locally indexed.  Calling makeColMap() to make the "
      "column Map requires that the graph be globally indexed.");

    // After the calling process is done going through all of the rows
    // it owns, myColumns will contain the list of indices owned by
    // this process in the column Map.
    Array<GO> myColumns;

    // If we reach this point, we don't have a column Map yet, so the
    // graph can't be locally indexed.  Thus, isGloballyIndexed() ==
    // false means that the graph is empty on this process, so
    // myColumns will be left empty.
    if (isGloballyIndexed ()) {
      // Go through all the rows, finding the populated column indices.
      //
      // Our final list of indices for the column Map constructor will
      // have the following properties (all of which are with respect
      // to the calling process):
      //
      // 1. Indices in the domain Map go first.
      // 2. Indices not in the domain Map follow, ordered first
      //    contiguously by their owning process rank (in the domain
      //    Map), then in increasing order within that.
      // 3. No duplicate indices.
      //
      // This imitates the ordering used by Aztec(OO) and Epetra.
      // Storing indices owned by the same process (in the domain Map)
      // contiguously permits the use of contiguous send and receive
      // buffers.
      //
      // We begin by partitioning the column indices into "local" GIDs
      // (owned by the domain Map) and "remote" GIDs (not owned by the
      // domain Map).  We use the same order for local GIDs as the
      // domain Map, so we track them in place in their array.  We use
      // an std::set (RemoteGIDSet) to keep track of remote GIDs, so
      // that we don't have to merge duplicates later.
      const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
      size_t numLocalColGIDs = 0, numRemoteColGIDs = 0;

      // GIDisLocal[lid] == 0 if and only if local index lid in the
      // domain Map is remote (not local).
      Array<char> GIDisLocal (domainMap_->getNodeNumElements (), 0);
      std::set<GO> RemoteGIDSet;
      // This preserves the not-sorted Epetra order of GIDs.
      std::vector<GO> RemoteGIDUnorderedVector;
      const LO myNumRows = this->getNodeNumRows ();
      for (LO r = 0; r < myNumRows; ++r) {
        const RowInfo rowinfo = this->getRowInfo (r);
        if (rowinfo.numEntries > 0) {
          // NOTE (mfh 02 Sep 2014) getGlobalView() returns a view of
          // all the space in the row, not just the occupied entries.
          // (This matters for the case of unpacked 1-D storage.  We
          // might not have packed it yet.)  That's why we need to
          // take a subview.
          ArrayView<const GO> rowGids = getGlobalView (rowinfo);
          rowGids = rowGids (0, rowinfo.numEntries);

          for (size_t k = 0; k < rowinfo.numEntries; ++k) {
            const GO gid = rowGids[k];
            const LO lid = domainMap_->getLocalElement (gid);
            if (lid != LINV) {
              const char alreadyFound = GIDisLocal[lid];
              if (alreadyFound == 0) {
                GIDisLocal[lid] = static_cast<char> (1);
                ++numLocalColGIDs;
              }
            }
            else {
              const bool notAlreadyFound = RemoteGIDSet.insert (gid).second;
              if (notAlreadyFound) { // gid did not exist in the set before
                if (! sortGhostsAssociatedWithEachProcessor_) {
                  // The user doesn't want to sort remote GIDs (for
                  // each remote process); they want us to keep remote
                  // GIDs in their original order.  We do this by
                  // stuffing each remote GID into an array as we
                  // encounter it for the first time.  The std::set
                  // helpfully tracks first encounters.
                  RemoteGIDUnorderedVector.push_back (gid);
                }
                ++numRemoteColGIDs;
              }
            }
          } // for each entry k in row r
        } // if row r contains > 0 entries
      } // for each locally owned row r

      // Possible short-circuit for serial scenario:
      //
      // If all domain GIDs are present as column indices, then set
      // ColMap=DomainMap.  By construction, LocalGIDs is a subset of
      // DomainGIDs.
      //
      // If we have
      //   * Number of remote GIDs is 0, so that ColGIDs == LocalGIDs,
      // and
      //   * Number of local GIDs is number of domain GIDs
      // then
      //   * LocalGIDs \subset DomainGIDs && size(LocalGIDs) ==
      //     size(DomainGIDs) => DomainGIDs == LocalGIDs == ColGIDs
      // on the calling process.
      //
      // We will concern ourselves only with the special case of a
      // serial DomainMap, obviating the need for communication.
      //
      // If
      //   * DomainMap has a serial communicator
      // then we can set the column Map as the domain Map
      // return. Benefit: this graph won't need an Import object
      // later.
      //
      // Note, for a serial domain map, there can be no RemoteGIDs,
      // because there are no remote processes.  Likely explanations
      // for this are:
      //  * user submitted erroneous column indices
      //  * user submitted erroneous domain Map
      if (domainMap_->getComm ()->getSize () == 1) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          numRemoteColGIDs != 0, std::runtime_error,
          ": " << numRemoteColGIDs << " column "
          << (numRemoteColGIDs != 1 ? "indices are" : "index is")
          << " not in the domain Map." << endl
          << "Either these indices are invalid or the domain Map is invalid."
          << endl << "Remember that nonsquare matrices, or matrices where the "
          "row and range Maps are different, require calling the version of "
          "fillComplete that takes the domain and range Maps as input.");
        if (numLocalColGIDs == domainMap_->getNodeNumElements()) {
          colMap_ = domainMap_;
          checkInternalState ();
          return;
        }
      }

      // Populate myColumns with a list of all column GIDs.  Put
      // locally owned (in the domain Map) GIDs at the front: they
      // correspond to "same" and "permuted" entries between the
      // column Map and the domain Map.  Put remote GIDs at the back.
      myColumns.resize (numLocalColGIDs + numRemoteColGIDs);
      // get pointers into myColumns for each part
      ArrayView<GO> LocalColGIDs  = myColumns (0, numLocalColGIDs);
      ArrayView<GO> RemoteColGIDs = myColumns (numLocalColGIDs, numRemoteColGIDs);

      // Copy the remote GIDs into myColumns
      if (sortGhostsAssociatedWithEachProcessor_) {
        // The std::set puts GIDs in increasing order.
        std::copy (RemoteGIDSet.begin(), RemoteGIDSet.end(),
                   RemoteColGIDs.begin());
      } else {
        // Respect the originally encountered order.
        std::copy (RemoteGIDUnorderedVector.begin(),
                   RemoteGIDUnorderedVector.end(), RemoteColGIDs.begin());
      }

      // Make a list of process ranks corresponding to the remote GIDs.
      Array<int> RemoteImageIDs (numRemoteColGIDs);
      // Look up the remote process' ranks in the domain Map.
      {
        const LookupStatus stat =
          domainMap_->getRemoteIndexList (RemoteColGIDs, RemoteImageIDs ());
#ifdef HAVE_TPETRA_DEBUG
        // If any process returns IDNotPresent, then at least one of
        // the remote indices was not present in the domain Map.  This
        // means that the Import object cannot be constructed, because
        // of incongruity between the column Map and domain Map.
        // This has two likely causes:
        //   - The user has made a mistake in the column indices
        //   - The user has made a mistake with respect to the domain Map
        const int missingID_lcl = (stat == IDNotPresent ? 1 : 0);
        int missingID_gbl = 0;
        reduceAll<int, int> (*getComm (), REDUCE_MAX, missingID_lcl,
                             outArg (missingID_gbl));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          missingID_gbl == 1, std::runtime_error,
          ": Some column indices are not in the domain Map." << endl
          << "Either these column indices are invalid or the domain Map is "
          "invalid." << endl << "Likely cause: For a nonsquare matrix, you "
          "must give the domain and range Maps as input to fillComplete.");
#else
        (void) stat; // forestall compiler warning for unused variable
#endif // HAVE_TPETRA_DEBUG
      }
      // Sort incoming remote column indices by their owning process
      // rank, so that all columns coming from a given remote process
      // are contiguous.  This means the Import's Distributor doesn't
      // need to reorder data.
      //
      // NOTE (mfh 02 Sep 2014) This needs to be a stable sort, so
      // that it respects either of the possible orderings of GIDs
      // (sorted, or original order) specified above.
      sort2 (RemoteImageIDs.begin(), RemoteImageIDs.end(), RemoteColGIDs.begin());

      // Copy the local GIDs into myColumns. Two cases:
      // 1. If the number of Local column GIDs is the same as the
      //    number of Local domain GIDs, we can simply read the domain
      //    GIDs into the front part of ColIndices (see logic above
      //    from the serial short circuit case)
      // 2. We step through the GIDs of the DomainMap, checking to see
      //    if each domain GID is a column GID.  We want to do this to
      //    maintain a consistent ordering of GIDs between the columns
      //    and the domain.

      const size_t numDomainElts = domainMap_->getNodeNumElements ();
      if (numLocalColGIDs == numDomainElts) {
        // If the number of locally owned GIDs are the same as the
        // number of local domain Map elements, then the local domain
        // Map elements are the same as the locally owned GIDs.
        if (domainMap_->isContiguous ()) {
          // NOTE (mfh 03 Mar 2013, 02 Sep 2014) In the common case
          // that the domain Map is contiguous, it's more efficient to
          // avoid calling getNodeElementList(), since that
          // permanently constructs and caches the GID list in the
          // contiguous Map.
          GO curColMapGid = domainMap_->getMinGlobalIndex ();
          for (size_t k = 0; k < numLocalColGIDs; ++k, ++curColMapGid) {
            LocalColGIDs[k] = curColMapGid;
          }
        }
        else {
          ArrayView<const GO> domainElts = domainMap_->getNodeElementList ();
          std::copy (domainElts.begin(), domainElts.end(), LocalColGIDs.begin());
        }
      }
      else {
        // Count the number of locally owned GIDs, both to keep track
        // of the current array index, and as a sanity check.
        size_t numLocalCount = 0;
        if (domainMap_->isContiguous ()) {
          // NOTE (mfh 03 Mar 2013, 02 Sep 2014) In the common case
          // that the domain Map is contiguous, it's more efficient to
          // avoid calling getNodeElementList(), since that
          // permanently constructs and caches the GID list in the
          // contiguous Map.
          GO curColMapGid = domainMap_->getMinGlobalIndex ();
          for (size_t i = 0; i < numDomainElts; ++i, ++curColMapGid) {
            if (GIDisLocal[i]) {
              LocalColGIDs[numLocalCount++] = curColMapGid;
            }
          }
        }
        else {
          ArrayView<const GO> domainElts = domainMap_->getNodeElementList ();
          for (size_t i = 0; i < numDomainElts; ++i) {
            if (GIDisLocal[i]) {
              LocalColGIDs[numLocalCount++] = domainElts[i];
            }
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          numLocalCount != numLocalColGIDs, std::logic_error,
          ": numLocalCount = " << numLocalCount << " != numLocalColGIDs = "
          << numLocalColGIDs << ".  This should never happen.  Please report "
          "this bug to the Tpetra developers.");
      }

      // FIXME (mfh 03 Apr 2013) Now would be a good time to use the
      // information we collected above to construct the Import.  In
      // particular, building an Import requires:
      //
      // 1. numSameIDs (length of initial contiguous sequence of GIDs
      //    on this process that are the same in both Maps; this
      //    equals the number of domain Map elements on this process)
      //
      // 2. permuteToLIDs and permuteFromLIDs (both empty in this
      //    case, since there's no permutation going on; the column
      //    Map starts with the domain Map's GIDs, and immediately
      //    after them come the remote GIDs)
      //
      // 3. remoteGIDs (exactly those GIDs that we found out above
      //    were not in the domain Map) and remoteLIDs (which we could
      //    have gotten above by using the three-argument version of
      //    getRemoteIndexList() that computes local indices as well
      //    as process ranks, instead of the two-argument version that
      //    was used above)
      //
      // 4. remotePIDs (which we have from the getRemoteIndexList()
      //    call above)
      //
      // 5. Sorting remotePIDs, and applying that permutation to
      //    remoteGIDs and remoteLIDs (by calling sort3 above instead
      //    of sort2)
      //
      // 6. Everything after the sort3 call in Import::setupExport():
      //    a. Create the Distributor via createFromRecvs(), which
      //       computes exportGIDs and exportPIDs
      //    b. Compute exportLIDs from exportGIDs (by asking the
      //       source Map, in this case the domain Map, to convert
      //       global to local)
      //
      // Steps 1-5 come for free, since we must do that work anyway in
      // order to compute the column Map.  In particular, Step 3 is
      // even more expensive than Step 6a, since it involves both
      // creating and using a new Distributor object.

    } // if the graph is globally indexed

    const global_size_t gstInv =
      Teuchos::OrdinalTraits<global_size_t>::invalid ();
    // FIXME (mfh 05 Mar 2014) Doesn't the index base of a Map have to
    // be the same as the Map's min GID? If the first column is empty
    // (contains no entries), then the column Map's min GID won't
    // necessarily be the same as the domain Map's index base.
    const GO indexBase = domainMap_->getIndexBase ();
    colMap_ = rcp (new map_type (gstInv, myColumns, indexBase,
                                 domainMap_->getComm (),
                                 domainMap_->getNode ()));

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  mergeAllIndices ()
  {
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed() ); // call only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( ! isSorted() ); // call only after sortIndices()
    if (! isMerged ()) {
      const LocalOrdinal lclNumRows =
        static_cast<LocalOrdinal> (this->getNodeNumRows ());
      for (LocalOrdinal row = 0; row < lclNumRows; ++row) {
        const RowInfo rowInfo = this->getRowInfo (row);
        mergeRowIndices(rowInfo);
      }
      // we just merged every row
      noRedundancies_ = true;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  makeImportExport ()
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(! hasColMap (), std::logic_error, "Tpetra::"
      "CrsGraph::makeImportExport: This method may not be called unless the "
      "graph has a column Map.");
    RCP<ParameterList> params = this->getNonconstParameterList (); // could be null

    // Don't do any checks to see if we need to create the Import, if
    // it exists already.
    //
    // FIXME (mfh 25 Mar 2013) This will become incorrect if we
    // change CrsGraph in the future to allow changing the column
    // Map after fillComplete.  For now, the column Map is fixed
    // after the first fillComplete call.
    if (importer_.is_null ()) {
      // Create the Import instance if necessary.
      if (domainMap_ != colMap_ && (! domainMap_->isSameAs (*colMap_))) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          importer_ = rcp (new import_type (domainMap_, colMap_));
        } else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          importer_ = rcp (new import_type (domainMap_, colMap_, importSublist));
        }
      }
    }

    // Don't do any checks to see if we need to create the Export, if
    // it exists already.
    if (exporter_.is_null ()) {
      // Create the Export instance if necessary.
      if (rangeMap_ != rowMap_ && ! rangeMap_->isSameAs (*rowMap_)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter_ = rcp (new export_type (rowMap_, rangeMap_));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter_ = rcp (new export_type (rowMap_, rangeMap_, exportSublist));
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  std::string
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  description () const
  {
    std::ostringstream oss;
    oss << dist_object_type::description ();
    if (isFillComplete ()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    const char tfecfFuncName[] = "describe()";
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, static_cast<size_t> (11)) + 2;
    Teuchos::OSTab tab (out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print graph indices
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << globalNumDiags_ << std::endl;
        out << "Global max number of row entries = " << globalMaxNumRowEntries_ << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        rowMap_->describe(out,vl);
        if (colMap_ != null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          colMap_->describe(out,vl);
        }
        if (domainMap_ != null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          domainMap_->describe(out,vl);
        }
        if (rangeMap_ != null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          rangeMap_->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << nodeNumEntries_ << std::endl
                << "Node number of diagonals = " << nodeNumDiags_ << std::endl
                << "Node max number of entries = " << nodeMaxNumRowEntries_ << std::endl;
            if (indicesAreAllocated()) {
              out << "Node number of allocated entries = " << nodeNumAllocated_ << std::endl;
            }
            else {
              out << "Indices are not allocated." << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          ! hasRowInfo (), std::runtime_error, ": reduce verbosity level; "
          "graph row information was deleted at fillComplete().");
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << " Entries";
            }
            out << std::endl;
            const LocalOrdinal lclNumRows =
              static_cast<LocalOrdinal> (this->getNodeNumRows ());
            for (LocalOrdinal r=0; r < lclNumRows; ++r) {
              const RowInfo rowinfo = this->getRowInfo (r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << rowinfo.numEntries;
              if (vl == VERB_EXTREME) {
                out << " ";
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowview = getGlobalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << rowview[j] << " ";
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowview = getLocalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << colMap_->getGlobalElement(rowview[j]) << " ";
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  checkSizes (const SrcDistObject& source)
  {
    (void) source; // forestall "unused variable" compiler warnings

    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.
    return true;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  copyAndPermute (const SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "copyAndPermute";
    typedef CrsGraph<LO, GO, node_type> this_type;
    typedef RowGraph<LO, GO, node_type> row_graph_type;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size.");
    // Make sure that the source object has the right type.  We only
    // actually need it to be a RowGraph, with matching first three
    // template parameters.  If it's a CrsGraph, we can use view mode
    // instead of copy mode to get each row's data.
    //
    // FIXME (mfh 07 Jul 2013) It should not be necessary for any of
    // the template parameters but GO to match.  GO has to match
    // because the graph has to send indices as global ordinals, if
    // the source and target graphs do not have the same column Map.
    // If LO doesn't match, the graphs could communicate using global
    // indices.  It could be possible that Node affects the graph's
    // storage format, but packAndPrepare should assume a common
    // communication format in any case.
    const row_graph_type* srcRowGraph = dynamic_cast<const row_graph_type*> (&source);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      srcRowGraph == NULL, std::invalid_argument,
      ": The source object must be a RowGraph with matching first three "
      "template parameters.");

    // If the source object is actually a CrsGraph, we can use view
    // mode instead of copy mode to access the entries in each row,
    // if the graph is not fill complete.
    const this_type* srcCrsGraph = dynamic_cast<const this_type*> (&source);

    const map_type& srcRowMap = * (srcRowGraph->getRowMap ());
    const map_type& tgtRowMap = * (this->getRowMap ());
    const bool src_filled = srcRowGraph->isFillComplete ();
    Array<GO> row_copy;
    LO myid = 0;

    //
    // "Copy" part of "copy and permute."
    //
    if (src_filled || srcCrsGraph == NULL) {
      // If the source graph is fill complete, we can't use view mode,
      // because the data might be stored in a different format not
      // compatible with the expectations of view mode.  Also, if the
      // source graph is not a CrsGraph, we can't use view mode,
      // because RowGraph only provides copy mode access to the data.
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        size_t row_length = srcRowGraph->getNumEntriesInGlobalRow (gid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph->getGlobalRowCopy (gid, row_copy (), check_row_length);
        this->insertGlobalIndices (gid, row_copy ());
      }
    } else {
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (gid, row);
        this->insertGlobalIndices (gid, row);
      }
    }

    //
    // "Permute" part of "copy and permute."
    //
    if (src_filled || srcCrsGraph == NULL) {
      for (LO i = 0; i < permuteToLIDs.size (); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs[i]);
        size_t row_length = srcRowGraph->getNumEntriesInGlobalRow (srcgid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph->getGlobalRowCopy (srcgid, row_copy (), check_row_length);
        this->insertGlobalIndices (mygid, row_copy ());
      }
    } else {
      for (LO i = 0; i < permuteToLIDs.size (); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs[i]);
        ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (srcgid, row);
        this->insertGlobalIndices (mygid, row);
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  packAndPrepare (const SrcDistObject& source,
                  const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                  Teuchos::Array<GlobalOrdinal> &exports,
                  const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor &distor)
  {
    using Teuchos::Array;
    const char tfecfFuncName[] = "packAndPrepare";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": exportLIDs and numPacketsPerLID must have the same size.");
    typedef RowGraph<LocalOrdinal, GlobalOrdinal, node_type> row_graph_type;
    const row_graph_type& srcGraph = dynamic_cast<const row_graph_type&> (source);

    // We don't check whether src_graph has had fillComplete called,
    // because it doesn't matter whether the *source* graph has been
    // fillComplete'd. The target graph can not be fillComplete'd yet.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      this->isFillComplete (), std::runtime_error,
      ": The target graph of an Import or Export must not be fill complete.");

    srcGraph.pack (exportLIDs, exports, numPacketsPerLID, constantNumPackets, distor);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<GlobalOrdinal>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Distributor& distor) const
  {
    using Teuchos::Array;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "pack";
    (void) distor; // forestall "unused argument" compiler warning

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": exportLIDs and numPacketsPerLID must have the same size.");

    const map_type& srcMap = * (this->getMap ());
    constantNumPackets = 0;

    // Set numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.  Compute the total number of
    // packets (that is, entries) owned by this process in all the
    // rows that the caller wants us to send out.
    size_t totalNumPackets = 0;
    Array<GO> row;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      numPacketsPerLID[i] = row_length;
      totalNumPackets += row_length;
    }

    exports.resize (totalNumPackets);

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.
    size_t exportsOffset = 0;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      row.resize (row_length);
      size_t check_row_length = 0;
      this->getGlobalRowCopy (GID, row (), check_row_length);
      typename Array<GO>::const_iterator row_iter = row.begin();
      typename Array<GO>::const_iterator row_end = row.end();
      size_t j = 0;
      for (; row_iter != row_end; ++row_iter, ++j) {
        exports[exportsOffset+j] = *row_iter;
      }
      exportsOffset += row.size();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                    const Teuchos::ArrayView<const GlobalOrdinal> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor& /* distor */,
                    CombineMode /* CM */)
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    // FIXME (mfh 02 Apr 2012) REPLACE combine mode has a perfectly
    // reasonable meaning, whether or not the matrix is fill complete.
    // It's just more work to implement.

    // We are not checking the value of the CombineMode input
    // argument.  For CrsGraph, we only support import/export
    // operations if fillComplete has not yet been called.  Any
    // incoming column-indices are inserted into the target graph. In
    // this context, CombineMode values of ADD vs INSERT are
    // equivalent. What is the meaning of REPLACE for CrsGraph? If a
    // duplicate column-index is inserted, it will be compressed out
    // when fillComplete is called.
    //
    // Note: I think REPLACE means that an existing row is replaced by
    // the imported row, i.e., the existing indices are cleared. CGB,
    // 6/17/2010

    const char tfecfFuncName[] = "unpackAndCombine: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      "importLIDs and numPacketsPerLID must have the same size.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete (), std::runtime_error,
      "Import or Export operations are not allowed on the destination "
      "CrsGraph if it is fill complete.");
    size_t importsOffset = 0;

    typedef typename ArrayView<const LO>::const_iterator iter_type;
    iter_type impLIDiter = importLIDs.begin();
    iter_type impLIDend = importLIDs.end();

    for (size_t i = 0; impLIDiter != impLIDend; ++impLIDiter, ++i) {
      LO row_length = numPacketsPerLID[i];

      const GO* const row_raw = (row_length == 0) ? NULL : &imports[importsOffset];
      ArrayView<const GlobalOrdinal> row (row_raw, row_length);
      insertGlobalIndicesFiltered (this->getMap ()->getGlobalElement (*impLIDiter), row);
      importsOffset += row_length;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // We'll set all the state "transactionally," so that this method
    // satisfies the strong exception guarantee.  This object's state
    // won't be modified until the end of this method.
    RCP<const map_type> rowMap, domainMap, rangeMap, colMap;
    RCP<import_type> importer;
    RCP<export_type> exporter;

    rowMap = newMap;
    RCP<const Comm<int> > newComm =
      (newMap.is_null ()) ? null : newMap->getComm ();

    if (! domainMap_.is_null ()) {
      if (domainMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original domain and row Maps are identical.
        // In that case, we need only replace the original domain Map
        // with the new Map.  This ensures that the new domain and row
        // Maps _stay_ identical.
        domainMap = newMap;
      } else {
        domainMap = domainMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! rangeMap_.is_null ()) {
      if (rangeMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original range and row Maps are identical.  In
        // that case, we need only replace the original range Map with
        // the new Map.  This ensures that the new range and row Maps
        // _stay_ identical.
        rangeMap = newMap;
      } else {
        rangeMap = rangeMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! colMap.is_null ()) {
      colMap = colMap_->replaceCommWithSubset (newComm);
    }

    // (Re)create the Export and / or Import if necessary.
    if (! newComm.is_null ()) {
      RCP<ParameterList> params = this->getNonconstParameterList (); // could be null
      //
      // The operations below are collective on the new communicator.
      //
      // (Re)create the Export object if necessary.  If I haven't
      // called fillComplete yet, I don't have a rangeMap, so I must
      // first check if the _original_ rangeMap is not null.  Ditto
      // for the Import object and the domain Map.
      if (! rangeMap_.is_null () &&
          rangeMap != rowMap &&
          ! rangeMap->isSameAs (*rowMap)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter = rcp (new export_type (rowMap, rangeMap));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter = rcp (new export_type (rowMap, rangeMap, exportSublist));
        }
      }
      // (Re)create the Import object if necessary.
      if (! domainMap_.is_null () &&
          domainMap != colMap &&
          ! domainMap->isSameAs (*colMap)) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          importer = rcp (new import_type (domainMap, colMap));
        } else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          importer = rcp (new import_type (domainMap, colMap, importSublist));
        }
      }
    } // if newComm is not null

    // Defer side effects until the end.  If no destructors throw
    // exceptions (they shouldn't anyway), then this method satisfies
    // the strong exception guarantee.
    exporter_ = exporter;
    importer_ = importer;
    rowMap_ = rowMap;
    // mfh 31 Mar 2013: DistObject's map_ is the row Map of a CrsGraph
    // or CrsMatrix.  CrsGraph keeps a redundant pointer (rowMap_) to
    // the same object.  We might want to get rid of this redundant
    // pointer sometime, but for now, we'll leave it alone and just
    // set map_ to the same object.
    this->map_ = rowMap;
    domainMap_ = domainMap;
    rangeMap_ = rangeMap;
    colMap_ = colMap;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  hasRowInfo () const
  {
    if (indicesAreAllocated () &&
        getProfileType () == StaticProfile &&
        k_rowPtrs_.dimension_0 () == 0) {
      return false;
    } else {
      return true;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalDiagOffsets (const Kokkos::View<size_t*, device_type, Kokkos::MemoryUnmanaged>& offsets) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "getLocalDiagOffsets: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! hasColMap (), std::runtime_error, "The graph must have a column Map.");
    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<LO> (offsets.dimension_0 ()) < lclNumRows,
       std::invalid_argument, "offsets.dimension_0() = " <<
       offsets.dimension_0 () << " < getNodeNumRows() = " << lclNumRows << ".");

    const map_type& rowMap = * (this->getRowMap ());
    const map_type& colMap = * (this->getColMap ());

#ifdef HAVE_TPETRA_DEBUG
    bool allRowMapDiagEntriesInColMap = true;
    bool allDiagEntriesFound = true;
    bool allOffsetsCorrect = true;
    bool noOtherWeirdness = true;
    std::vector<std::pair<LO, size_t> > wrongOffsets;
#endif // HAVE_TPETRA_DEBUG

    // mfh 12 Mar 2016: LocalMap works on (CUDA) device.  It has just
    // the subset of Map functionality that we need below.
    auto lclRowMap = rowMap.getLocalMap ();
    auto lclColMap = colMap.getLocalMap ();

    // FIXME (mfh 16 Dec 2015) It's easy to thread-parallelize this
    // setup, at least on the host.  For CUDA, we have to use LocalMap
    // (that comes from each of the two Maps).

    const bool sorted = this->isSorted ();
    if (isFillComplete ()) {
      auto lclGraph = this->getLocalGraph ();
      // This actually invokes the parallel kernel to do the work.
      Details::GetLocalDiagOffsets<LO, GO, Node> doIt (offsets,
                                                       lclRowMap,
                                                       lclColMap,
                                                       lclGraph.row_map,
                                                       lclGraph.entries,
                                                       sorted);
    }
    else {
      for (LO lclRowInd = 0; lclRowInd < lclNumRows; ++lclRowInd) {
        const GO gblRowInd = lclRowMap.getGlobalElement (lclRowInd);
        const GO gblColInd = gblRowInd;
        const LO lclColInd = lclColMap.getLocalElement (gblColInd);

        if (lclColInd == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
#ifdef HAVE_TPETRA_DEBUG
          allRowMapDiagEntriesInColMap = false;
#endif // HAVE_TPETRA_DEBUG
          offsets[lclRowInd] = Tpetra::Details::OrdinalTraits<size_t>::invalid ();
        }
        else {
          const RowInfo rowInfo = this->getRowInfo (lclRowInd);
          if (static_cast<LO> (rowInfo.localRow) == lclRowInd &&
              rowInfo.numEntries > 0) {

            auto colInds = this->getLocalKokkosRowView (rowInfo);
            const size_t hint = 0; // not needed for this algorithm
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           lclColInd, hint, sorted);
            offsets(lclRowInd) = offset;

#ifdef HAVE_TPETRA_DEBUG
            // Now that we have what we think is an offset, make sure
            // that it really does point to the diagonal entry.  Offsets
            // are _relative_ to each row, not absolute (for the whole
            // (local) graph).
            Teuchos::ArrayView<const LO> lclColInds;
            try {
              this->getLocalRowView (lclRowInd, lclColInds);
            }
            catch (...) {
              noOtherWeirdness = false;
            }
            // Don't continue with error checking if the above failed.
            if (noOtherWeirdness) {
              const size_t numEnt = lclColInds.size ();
              if (offset >= numEnt) {
                // Offsets are relative to each row, so this means that
                // the offset is out of bounds.
                allOffsetsCorrect = false;
                wrongOffsets.push_back (std::make_pair (lclRowInd, offset));
              } else {
                const LO actualLclColInd = lclColInds[offset];
                const GO actualGblColInd = lclColMap.getGlobalElement (actualLclColInd);
                if (actualGblColInd != gblColInd) {
                  allOffsetsCorrect = false;
                  wrongOffsets.push_back (std::make_pair (lclRowInd, offset));
                }
              }
            }
#endif // HAVE_TPETRA_DEBUG
          }
          else {
            offsets(lclRowInd) = Tpetra::Details::OrdinalTraits<size_t>::invalid ();
#ifdef HAVE_TPETRA_DEBUG
            allDiagEntriesFound = false;
#endif // HAVE_TPETRA_DEBUG
          }
        }
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    if (wrongOffsets.size () != 0) {
      std::ostringstream os;
      os << "Proc " << this->getComm ()->getRank () << ": Wrong offsets: [";
      for (size_t k = 0; k < wrongOffsets.size (); ++k) {
        os << "(" << wrongOffsets[k].first << ","
           << wrongOffsets[k].second << ")";
        if (k + 1 < wrongOffsets.size ()) {
          os << ", ";
        }
      }
      os << "]" << std::endl;
      std::cerr << os.str ();
    }
#endif // HAVE_TPETRA_DEBUG

#ifdef HAVE_TPETRA_DEBUG
    using Teuchos::reduceAll;
    using std::endl;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm ();
    const bool localSuccess =
      allRowMapDiagEntriesInColMap && allDiagEntriesFound && allOffsetsCorrect;
    const int numResults = 5;
    int lclResults[5];
    lclResults[0] = allRowMapDiagEntriesInColMap ? 1 : 0;
    lclResults[1] = allDiagEntriesFound ? 1 : 0;
    lclResults[2] = allOffsetsCorrect ? 1 : 0;
    lclResults[3] = noOtherWeirdness ? 1 : 0;
    // min-all-reduce will compute least rank of all the processes
    // that didn't succeed.
    lclResults[4] = ! localSuccess ? comm->getRank () : comm->getSize ();

    int gblResults[5];
    gblResults[0] = 0;
    gblResults[1] = 0;
    gblResults[2] = 0;
    gblResults[3] = 0;
    gblResults[4] = 0;
    reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN,
                         numResults, lclResults, gblResults);

    if (gblResults[0] != 1 || gblResults[1] != 1 || gblResults[2] != 1
        || gblResults[3] != 1) {
      std::ostringstream os; // build error message
      os << "Issue(s) that we noticed (on Process " << gblResults[4] << ", "
        "possibly among others): " << endl;
      if (gblResults[0] == 0) {
        os << "  - The column Map does not contain at least one diagonal entry "
          "of the graph." << endl;
      }
      if (gblResults[1] == 0) {
        os << "  - On one or more processes, some row does not contain a "
          "diagonal entry." << endl;
      }
      if (gblResults[2] == 0) {
        os << "  - On one or more processes, some offsets are incorrect."
           << endl;
      }
      if (gblResults[3] == 0) {
        os << "  - One or more processes had some other error."
           << endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error, os.str());
    }
#endif // HAVE_TPETRA_DEBUG
  }

  namespace { // (anonymous)

    // mfh 21 Jan 2016: This is useful for getLocalDiagOffsets (see
    // below).  The point is to avoid the deep copy between the input
    // Teuchos::ArrayRCP and the internally used Kokkos::View.  We
    // can't use UVM to avoid the deep copy with CUDA, because the
    // ArrayRCP is a host pointer, while the input to the graph's
    // getLocalDiagOffsets method is a device pointer.  Assigning a
    // host pointer to a device pointer is incorrect unless the host
    // pointer points to host pinned memory.  The goal is to get rid
    // of the Teuchos::ArrayRCP overload anyway, so we accept the deep
    // copy for backwards compatibility.
    //
    // We have to use template magic because
    // "staticGraph_->getLocalDiagOffsets(offsetsHosts)" won't compile
    // if device_type::memory_space is not Kokkos::HostSpace (as is
    // the case with CUDA).

    template<class DeviceType,
             const bool memSpaceIsHostSpace =
               std::is_same<typename DeviceType::memory_space,
                            Kokkos::HostSpace>::value>
    struct HelpGetLocalDiagOffsets {};

    template<class DeviceType>
    struct HelpGetLocalDiagOffsets<DeviceType, true> {
      typedef DeviceType device_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> device_offsets_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> host_offsets_type;

      static device_offsets_type
      getDeviceOffsets (const host_offsets_type& hostOffsets)
      {
        // Host and device are the same; no need to allocate a
        // temporary device View.
        return hostOffsets;
      }

      static void
      copyBackIfNeeded (const host_offsets_type& /* hostOffsets */,
                        const device_offsets_type& /* deviceOffsets */)
      { /* copy back not needed; host and device are the same */ }
    };

    template<class DeviceType>
    struct HelpGetLocalDiagOffsets<DeviceType, false> {
      typedef DeviceType device_type;
      // We have to do a deep copy, since host memory space != device
      // memory space.  Thus, the device View is managed (we need to
      // allocate a temporary device View).
      typedef Kokkos::View<size_t*, device_type> device_offsets_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> host_offsets_type;

      static device_offsets_type
      getDeviceOffsets (const host_offsets_type& hostOffsets)
      {
        // Host memory space != device memory space, so we must
        // allocate a temporary device View for the graph.
        return device_offsets_type ("offsets", hostOffsets.dimension_0 ());
      }

      static void
      copyBackIfNeeded (const host_offsets_type& hostOffsets,
                        const device_offsets_type& deviceOffsets)
      {
        Kokkos::deep_copy (hostOffsets, deviceOffsets);
      }
    };
  } // namespace (anonymous)


  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "getLocalDiagOffsets: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::runtime_error,
       "The graph does not yet have a column Map.");
    const LO myNumRows = static_cast<LO> (this->getNodeNumRows ());
    if (static_cast<LO> (offsets.size ()) != myNumRows) {
      // NOTE (mfh 21 Jan 2016) This means that the method does not
      // satisfy the strong exception guarantee (no side effects
      // unless successful).
      offsets.resize (myNumRows);
    }

    // mfh 21 Jan 2016: This method unfortunately takes a
    // Teuchos::ArrayRCP, which is host memory.  The graph wants a
    // device pointer.  We can't access host memory from the device;
    // that's the wrong direction for UVM.  (It's the right direction
    // for inefficient host pinned memory, but we don't want to use
    // that here.)  Thus, if device memory space != host memory space,
    // we allocate and use a temporary device View to get the offsets.
    // If the two spaces are equal, the template magic makes the deep
    // copy go away.
    typedef HelpGetLocalDiagOffsets<device_type> helper_type;
    typedef typename helper_type::host_offsets_type host_offsets_type;
    // Unmanaged host View that views the output array.
    host_offsets_type hostOffsets (offsets.getRawPtr (), myNumRows);
    // Allocate temp device View if host != device, else reuse host array.
    auto deviceOffsets = helper_type::getDeviceOffsets (hostOffsets);
    // NOT recursion; this calls the overload that takes a device View.
    this->getLocalDiagOffsets (deviceOffsets);
    helper_type::copyBackIfNeeded (hostOffsets, deviceOffsets);
  }

} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_CRSGRAPH_GRAPH_INSTANT(LO,GO,NODE) template class CrsGraph< LO , GO , NODE >;

// WARNING: These macros exist only for backwards compatibility.
// We will remove them at some point.
#define TPETRA_CRSGRAPH_SORTROWINDICESANDVALUES_INSTANT(S,LO,GO,NODE)
#define TPETRA_CRSGRAPH_MERGEROWINDICESANDVALUES_INSTANT(S,LO,GO,NODE)
#define TPETRA_CRSGRAPH_ALLOCATEVALUES1D_INSTANT(S,LO,GO,NODE)
#define TPETRA_CRSGRAPH_ALLOCATEVALUES2D_INSTANT(S,LO,GO,NODE)

#define TPETRA_CRSGRAPH_INSTANT(S,LO,GO,NODE)                    \
  TPETRA_CRSGRAPH_SORTROWINDICESANDVALUES_INSTANT(S,LO,GO,NODE)  \
  TPETRA_CRSGRAPH_MERGEROWINDICESANDVALUES_INSTANT(S,LO,GO,NODE) \
  TPETRA_CRSGRAPH_ALLOCATEVALUES1D_INSTANT(S,LO,GO,NODE)         \
  TPETRA_CRSGRAPH_ALLOCATEVALUES2D_INSTANT(S,LO,GO,NODE)

#endif // TPETRA_CRSGRAPH_DEF_HPP
