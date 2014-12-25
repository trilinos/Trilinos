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

#ifndef TPETRA_KOKKOSREFACTOR_CRSGRAPH_DEF_HPP
#define TPETRA_KOKKOSREFACTOR_CRSGRAPH_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_CrsGraph_decl.hpp"
#endif
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

    // Allocate k_numAllocPerRow_ (upper bound on number of entries
    // per row).  We get a host view of the data, as an ArrayRCP.
    // Create a nonowning Kokkos::View of it, copy into
    // k_numAllocPerRow, and sync.  Then assign to k_numAllocPerRow_
    // (which is a const view, so we can't copy into it directly).
    typedef Kokkos::DualView<size_t*, Kokkos::LayoutLeft, device_type>
      dual_view_type;
    typedef typename dual_view_type::host_mirror_space host_type;
    typedef Kokkos::View<const size_t*, Kokkos::LayoutLeft, host_type,
      Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    dual_view_type k_numAllocPerRow ("Tpetra::CrsGraph::numAllocPerRow",
                                     lclNumRows);
    k_numAllocPerRow.template modify<host_type> ();
    Kokkos::deep_copy (k_numAllocPerRow.h_view, numAllocPerRowIn);
    k_numAllocPerRow.template sync<device_type> ();
    k_numAllocPerRow_ = k_numAllocPerRow;

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , k_numAllocPerRow_ (numEntPerRow)
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Kokkos::DualView<const size_t*, device_type>& numEntPerRow,
            const ProfileType pftype,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , nodeNumEntries_ (0)
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (pftype)
    , k_numAllocPerRow_ (numEntPerRow)
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

    // Allocate k_numAllocPerRow_ (upper bound on number of entries
    // per row).  We get a host view of the data, as an ArrayRCP.
    // Create a nonowning Kokkos::View of it, copy into
    // k_numAllocPerRow, and sync.  Then assign to k_numAllocPerRow_
    // (which is a const view, so we can't copy into it directly).
    typedef Kokkos::DualView<size_t*, Kokkos::LayoutLeft, device_type>
      dual_view_type;
    typedef typename dual_view_type::host_mirror_space host_type;
    typedef Kokkos::View<const size_t*, Kokkos::LayoutLeft, host_type,
      Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    dual_view_type k_numAllocPerRow ("Tpetra::CrsGraph::numAllocPerRow",
                                     lclNumRows);
    k_numAllocPerRow.template modify<host_type> ();
    Kokkos::deep_copy (k_numAllocPerRow.h_view, numAllocPerRowIn);
    k_numAllocPerRow.template sync<device_type> ();
    k_numAllocPerRow_ = k_numAllocPerRow;

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const t_RowPtrs& rowPointers,
            const t_LocalOrdinal_1D& columnIndices,
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const LocalStaticCrsGraphType& k_local_graph_,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
    : DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, node_type> (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , k_lclGraph_ (k_local_graph_)
    , globalNumEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalNumDiags_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalMaxNumRowEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , nodeNumEntries_ (0) // FIXME (mfh 17 Mar 2014) should get from k_lclGraph_ right now
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
    using Teuchos::as;
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

    RCP<ParameterList> lclparams;
    if (params.is_null ()) {
      lclparams = parameterList ();
    } else {
      lclparams = sublist (params, "Local Graph");
    }
    // FIXME (mfh 28 Aug 2014) "Local Graph" sublist not used.

    k_lclInds1D_ = k_lclGraph_.entries;
    k_rowPtrs_ = k_lclGraph_.row_map;

    typename LocalStaticCrsGraphType::row_map_type d_ptrs = k_lclGraph_.row_map;
    typename LocalStaticCrsGraphType::entries_type d_inds = k_lclGraph_.entries;

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
          std::max (d_ptrs(rlcid + 1) - d_ptrs(rlcid), nodeMaxNumRowEntries_);
      }
    }

    haveLocalConstants_ = true;
    computeGlobalConstants ();

    fillComplete_ = true;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  ~CrsGraph ()
  {}


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Teuchos::ParameterList>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<const Teuchos::ParameterList> validParams =
      getValidParameters ();
    params->validateParametersAndSetDefaults (*validParams);
    this->setMyParamList (params);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalNumRows () const
  {
    return rowMap_->getGlobalNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalNumCols () const
  {
    const char tfecfFuncName[] = "getGlobalNumCols: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete () || getDomainMap ().is_null (), std::runtime_error,
      "The graph does not have a domain Map.  You may not call this method in "
      "that case.");
    return getDomainMap ()->getGlobalNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeNumRows () const
  {
    return rowMap_.is_null () ? static_cast<size_t> (0) :
      rowMap_->getNodeNumElements ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeNumDiags () const
  {
    return nodeNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalNumDiags () const
  {
    return globalNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getNode () const
  {
    return rowMap_.is_null () ? Teuchos::null : rowMap_->getNode ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRowMap () const
  {
    return rowMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getColMap () const
  {
    return colMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getDomainMap () const
  {
    return domainMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRangeMap () const
  {
    return rangeMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Import<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getImporter () const
  {
    return importer_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Export<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getExporter () const
  {
    return exporter_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  hasColMap () const
  {
    return ! colMap_.is_null ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ProfileType
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getProfileType () const
  {
    return pftype_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalNumEntries () const
  {
    return globalNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeNumEntries () const
  {
    return nodeNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalMaxNumRowEntries () const
  {
    return globalMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeMaxNumRowEntries () const
  {
    return nodeMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isFillComplete () const
  {
    return fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isFillActive () const
  {
    return ! fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isUpperTriangular () const
  {
    return upperTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isLowerTriangular () const
  {
    return lowerTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isLocallyIndexed () const
  {
    return indicesAreLocal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isGloballyIndexed () const
  {
    return indicesAreGlobal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeAllocationSize () const
  {
    return nodeNumAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getComm () const
  {
    return rowMap_.is_null () ? Teuchos::null : rowMap_->getComm ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getIndexBase () const
  {
    return rowMap_->getIndexBase ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  indicesAreAllocated () const
  {
    return indicesAreAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isSorted () const
  {
    return indicesAreSorted_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isMerged () const
  {
    return noRedundancies_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  allocateIndices (const ELocalGlobal lg)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ArrayRCP<size_t>::size_type size_type;
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
      t_RowPtrsNC k_rowPtrs ("Tpetra::CrsGraph::ptr", numRows + 1);

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
        // FIXME hack until we get parallel_scan in kokkos
        //
        // FIXME (mfh 11 Aug 2014) This assumes UVM, since k_rowPtrs_
        // is currently a device View.  Should instead use a DualView.
        typename Kokkos::DualView<const size_t*, device_type>::t_host h_numAllocPerRow =
          k_numAllocPerRow_.h_view; // use a host view for now, since we compute on host
        bool anyInvalidAllocSizes = false;
        for (size_t i = 0; i < numRows; ++i) {
          size_t allocSize = h_numAllocPerRow(i);
          if (allocSize == Teuchos::OrdinalTraits<size_t>::invalid ()) {
            anyInvalidAllocSizes = true;
            allocSize = 0;
          }
          k_rowPtrs(i+1) = k_rowPtrs(i) + allocSize;
        }
        // It's OK to throw std::invalid_argument here, because we
        // haven't incurred any side effects yet.  Throwing that
        // exception (and not, say, std::logic_error) implies that the
        // instance can recover.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          anyInvalidAllocSizes, std::invalid_argument, "The input array of "
          "allocation sizes per row had at least one invalid (== "
          "Teuchos::OrdinalTraits<size_t>::invalid()) entry.");
      }
      else {
        // It's OK to throw std::invalid_argument here, because we
        // haven't incurred any side effects yet.  Throwing that
        // exception (and not, say, std::logic_error) implies that the
        // instance can recover.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          numAllocForAllRows_ == Teuchos::OrdinalTraits<size_t>::invalid (),
          std::invalid_argument, "numAllocForAllRows_ has an invalid value, "
          "namely Teuchos::OrdinalTraits<size_t>::invalid() = " <<
          Teuchos::OrdinalTraits<size_t>::invalid () << ".");
        // FIXME hack until we get parallel_scan in kokkos
        //
        // FIXME (mfh 11 Aug 2014) This assumes UVM, since k_rowPtrs_
        // is currently a device View.  Should instead use a DualView.
        for (size_t i = 0; i < numRows; ++i) {
          k_rowPtrs(i+1) = k_rowPtrs(i) + numAllocForAllRows_;
        }
      }

      // "Commit" the resulting row offsets.
      k_rowPtrs_ = k_rowPtrs;

      // FIXME (mfh 05,11 Aug 2014) This assumes UVM, since k_rowPtrs_
      // is currently a device View.  Should instead use a DualView.
      const size_type numInds = static_cast<size_type> (k_rowPtrs_(numRows));
      if (lg == LocalIndices) {
        k_lclInds1D_ = t_LocalOrdinal_1D ("Tpetra::CrsGraph::ind", numInds);
      }
      else {
        k_gblInds1D_ = t_GlobalOrdinal_1D ("Tpetra::CrsGraph::ind", numInds);
        gblInds1D_ = Kokkos::Compat::persistingView (k_gblInds1D_);
      }
      nodeNumAllocated_ = numInds;
      storageStatus_ = Details::STORAGE_1D_UNPACKED;
    }
    else {
      //
      //  DYNAMIC ALLOCATION PROFILE
      //

      // Use the host view of k_numAllocPerRow_, since we have to
      // allocate 2-D storage on the host.
      typename Kokkos::DualView<const size_t*, device_type>::t_host h_numAllocPerRow =
        k_numAllocPerRow_.h_view;
      const bool useNumAllocPerRow = (k_numAllocPerRow_.dimension_0 () != 0);

      if (lg == LocalIndices) {
        lclInds2D_ = arcp<Array<LocalOrdinal> > (numRows);
        nodeNumAllocated_ = 0;
        for (size_t i = 0; i < numRows; ++i) {
          const size_t howMany = useNumAllocPerRow ?
            h_numAllocPerRow(i) : numAllocForAllRows_;
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
            h_numAllocPerRow(i) : numAllocForAllRows_;
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

    if (numRows > 0) {
      k_numRowEntries_ =
        t_numRowEntries_ ("Tpetra::CrsGraph::numRowEntries", numRows);

      // Fill with zeros on the host.  k_numRowEntries_ is a DualView.
      //
      // TODO (mfh 05 Aug 2014) Write a device kernel for this.
      k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
      size_t* const hostPtr = k_numRowEntries_.h_view.ptr_on_device ();
      std::fill (hostPtr, hostPtr + numRows, static_cast<size_t> (0));
      k_numRowEntries_.template sync<device_type> ();
    }

    // done with these
    numAllocForAllRows_ = 0;
    k_numAllocPerRow_ = Kokkos::DualView<const size_t*, device_type> ();
    indicesAreAllocated_ = true;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  Teuchos::ArrayRCP<Teuchos::Array<T> >
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  allocateValues2D () const
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    const char tfecfFuncName[] = "allocateValues2D: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! indicesAreAllocated (), std::runtime_error,
      "Graph indices must be allocated before values.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getProfileType () != DynamicProfile, std::runtime_error,
      "Graph indices must be allocated in a dynamic profile.");

    ArrayRCP<Array<T> > values2D;
    values2D = arcp<Array<T> > (getNodeNumRows ());
    if (lclInds2D_ != null) {
      const size_t numRows = lclInds2D_.size ();
      for (size_t r = 0; r < numRows; ++r) {
        values2D[r].resize (lclInds2D_[r].size ());
      }
    }
    else if (gblInds2D_ != null) {
      const size_t numRows = gblInds2D_.size ();
      for (size_t r = 0; r < numRows; ++r) {
        values2D[r].resize (gblInds2D_[r].size ());
      }
    }
    return values2D;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayView<const LocalOrdinal>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalView (const RowInfo rowinfo) const
  {
    using Kokkos::subview;
    using Kokkos::View;
    typedef LocalOrdinal LO;
    typedef View<const LO*, device_type, Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return Teuchos::ArrayView<const LO> ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        row_view_type rowView = subview<row_view_type> (k_lclInds1D_, rng);
        return Teuchos::ArrayView<const LO> (rowView.ptr_on_device (), len,
                                             Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        return lclInds2D_[rowinfo.localRow] ();
      }
      else {
        return Teuchos::ArrayView<const LO> (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayView<LocalOrdinal>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalViewNonConst (const RowInfo rowinfo)
  {
    using Kokkos::subview;
    using Kokkos::View;
    typedef LocalOrdinal LO;
    typedef View<LO*, device_type, Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) { // nothing in the row to view
      return Teuchos::ArrayView<LO> ();
    }
    else {
      if (k_lclInds1D_.dimension_0 () != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        row_view_type rowView = subview<row_view_type> (k_lclInds1D_, rng);
        return Teuchos::ArrayView<LO> (rowView.ptr_on_device (), len,
                                       Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty ()) { // 2-D storage
        return lclInds2D_[rowinfo.localRow] ();
      }
      else {
        return Teuchos::ArrayView<LO> (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayView<const GlobalOrdinal>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalView (const RowInfo rowinfo) const
  {
    Teuchos::ArrayView<const GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (! gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayView<GlobalOrdinal>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalViewNonConst (const RowInfo rowinfo)
  {
    Teuchos::ArrayView<GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (!gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RowInfo
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRowInfo (const size_t myRow) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getRowInfo";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! rowMap_->isNodeLocalElement (myRow), std::logic_error,
      ": The given (local) row index myRow = " << myRow
      << " does not belong to the graph's row Map.  "
      "This probably indicates a bug in Tpetra::CrsGraph or Tpetra::CrsMatrix.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::logic_error,
      ": Late catch! Graph does not have row info anymore.  "
      "Error should have been caught earlier.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    ret.localRow = myRow;
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
          ret.numEntries = k_numRowEntries_.h_view(myRow);
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
        ret.numEntries = k_numRowEntries_.h_view(myRow);
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
        ret.allocSize = k_numAllocPerRow_.h_view(myRow);
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  staticAssertions () const
  {
    using Teuchos::OrdinalTraits;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;

    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal):
    //     This is so that we can store local indices in the memory
    //     formerly occupied by global indices.
    Teuchos::CompileTimeAssert< sizeof(GO) < sizeof(LO) > cta_size1;
    (void) cta_size1;

    // Assumption: max(GlobalOrdinal) >= max(LocalOrdinal) and
    //   max(size_t) >= max(LocalOrdinal)
    //     This is so that we can represent any LocalOrdinal as a
    //     size_t, and any LocalOrdinal as a GlobalOrdinal
    Teuchos::CompileTimeAssert< sizeof(GST) < sizeof(size_t) > cta_size2;
    (void) cta_size2;

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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

    // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
    // but for now, for correctness, do the modify-sync cycle here.
    k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
    k_numRowEntries_.h_view(rowinfo.localRow) += numNewInds;
    k_numRowEntries_.template sync<device_type> ();

    nodeNumEntries_ += numNewInds;
    setLocallyModified ();
    return numNewInds;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
    if (gblInds1D_ != null)
      std::copy(indices.begin(), indices.end(),
                gblInds1D_.begin()+rowInfo.offset1D+rowInfo.numEntries);
    else
      std::copy(indices.begin(), indices.end(),
                gblInds2D_[myRow].begin()+rowInfo.numEntries);

    // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
    // but for now, for correctness, do the modify-sync cycle here.
    k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
    k_numRowEntries_.h_view(myRow) += numNewInds;
    k_numRowEntries_.template sync<device_type> ();

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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  insertLocalIndicesImpl (const LocalOrdinal myRow,
                          const Teuchos::ArrayView<const LocalOrdinal>& indices)
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::subview;
    using Kokkos::View;
    typedef LocalOrdinal LO;
    const char* tfecfFuncName ("insertLocallIndicesImpl");

    RowInfo rowInfo = getRowInfo(myRow);
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
      typedef View<const LO*, device_type, MemoryUnmanaged> input_view_type;
      typedef View<LO*, device_type> row_view_type;

      input_view_type inputInds (indices.getRawPtr (), indices.size ());
      const size_t start = rowInfo.offset1D + rowInfo.numEntries; // end of row
      const std::pair<size_t, size_t> rng (start, start + newNumEntries);
      row_view_type myInds = subview<row_view_type> (k_lclInds1D_, rng);
      Kokkos::deep_copy (myInds, inputInds);
    }
    else {
      std::copy (indices.begin (), indices.end (),
                 lclInds2D_[myRow].begin () + rowInfo.numEntries);
    }

    // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
    // but for now, for correctness, do the modify-sync cycle here.
    k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
    k_numRowEntries_.h_view(myRow) += numNewInds;
    k_numRowEntries_.template sync<device_type> ();

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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class Scalar>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  insertIndicesAndValues (const RowInfo& rowInfo,
                          const SLocalGlobalViews& newInds,
                          const Teuchos::ArrayView<Scalar>& oldRowVals,
                          const Teuchos::ArrayView<const Scalar>& newRowVals,
                          const ELocalGlobal lg,
                          const ELocalGlobal I)
  {
#ifdef HAVE_TPETRA_DEBUG
    size_t numNewInds = 0;
    try {
      numNewInds = insertIndices (rowInfo, newInds, lg, I);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsGraph::insertIndicesAndValues: "
        "insertIndices threw an exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      numNewInds > oldRowVals.size (), std::runtime_error,
      "Tpetra::CrsGraph::insertIndicesAndValues: numNewInds (" << numNewInds
      << ") > oldRowVals.size() (" << oldRowVals.size () << ".");
#else
    const size_t numNewInds = insertIndices (rowInfo, newInds, lg, I);
#endif // HAVE_TPETRA_DEBUG

    typedef typename Teuchos::ArrayView<Scalar>::size_type size_type;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      rowInfo.numEntries + numNewInds > oldRowVals.size (), std::runtime_error,
      "Tpetra::CrsGraph::insertIndicesAndValues: rowInfo.numEntries (" <<
      rowInfo.numEntries << ") + numNewInds (" << numNewInds <<
      ") > oldRowVals.size() (" << oldRowVals.size () << ").");
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_type> (numNewInds) > newRowVals.size (),
      std::runtime_error, "Tpetra::CrsGraph::insertIndicesAndValues: "
      "numNewInds (" << numNewInds << ") > newRowVals.size() ("
      << newRowVals.size () << ").");
#endif // HAVE_TPETRA_DEBUG

    size_type oldInd = static_cast<size_type> (rowInfo.numEntries);
#ifdef HAVE_TPETRA_DEBUG
    try {
#endif // HAVE_TPETRA_DEBUG
      for (size_type newInd = 0; newInd < static_cast<size_type> (numNewInds);
           ++newInd, ++oldInd) {
        oldRowVals[oldInd] = newRowVals[newInd];
      }
#ifdef HAVE_TPETRA_DEBUG
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsGraph::insertIndicesAndValues: "
        "for loop for copying values threw an exception: " << e.what ());
    }
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  sortRowIndices (const RowInfo rowinfo)
  {
    if (rowinfo.numEntries > 0) {
      Teuchos::ArrayView<LocalOrdinal> inds_view =
        this->getLocalViewNonConst (rowinfo);
      std::sort (inds_view.begin (), inds_view.begin () + rowinfo.numEntries);
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class Scalar>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  sortRowIndicesAndValues (const RowInfo rowinfo,
                           const Teuchos::ArrayView<Scalar>& values)
  {
    if (rowinfo.numEntries > 0) {
      Teuchos::ArrayView<LocalOrdinal> inds_view =
        this->getLocalViewNonConst (rowinfo);
      sort2 (inds_view.begin (), inds_view.begin () + rowinfo.numEntries,
             values.begin ());
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

    // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
    // but for now, for correctness, do the modify-sync cycle here.
    k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
    k_numRowEntries_.h_view(rowinfo.localRow) = mergedEntries;
    k_numRowEntries_.template sync<device_type> ();

    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template<class Scalar>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  mergeRowIndicesAndValues (RowInfo rowinfo,
                            const Teuchos::ArrayView<Scalar>& rowValues)
  {
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "mergeRowIndicesAndValues: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized(), std::logic_error, "It is invalid to call this "
      "method if the graph's storage has already been optimized.  Please "
      "report this bug to the Tpetra developers.");

    typedef typename ArrayView<Scalar>::iterator Iter;
    Iter rowValueIter = rowValues.begin ();
    ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst (rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;

    // beg,end define a half-exclusive interval over which to iterate.
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = beg;
    if (beg != end) {
      typename ArrayView<LocalOrdinal>::iterator cur = beg + 1;
      Iter vcur = rowValueIter + 1;
      Iter vend = rowValueIter;
      cur = beg+1;
      while (cur != end) {
        if (*cur != *newend) {
          // new entry; save it
          ++newend;
          ++vend;
          (*newend) = (*cur);
          (*vend) = (*vcur);
        }
        else {
          // old entry; merge it
          //(*vend) = f (*vend, *vcur);
          (*vend) += *vcur;
        }
        ++cur;
        ++vcur;
      }
      ++newend; // one past the last entry, per typical [beg,end) semantics
    }
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the
    // assignment below will destroy the packed structure
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized() && mergedEntries != rowinfo.numEntries,
      std::logic_error,
      ": Merge was incorrect; it eliminated entries from the graph.  "
      << std::endl << "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
    // but for now, for correctness, do the modify-sync cycle here.
    k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
    k_numRowEntries_.h_view(rowinfo.localRow) = mergedEntries;
    k_numRowEntries_.template sync<device_type> ();

    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  findLocalIndex (RowInfo rowinfo, LocalOrdinal ind, size_t hint) const
  {
    using Teuchos::ArrayView;
    ArrayView<const LocalOrdinal> colInds = this->getLocalView (rowinfo);
    return this->findLocalIndex (rowinfo, ind, colInds, hint);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  findLocalIndex (RowInfo rowinfo,
                  LocalOrdinal ind,
                  Teuchos::ArrayView<const LocalOrdinal> colInds,
                  size_t hint) const
  {
    typedef typename Teuchos::ArrayView<const LocalOrdinal>::iterator IT;

    // If the hint was correct, then the hint is the offset to return.
    if (hint < rowinfo.numEntries && colInds[hint] == ind) {
      return hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.  How we do the
    // search depends on whether the graph's column indices are
    // sorted.
    IT beg = colInds.begin ();
    IT end = beg + rowinfo.numEntries;
    IT ptr = beg + rowinfo.numEntries; // "null"
    bool found = true;

    if (isSorted ()) {
      std::pair<IT,IT> p = std::equal_range (beg, end, ind); // binary search
      if (p.first == p.second) {
        found = false;
      } else {
        ptr = p.first;
      }
    }
    else {
      ptr = std::find (beg, end, ind); // direct search
      if (ptr == end) {
        found = false;
      }
    }

    if (found) {
      return static_cast<size_t> (ptr - beg);
    }
    else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  findGlobalIndex (RowInfo rowinfo, GlobalOrdinal ind, size_t hint) const
  {
    using Teuchos::ArrayView;
    typedef typename ArrayView<const GlobalOrdinal>::iterator IT;

    // Don't let an invalid global column index through.
    if (ind == Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()) {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }

    ArrayView<const GlobalOrdinal> indices = getGlobalView (rowinfo);

    // We don't actually require that the hint be a valid index.
    // If it is not in range, we just ignore it.
    if (hint < rowinfo.numEntries && indices[hint] == ind) {
      return hint;
    }

    IT beg = indices.begin ();
    IT end = indices.begin () + rowinfo.numEntries; // not indices.end()
    if (isSorted ()) { // use binary search
      const std::pair<IT,IT> p = std::equal_range (beg, end, ind);
      if (p.first == p.second) { // range of matching entries is empty
        return Teuchos::OrdinalTraits<size_t>::invalid ();
      } else {
        return p.first - beg;
      }
    }
    else { // not sorted; must use linear search
      const IT loc = std::find (beg, end, ind);
      if (loc == end) {
        return Teuchos::OrdinalTraits<size_t>::invalid ();
      } else {
        return loc - beg;
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  clearGlobalConstants ()
  {
    globalNumEntries_       = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    globalNumDiags_         = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    globalMaxNumRowEntries_ = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    haveGlobalConstants_    = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
      (k_numRowEntries_.dimension_0 () == 0 || (lclInds2D_.is_null () && gblInds2D_.is_null ())),
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    const LocalOrdinal lrow = rowMap_->getLocalElement (globalRow);
    if (hasRowInfo () && lrow != OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowinfo = getRowInfo (lrow);
      return rowinfo.numEntries;
    } else {
      return OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNumEntriesInLocalRow (LocalOrdinal localRow) const
  {
    if (hasRowInfo () && rowMap_->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = getRowInfo (localRow);
      return rowinfo.numEntries;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNumAllocatedEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement (globalRow);
    if (hasRowInfo () && lrow != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowinfo = getRowInfo (lrow);
      return rowinfo.allocSize;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNumAllocatedEntriesInLocalRow (LocalOrdinal localRow) const
  {
    if (hasRowInfo () && rowMap_->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = getRowInfo (localRow);
      return rowinfo.allocSize;
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const size_t>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeRowPtrs () const
  {
    return Kokkos::Compat::persistingView (k_rowPtrs_);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const LocalOrdinal>
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodePackedIndices () const
  {
    return Kokkos::Compat::persistingView (k_lclInds1D_);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalRowCopy (LocalOrdinal localRow,
                   const Teuchos::ArrayView<LocalOrdinal>&indices,
                   size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed () && ! hasColMap (), std::runtime_error,
      "Tpetra::CrsGraph::getLocalRowCopy: The graph is globally indexed and "
      "does not have a column Map yet.  That means we don't have local indices "
      "for columns yet, so it doesn't make sense to call this method.  If the "
      "graph doesn't have a column Map yet, you should call fillComplete on "
      "it first.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! hasRowInfo(), std::runtime_error,
      "Tpetra::CrsGraph::getLocalRowCopy: graph row information was deleted "
      "at fillComplete().");

    if (! getRowMap ()->isNodeLocalElement (localRow)) {
      numEntries = 0;
      return;
    }

    const RowInfo rowinfo = getRowInfo (localRow);
    const size_t theNumEntries = rowinfo.numEntries;

    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (indices.size ()) < theNumEntries,
      std::runtime_error,
      "Tpetra::CrsGraph::getLocalRowCopy: The given row " << localRow << " has "
      << theNumEntries << " entries, but indices.size() = " << indices.size ()
      << ", which does not suffice to store the row's indices.");

    numEntries = theNumEntries;

    if (isLocallyIndexed ()) {
      ArrayView<const LO> lview = getLocalView (rowinfo);
      std::copy (lview.begin (), lview.begin () + numEntries, indices.begin ());
    }
    else if (isGloballyIndexed ()) {
      ArrayView<const GO> gview = getGlobalView (rowinfo);
      const map_type& colMap = * (this->getColMap ());
      for (size_t j = 0; j < numEntries; ++j) {
        indices[j] = colMap.getLocalElement (gview[j]);
      }
    }
    else {
      // If the graph on the calling process is neither locally nor
      // globally indexed, that means it owns no column indices.
      //
      // FIXME (mfh 21 Oct 2013) It's not entirely clear to me whether
      // we can reach this branch, given the checks above.  However,
      // if that is the case, it should still be correct to call this
      // function if the calling process owns no column indices.
      numEntries = 0;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalRowCopy (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<GlobalOrdinal>& indices,
                    size_t& NumIndices) const
  {
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "getGlobalRowCopy: ";
    // we either currently store global indices, or we have a column
    // map with which to transcribe our local indices for the user
    const LocalOrdinal lrow = rowMap_->getLocalElement (globalRow);

    // FIXME (mfh 22 Aug 2014) Instead of throwing an exception,
    // should just set NumIndices=0 and return.  In that case, the
    // calling process owns no entries in that row, so the right thing
    // to do is to return an empty copy.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid (),
      std::runtime_error,
      "GlobalRow (== " << globalRow << ") does not belong to this process.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      "Graph row information was deleted at fillComplete().");
    const RowInfo rowinfo = this->getRowInfo (static_cast<size_t> (lrow));
    NumIndices = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < NumIndices, std::runtime_error,
      "Specified storage (size==" << indices.size () << ") does not suffice "
      "to hold all entries for this row (NumIndices == " << NumIndices << ").");
    if (isLocallyIndexed ()) {
      ArrayView<const LocalOrdinal> lview = this->getLocalView (rowinfo);
      for (size_t j = 0; j < NumIndices; ++j) {
        indices[j] = colMap_->getGlobalElement (lview[j]);
      }
    }
    else if (isGloballyIndexed ()) {
      ArrayView<const GlobalOrdinal> gview = this->getGlobalView (rowinfo);
      std::copy (gview.begin (), gview.begin () + NumIndices, indices.begin ());
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalRowView (LocalOrdinal localRow,
                   Teuchos::ArrayView<const LocalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getLocalRowView: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "The graph's indices are "
      "currently stored as global indices, so we cannot return a view with "
      "local column indices, whether or not the graph has a column Map.  If "
      "the graph _does_ have a column Map, use getLocalRowCopy() instead.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error, "Graph row information was "
      "deleted at fillComplete().");
    indices = Teuchos::null;
    if (rowMap_->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = this->getRowInfo (localRow);
      if (rowinfo.numEntries > 0) {
        indices = this->getLocalView (rowinfo);
        // getLocalView returns a view of the _entire_ row, including
        // any extra space at the end (which 1-D unpacked storage
        // might have, for example).  That's why we have to take a
        // subview of the returned view.
        indices = indices (0, rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) != this->getNumEntriesInLocalRow (localRow),
      std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalRowView (GlobalOrdinal globalRow,
                    Teuchos::ArrayView<const GlobalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getGlobalRowView: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error, "The graph's indices are "
      "currently stored as local indices, so we cannot return a view with "
      "global column indices.  Use getGlobalRowCopy() instead.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      "Graph row information was deleted at fillComplete().");

    // isNodeGlobalElement() requires a global to local lookup anyway,
    // and getLocalElement() returns invalid() if the element wasn't found.
    const LocalOrdinal localRow = rowMap_->getLocalElement (globalRow);
    indices = Teuchos::null;
    if (localRow != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowInfo = getRowInfo (static_cast<size_t> (localRow));
      if (rowInfo.numEntries > 0) {
        indices = (this->getGlobalView (rowInfo)) (0, rowInfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) != this->getNumEntriesInGlobalRow (globalRow),
      std::logic_error,
      "Violated stated postconditions: indices.size() = " << indices.size ()
      << " != getNumEntriesInGlobalRow(globalRow=" << globalRow
      << ") = " << this->getNumEntriesInGlobalRow (globalRow)
      << ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
      const size_t oldNumEntries = k_numRowEntries_.h_view (lrow);
      nodeNumEntries_ -= oldNumEntries;

      // FIXME (mfh 07 Aug 2014) We should just sync at fillComplete,
      // but for now, for correctness, do the modify-sync cycle here.
      k_numRowEntries_.template modify<typename t_numRowEntries_::host_mirror_space> ();
      k_numRowEntries_.h_view(lrow) = 0;
      k_numRowEntries_.template sync<device_type> ();
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNumEntriesInLocalRow (lrow) != 0 ||
      ! indicesAreAllocated () ||
      ! isLocallyIndexed (), std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  setAllIndices (const t_RowPtrs& rowPointers,
                 const t_LocalOrdinal_1D& columnIndices)
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

    // These normally get cleared out at the end of allocateIndices.
    // It makes sense to clear them out here, because at the end of
    // this method, the graph is allocated on the calling process.
    numAllocForAllRows_ = 0;
    k_numAllocPerRow_ = Kokkos::DualView<const size_t*, device_type> ();

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  setAllIndices (const Teuchos::ArrayRCP<size_t>& rowPointers,
                 const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices)
  {
    Kokkos::View<size_t*, device_type> k_ptr =
      Kokkos::Compat::getKokkosViewDeepCopy<DeviceType> (rowPointers ());
    Kokkos::View<LocalOrdinal*, device_type> k_ind =
      Kokkos::Compat::getKokkosViewDeepCopy<DeviceType> (columnIndices ());

    setAllIndices (k_ptr, k_ind);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
        numEntriesPerRow = Kokkos::Compat::persistingView (k_numAllocPerRow_.h_view);
        allRowsSame = false; // conservatively; we don't check the array
      }
      else {
        numEntriesForAll = numAllocForAllRows_;
        allRowsSame = true;
      }
    }
    else if (k_numRowEntries_.dimension_0 () != 0) {
      numEntriesPerRow = Kokkos::Compat::persistingView (k_numRowEntries_.h_view);
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
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  globalAssemble ()
  {
    using Teuchos::Array;
    using Teuchos::as;
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
      as<typename Array<int>::size_type> (numSends) != sendIDs.size (),
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                            const Teuchos::RCP<const map_type>& rangeMap,
                            const Teuchos::RCP<const import_type>& importer,
                            const Teuchos::RCP<const export_type>& exporter,
                            const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "expertStaticFillComplete: ";

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
    k_numAllocPerRow_    = Kokkos::DualView<const size_t*, device_type> ();
    indicesAreAllocated_ = true;

    // Constants from makeIndicesLocal
    //
    // The graph has a column Map, so its indices had better be local.
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;

    // set domain/range map: may clear the import/export objects
    setDomainRangeMaps (domainMap, rangeMap);

    // Presume the user sorted and merged the arrays first
    indicesAreSorted_ = true;
    noRedundancies_ = true;

    // makeImportExport won't create a new importer/exporter if I set one here first.
    importer_ = Teuchos::null;
    exporter_ = Teuchos::null;
    if (importer != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! importer->getSourceMap ()->isSameAs (*getDomainMap ()) ||
        ! importer->getTargetMap ()->isSameAs (*getColMap ()),
        std::invalid_argument,": importer does not match matrix maps.");
      importer_ = importer;

    }
    if (exporter != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! exporter->getSourceMap ()->isSameAs (*getRowMap ()) ||
        ! exporter->getTargetMap ()->isSameAs (*getRangeMap ()),
        std::invalid_argument,": exporter does not match matrix maps.");
      exporter_ = exporter;
    }
    makeImportExport ();

    // Compute the constants
    computeGlobalConstants ();

    // Since we have a StaticProfile, fillLocalGraph will do the right thing...
    fillLocalGraph (params);
    fillComplete_ = true;

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Kokkos::create_mirror_view;
    using Teuchos::ArrayRCP;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef ArrayRCP<size_t>::size_type size_type;
    typedef t_numRowEntries_ row_entries_type;
    typedef t_RowPtrsNC row_offsets_type;
    typedef t_LocalOrdinal_1D lclinds_1d_type;

    const size_t lclNumRows = this->getNodeNumRows ();

    // This method's goal is to fill in the two arrays (compressed
    // sparse row format) that define the sparse graph's structure.
    //
    // Use t_RowPtrs and not LocalStaticCrsGraphType::row_map_type for
    // k_ptrs, because the latter is const and we need to modify
    // k_ptrs here.
    row_offsets_type k_ptrs;
    t_RowPtrs k_ptrs_const;
    lclinds_1d_type k_inds;

    // The number of entries in each locally owned row.  This is a
    // DualView.  2-D storage lives on host and is currently not
    // thread-safe for parallel kernels even on host, so we have to
    // work sequentially with host storage in that case.
    row_entries_type k_numRowEnt = k_numRowEntries_;
    typename row_entries_type::t_host h_numRowEnt = k_numRowEnt.h_view;

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
      // (k_inds) and then copy from 2-D storage (lclInds2D_) into 1-D
      // storage (k_inds).
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_numRowEnt.dimension_0 ()) != lclNumRows,
        std::logic_error, "Tpetra::CrsGraph::fillLocalGraph (called from "
        "fillComplete or expertStaticFillComplete): For the DynamicProfile "
        "branch, k_numRowEnt has the wrong length.  "
        "k_numRowEnt.dimension_0() = " << k_numRowEnt.dimension_0 ()
        << " != getNodeNumRows() = " << lclNumRows << "");

      // Pack the row offsets into k_ptrs, by doing a sum-scan of
      // the array of valid entry counts per row (h_numRowEnt).
      //
      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      // This will be a host view of packed row offsets.
      typename row_offsets_type::HostMirror h_ptrs;
      {
        // Allocate the packed row offsets array.
        k_ptrs = row_offsets_type ("Tpetra::CrsGraph::ptr", lclNumRows+1);
        k_ptrs_const = k_ptrs;
        //
        // FIXME hack until we get parallel_scan in kokkos
        //
        h_ptrs = create_mirror_view (k_ptrs);
        h_ptrs(0) = 0;
        for (size_type i = 0; i < static_cast<size_type> (lclNumRows); ++i) {
          const size_t numEnt = h_numRowEnt(i);
          lclTotalNumEntries += numEnt;
          h_ptrs(i+1) = h_ptrs(i) + numEnt;
        }
        Kokkos::deep_copy (k_ptrs, h_ptrs);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: In DynamicProfile "
        "branch, after packing k_ptrs, k_ptrs.dimension_0() = "
        << k_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (h_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: In DynamicProfile "
        "branch, after packing h_ptrs, h_ptrs.dimension_0() = "
        << h_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      // FIXME (mfh 08 Aug 2014) This assumes UVM.
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_ptrs(lclNumRows) != lclTotalNumEntries, std::logic_error,
        "Tpetra::CrsGraph::fillLocalGraph: In DynamicProfile branch, after "
        "packing k_ptrs, k_ptrs(lclNumRows = " << lclNumRows << ") = " <<
        k_ptrs(lclNumRows) << " != total number of entries on the calling "
        "process = " << lclTotalNumEntries << ".");

      // Allocate the array of packed column indices.
      k_inds = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);
      // We need a host view of the above, since 2-D storage lives on host.
      typename lclinds_1d_type::HostMirror h_inds = create_mirror_view (k_inds);
      // Pack the column indices.
      for (size_t row = 0; row < lclNumRows; ++row) {
        const size_t numEnt = h_numRowEnt(row);
        std::copy (lclInds2D_[row].begin (),
                   lclInds2D_[row].begin () + numEnt,
                   h_inds.ptr_on_device () + h_ptrs(row));
      }
      Kokkos::deep_copy (k_inds, h_inds);

      // Sanity check of packed row offsets.
      if (k_ptrs.dimension_0 () != 0) {
        const size_t numOffsets = static_cast<size_t> (k_ptrs.dimension_0 ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs(numOffsets-1)) != k_inds.dimension_0 (),
          std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: "
          "In DynamicProfile branch, after packing, k_ptrs(" << (numOffsets-1)
          << ") = " << k_ptrs(numOffsets-1) << " != k_inds.dimension_0() = "
          << k_inds.dimension_0 () << ".");
      }
    }
    else if (getProfileType () == StaticProfile) {
      // StaticProfile means that the graph's column indices are
      // currently stored in a 1-D format, with row offsets in
      // k_rowPtrs_ and local column indices in k_lclInds1D_.

      // StaticProfile also means that the graph's array of row
      // offsets must already be allocated.
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_rowPtrs_.dimension_0 () == 0, std::logic_error,
        "k_rowPtrs_ has size zero, but shouldn't");
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_rowPtrs_.dimension_0 () != lclNumRows + 1, std::logic_error,
        "Tpetra::CrsGraph::fillLocalGraph: k_rowPtrs_ has size "
        << k_rowPtrs_.dimension_0 () << " != (lclNumRows + 1) = "
        << (lclNumRows + 1) << ".")
      {
        const size_t numOffsets = k_rowPtrs_.dimension_0 ();
        // FIXME (mfh 08 Aug 2014) This relies on UVM.
        TEUCHOS_TEST_FOR_EXCEPTION(
          numOffsets != 0 &&
          k_lclInds1D_.dimension_0 () != k_rowPtrs_(numOffsets-1),
          std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: "
          "numOffsets = " << numOffsets << " != 0 and "
          "k_lclInds1D_.dimension_0() = " << k_lclInds1D_.dimension_0 ()
          << " != k_rowPtrs_(" << numOffsets << ") = "
          << k_rowPtrs_(numOffsets-1) << ".");
      }

      if (nodeNumEntries_ != nodeNumAllocated_) {
        // The graph's current 1-D storage is "unpacked."  This means
        // the row offsets may differ from what the final row offsets
        // should be.  This could happen, for example, if the user
        // specified StaticProfile in the constructor and set an upper
        // bound on the number of entries in each row, but didn't fill
        // all those entries.
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_numRowEnt.dimension_0 ()) != lclNumRows,
          std::logic_error, "Tpetra::CrsGraph::fillLocalGraph (called from "
          "fillComplete or expertStaticFillComplete): In StaticProfile "
          "unpacked branch, k_numRowEnt has the wrong length.  "
          "k_numRowEnt.dimension_0() = " << k_numRowEnt.dimension_0 ()
          << " != getNodeNumRows() = " << lclNumRows << "");

        if (k_rowPtrs_.dimension_0 () != 0) {
          const size_t numOffsets =
            static_cast<size_t> (k_rowPtrs_.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            k_rowPtrs_(numOffsets-1) != static_cast<size_t> (k_lclInds1D_.dimension_0 ()),
            std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: "
            "In StaticProfile branch, before allocating or packing, "
            "k_rowPtrs_(" << (numOffsets-1) << ") = "
            << k_rowPtrs_(numOffsets-1) << " != k_lclInds1D_.dimension_0() = "
            << k_lclInds1D_.dimension_0 () << ".");
        }

        // Pack the row offsets into k_ptrs, by doing a sum-scan of
        // the array of valid entry counts per row (h_numRowEnt).

        // Total number of entries in the matrix on the calling
        // process.  We will compute this in the loop below.  It's
        // cheap to compute and useful as a sanity check.
        size_t lclTotalNumEntries = 0;
        {
          // Allocate the packed row offsets array.
          k_ptrs = row_offsets_type ("Tpetra::CrsGraph::ptr", lclNumRows + 1);
          k_ptrs_const = k_ptrs;
          //
          // FIXME hack until we get parallel_scan in kokkos
          //
          // Unlike in the 2-D storage case above, we don't need the
          // host view of the packed row offsets array after packing
          // the row offsets.
          typename row_offsets_type::HostMirror h_k_ptrs =
            create_mirror_view (k_ptrs);
          h_k_ptrs(0) = 0;
          for (size_t i = 0; i < lclNumRows; ++i) {
            const size_t numEnt = h_numRowEnt(i);
            lclTotalNumEntries += numEnt;
            h_k_ptrs(i+1) = h_k_ptrs(i) + numEnt;
          }
          Kokkos::deep_copy (k_ptrs, h_k_ptrs);
        }

        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs.dimension_0 ()) != lclNumRows + 1,
          std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: In "
          "StaticProfile unpacked branch, after allocating k_ptrs, "
          "k_ptrs.dimension_0() = " << k_ptrs.dimension_0 () << " != "
          "lclNumRows+1 = " << (lclNumRows+1) << ".");
        // FIXME (mfh 08 Aug 2014) This assumes UVM.
        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs(lclNumRows) != lclTotalNumEntries, std::logic_error,
          "Tpetra::CrsGraph::fillLocalGraph: In StaticProfile unpacked "
          "branch, after filling k_ptrs, k_ptrs(lclNumRows=" << lclNumRows
          << ") = " << k_ptrs(lclNumRows) << " != total number of entries on "
          "the calling process = " << lclTotalNumEntries << ".");

        // Allocate the array of packed column indices.
        k_inds = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);

        // k_rowPtrs_ and k_lclInds1D_ are currently unpacked.  Pack
        // them, using the packed row offsets array k_ptrs that we
        // created above.
        //
        // FIXME (mfh 08 Aug 2014) If "Optimize Storage" is false (in
        // CrsMatrix?), we need to keep around the unpacked row
        // offsets and column indices.

        // Pack the column indices from unpacked k_lclInds1D_ into
        // packed k_inds.  We will replace k_lclInds1D_ below.
        typedef pack_functor<t_LocalOrdinal_1D, t_RowPtrs> inds_packer_type;
        inds_packer_type f (k_inds, k_lclInds1D_, k_ptrs, k_rowPtrs_);
        Kokkos::parallel_for (lclNumRows, f);

        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs.dimension_0 () == 0, std::logic_error, "Tpetra::CrsGraph::"
          "fillLocalGraph: In StaticProfile \"Optimize Storage\" = true branch,"
          " after packing, k_ptrs.dimension_0() = 0.  This probably means that "
          "k_rowPtrs_ was never allocated.");
        if (k_ptrs.dimension_0 () != 0) {
          const size_t numOffsets = static_cast<size_t> (k_ptrs.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs(numOffsets - 1)) != k_inds.dimension_0 (),
            std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: "
            "In StaticProfile \"Optimize Storage\"=true branch, after packing, "
            "k_ptrs(" << (numOffsets-1) << ") = " << k_ptrs(numOffsets-1) <<
            " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
        }
      }
      else { // We don't have to pack, so just set the pointers.
        k_ptrs_const = k_rowPtrs_;
        k_inds = k_lclInds1D_;

        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs_const.dimension_0 () == 0, std::logic_error, "Tpetra::CrsGraph::"
          "fillLocalGraph: In StaticProfile \"Optimize Storage\" = "
          "false branch, k_ptrs_const.dimension_0() = 0.  This probably means that "
          "k_rowPtrs_ was never allocated.");
        if (k_ptrs_const.dimension_0 () != 0) {
          const size_t numOffsets = static_cast<size_t> (k_ptrs_const.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_inds.dimension_0 (),
            std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: "
            "In StaticProfile \"Optimize Storage\" = false branch, "
            "k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets - 1)
            << " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
        }
      }
    }

    // Extra sanity checks.
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (k_ptrs_const.dimension_0 ()) != lclNumRows + 1,
      std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: After packing, "
      "k_ptrs_const.dimension_0() = " << k_ptrs_const.dimension_0 ()
      << " != lclNumRows+1 = " << (lclNumRows+1) << ".");
    if (k_ptrs_const.dimension_0 () != 0) {
      const size_t numOffsets = static_cast<size_t> (k_ptrs_const.dimension_0 ());
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_inds.dimension_0 (),
        std::logic_error, "Tpetra::CrsGraph::fillLocalGraph: After packing, "
        "k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets-1)
        << " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
    }

    if (requestOptimizedStorage) {
      // With optimized storage, we don't need to store the 2-D column
      // indices array-of-arrays, or the array of row entry counts.

      // Free graph data structures that are only needed for 2-D or
      // unpacked 1-D storage.
      lclInds2D_ = null;
      k_numRowEntries_ = row_entries_type ();

      // Keep the new 1-D packed allocations.
      k_rowPtrs_   = k_ptrs_const;
      k_lclInds1D_ = k_inds;

      // Storage is packed now, so the number of allocated entries is
      // the same as the actual number of entries.
      nodeNumAllocated_ = nodeNumEntries_;
      // The graph is definitely StaticProfile now, whether or not it
      // was before.
      pftype_ = StaticProfile;
    }

    // FIXME (mfh 28 Aug 2014) "Local Graph" sublist no longer used.

    // Build the local graph.
    k_lclGraph_ = LocalStaticCrsGraphType (k_inds, k_ptrs_const);

    // TODO (mfh 13 Mar 2014) getNodeNumDiags(), isUpperTriangular(),
    // and isLowerTriangular() depend on computeGlobalConstants(), in
    // particular the part where it looks at the local matrix.  You
    // have to use global indices to determine which entries are
    // diagonal, or above or below the diagonal.  However, lower or
    // upper triangularness is a local property.
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  reindexColumns (const Teuchos::RCP<const map_type>& newColMap,
                  const Teuchos::RCP<const import_type>& newImport,
                  const bool sortIndicesInEachRow)
  {
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
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

    const size_t lclNumRows = getNodeNumRows ();

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
    t_LocalOrdinal_1D newLclInds1D;
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
            newLclInds1D = t_LocalOrdinal_1D ("Tpetra::CrsGraph::ind",
                                              nodeNumAllocated_);
            // Attempt to convert the new indices locally.
            for (size_t lclRow = 0; lclRow < lclNumRows; ++lclRow) {
              const RowInfo rowInfo = getRowInfo (lclRow);
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
            for (size_t lclRow = 0; lclRow < lclNumRows; ++lclRow) {
              const RowInfo rowInfo = getRowInfo (lclRow);
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
        for (size_t lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const RowInfo rowInfo = getRowInfo (lclRow);
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::LocalStaticCrsGraphType
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalGraph_Kokkos () const
  {
    return k_lclGraph_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  computeGlobalConstants ()
  {
    using Teuchos::as;
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
        const size_t numLocalRows = getNodeNumRows ();
        for (size_t localRow = 0; localRow < numLocalRows; ++localRow) {
          const GO globalRow = rowMap.getGlobalElement (localRow);
          // Find the local (column) index for the diagonal entry.
          // This process might not necessarily own _any_ entries in
          // the current row.  If it doesn't, skip this row.  It won't
          // affect any of the attributes (nodeNumDiagons_,
          // upperTriangular_, lowerTriangular_, or
          // nodeMaxNumRowEntries_) which this loop sets.
          const LO rlcid = colMap.getLocalElement (globalRow);
            // This process owns one or more entries in the current row.
            RowInfo rowInfo = getRowInfo (localRow);
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

              if (smallestCol < localRow) {
                upperTriangular_ = false;
              }
              if (localRow < largestCol) {
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  makeIndicesLocal ()
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "makeIndicesLocal";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! this->hasColMap (), std::logic_error, ": The graph does not have a "
      "column Map yet.  This method should never be called in that case.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      this->getColMap ().is_null (), std::logic_error, ": The graph claims "
      "that it has a column Map, because hasColMap() returns true.  However, "
      "the result of getColMap() is null.  This should never happen.  Please "
      "report this bug to the Tpetra developers.");

    const size_t lclNumRows = this->getNodeNumRows ();
    const map_type& colMap = * (this->getColMap ());

    if (isGloballyIndexed () && lclNumRows > 0) {
      typename t_numRowEntries_::t_host h_numRowEnt = k_numRowEntries_.h_view;

      // allocate data for local indices
      if (getProfileType () == StaticProfile) {
        // If GO and LO are the same size, we can reuse the existing
        // array of 1-D index storage to convert column indices from
        // GO to LO.  Otherwise, we'll just allocate a new buffer.
        if (nodeNumAllocated_ && Kokkos::Impl::is_same<LO,GO>::value) {
          k_lclInds1D_ = Kokkos::Impl::if_c<Kokkos::Impl::is_same<LO,GO>::value,
            t_GlobalOrdinal_1D,
            t_LocalOrdinal_1D >::select (k_gblInds1D_, k_lclInds1D_);
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            k_rowPtrs_.dimension_0 () == 0, std::logic_error, ": This should "
            "never happen at this point.  Please report this bug to the Tpetra "
            "developers.");
          const size_t numEnt = k_rowPtrs_[lclNumRows];

          k_lclInds1D_ = t_LocalOrdinal_1D ("Tpetra::CrsGraph::lclind", numEnt);
        }

        for (size_t r = 0; r < lclNumRows; ++r) {
          const size_t offset   = k_rowPtrs_(r);
          const size_t numentry = h_numRowEnt(r);
          for (size_t j = 0; j < numentry; ++j) {
            const GO gid = k_gblInds1D_(offset + j);
            const LO lid = colMap.getLocalElement (gid);
            k_lclInds1D_(offset + j) = lid;
#ifdef HAVE_TPETRA_DEBUG
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
              k_lclInds1D_(offset + j) == Teuchos::OrdinalTraits<LO>::invalid(),
              std::logic_error,
              ": In local row r=" << r << ", global column " << gid << " is "
              "not in the column Map.  This should never happen.  Please "
              "report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
          }
        }
        // We've converted column indices from global to local, so we
        // can deallocate the global column indices (which we know are
        // in 1-D storage, because the graph has static profile).
        k_gblInds1D_ = t_GlobalOrdinal_1D ();
        gblInds1D_ = Teuchos::null;
      }
      else {  // the graph has dynamic profile (2-D index storage)
        lclInds2D_ = arcp<Array<LO> > (lclNumRows);
        for (size_t r = 0; r < lclNumRows; ++r) {
          if (! gblInds2D_[r].empty ()) {
            const GO* const ginds = gblInds2D_[r].getRawPtr ();
            const size_t rna = gblInds2D_[r].size ();
            const size_t numentry = h_numRowEnt(r);
            lclInds2D_[r].resize (rna);
            LO* const linds = lclInds2D_[r].getRawPtr ();
            for (size_t j = 0; j < numentry; ++j) {
              const GO gid = ginds[j];
              const LO lid = colMap.getLocalElement (gid);
              linds[j] = lid;
#ifdef HAVE_TPETRA_DEBUG
              TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
                linds[j] == Teuchos::OrdinalTraits<LO>::invalid(),
                std::logic_error,
                ": Global column ginds[j=" << j << "]=" << ginds[j]
                << " of local row r=" << r << " is not in the column Map.  "
                "This should never happen.  Please report this bug to the "
                "Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
            }
          }
        }
        gblInds2D_ = Teuchos::null;
      }
    } // globallyIndexed() && lclNumRows > 0

    k_lclGraph_ = LocalStaticCrsGraphType (k_lclInds1D_, k_rowPtrs_);
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  sortAllIndices ()
  {
    // this should be called only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed () );
    if (isSorted () == false) {
      // FIXME (mfh 06 Mar 2014) This would be a good place for a
      // thread-parallel kernel.
      for (size_t row = 0; row < getNodeNumRows (); ++row) {
        RowInfo rowInfo = getRowInfo (row);
        sortRowIndices (rowInfo);
      }
    }
    indicesAreSorted_ = true; // we just sorted every row
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
      const size_t myNumRows = getNodeNumRows ();
      for (size_t r = 0; r < myNumRows; ++r) {
        RowInfo rowinfo = getRowInfo (r);
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  mergeAllIndices ()
  {
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed() ); // call only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( ! isSorted() ); // call only after sortIndices()
    if (! isMerged ()) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = getRowInfo(row);
        mergeRowIndices(rowInfo);
      }
      // we just merged every row
      noRedundancies_ = true;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
              out << "Entries";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              RowInfo rowinfo = getRowInfo(r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << rowinfo.numEntries;
              if (vl == VERB_EXTREME) {
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  checkSizes (const SrcDistObject& source)
  {
    (void) source; // forestall "unused variable" compiler warnings

    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.
    return true;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                    const Teuchos::ArrayView<const GlobalOrdinal> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor& /* distor */,
                    CombineMode /* CM */)
  {
    using Teuchos::ArrayView;

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

    const char tfecfFuncName[] = "unpackAndCombine";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": importLIDs and numPacketsPerLID must have the same size.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete (), std::runtime_error,
      ": Import or Export operations are not allowed on the destination "
      "CrsGraph if it is fill complete.");
    size_t importsOffset = 0;

    typedef typename ArrayView<const LocalOrdinal>::const_iterator iter_type;
    iter_type impLIDiter = importLIDs.begin();
    iter_type impLIDend = importLIDs.end();

    for (size_t i = 0; impLIDiter != impLIDend; ++impLIDiter, ++i) {
      LocalOrdinal row_length = numPacketsPerLID[i];
      ArrayView<const GlobalOrdinal> row (&imports[importsOffset], row_length);
      insertGlobalIndicesFiltered (this->getMap ()->getGlobalElement (*impLIDiter), row);
      importsOffset += row_length;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

} // namespace Tpetra

#endif // TPETRA_CRSGRAPH_DEF_HPP
