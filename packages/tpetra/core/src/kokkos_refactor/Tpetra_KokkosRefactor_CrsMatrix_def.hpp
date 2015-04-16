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

#ifndef TPETRA_KOKKOSREFACTOR_CRSMATRIX_DEF_HPP
#define TPETRA_KOKKOSREFACTOR_CRSMATRIX_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_KokkosRefactor_CrsMatrix_decl.hpp"
#endif
#include <Kokkos_Sequential_SparseKernels.hpp>

namespace Tpetra {

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             size_t maxNumEntriesPerRow,
             ProfileType pftype,
             const RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (pftype == StaticProfile ?
                    Details::STORAGE_1D_UNPACKED :
                    Details::STORAGE_2D),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    try {
      myGraph_ = rcp (new crs_graph_type (rowMap, maxNumEntriesPerRow,
                                          pftype, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
             ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (pftype == StaticProfile ?
                    Details::STORAGE_1D_UNPACKED :
                    Details::STORAGE_2D),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    try {
      myGraph_ = rcp (new Graph (rowMap, NumEntriesPerRowToAlloc, pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             size_t maxNumEntriesPerRow,
             ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (pftype == StaticProfile ?
                    Details::STORAGE_1D_UNPACKED :
                    Details::STORAGE_2D),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    TEUCHOS_TEST_FOR_EXCEPTION(! staticGraph_.is_null(), std::logic_error,
      "Tpetra::CrsMatrix ctor (row Map, col Map, maxNumEntriesPerRow, ...): "
      "staticGraph_ is not null at the beginning of the constructor.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(! myGraph_.is_null(), std::logic_error,
      "Tpetra::CrsMatrix ctor (row Map, col Map, maxNumEntriesPerRow, ...): "
      "myGraph_ is not null at the beginning of the constructor.  "
      "Please report this bug to the Tpetra developers.");
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, maxNumEntriesPerRow,
                                 pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
             ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (pftype == StaticProfile ?
                    Details::STORAGE_1D_UNPACKED :
                    Details::STORAGE_2D),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, numEntPerRow, pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (graph->getRowMap ()),
    staticGraph_ (graph),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    const char tfecfFuncName[] = "CrsMatrix(graph[,params])";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(staticGraph_.is_null (),
      std::runtime_error, ": When calling the CrsMatrix constructor that "
      "accepts a static graph, the pointer to the graph must not be null.");
    // We prohibit the case where the graph is not yet filled.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! staticGraph_->isFillComplete (),
      std::runtime_error, ": The specified graph is not fill-complete. You "
      "must invoke fillComplete() on the graph before using it to construct a "
      "CrsMatrix.  Note that calling resumeFill() makes the graph not fill-"
      "complete, even if you had previously called fillComplete().  In that "
      "case, you must call fillComplete() on the graph again.");
    // the graph has entries, and the matrix should have entries as well, set to zero. no need or point in lazy allocating in this case.
    // first argument LocalIndices is ignored; the graph is already allocated (local or global, we don't care here)
    allocateValues (LocalIndices, GraphAlreadyAllocated);
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const typename local_matrix_type::row_map_type& rowPointers,
             const typename local_graph_type::entries_type::non_const_type& columnIndices,
             const typename local_matrix_type::values_type& values,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, rowPointers,
                                 columnIndices, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    k_values1D_  = values;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::ArrayRCP<size_t> & rowPointers,
             const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices,
             const Teuchos::ArrayRCP<Scalar> & values,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::rcp;
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, rowPointers,
                                 columnIndices, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix constructor: Caught "
        "exception while allocating CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    // FIXME (mfh 05 Aug 2014) It should be possible to convince the
    // ArrayRCP to relinquish its allocation, but that might require
    // passing the ArrayRCP in by nonconst reference.
    Teuchos::ArrayRCP<impl_scalar_type> vals =
      Teuchos::arcp_reinterpret_cast<impl_scalar_type> (values);
    k_values1D_ = Kokkos::Compat::getKokkosViewDeepCopy<DeviceType> (vals ());
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const local_matrix_type& lclMatrix,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    lclMatrix_ (lclMatrix),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (false),
    frobNorm_ (-STM::one ())
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    using Teuchos::rcp;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(rowMap,colMap,lclMatrix,params): ";

    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, lclMatrix.graph, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error, "Caught exception while allocating "
        "CrsGraph: " << e.what ());
    }
    staticGraph_ = myGraph_;
    computeGlobalConstants ();

    k_values1D_ = lclMatrix_.values;

    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    // Now we're fill complete!
    fillComplete_ = true;

    // Sanity checks at the end.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillActive (), std::logic_error,
      "We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillComplete (), std::logic_error,
      "We're at the end of fillComplete(), but isFillComplete() is false.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    checkInternalState ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  ~CrsMatrix ()
  {}

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getComm () const {
    return getCrsGraph ()->getComm ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNode () const {
    return getCrsGraph ()->getNode ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ProfileType
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getProfileType () const {
    return getCrsGraph ()->getProfileType ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isFillComplete () const {
    return fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isFillActive () const {
    return ! fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isStorageOptimized () const {
    return getCrsGraph()->isStorageOptimized();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isLocallyIndexed () const {
    return getCrsGraph ()->isLocallyIndexed ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isGloballyIndexed () const {
    return getCrsGraph ()->isGloballyIndexed ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  hasColMap () const {
    return getCrsGraph ()->hasColMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalNumEntries () const {
    return getCrsGraph ()->getGlobalNumEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNodeNumEntries () const {
    return getCrsGraph ()->getNodeNumEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalNumRows () const {
    return getCrsGraph ()->getGlobalNumRows ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalNumCols () const {
    return getCrsGraph ()->getGlobalNumCols ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNodeNumRows () const {
    return getCrsGraph ()->getNodeNumRows ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNodeNumCols () const {
    return getCrsGraph ()->getNodeNumCols ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,Kokkos
    ::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalNumDiags () const {
    return getCrsGraph ()->getGlobalNumDiags ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNodeNumDiags () const {
    return getCrsGraph ()->getNodeNumDiags ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {
    return getCrsGraph ()->getNumEntriesInGlobalRow (globalRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNumEntriesInLocalRow (LocalOrdinal localRow) const {
    return getCrsGraph ()->getNumEntriesInLocalRow (localRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalMaxNumRowEntries () const {
    return getCrsGraph ()->getGlobalMaxNumRowEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNodeMaxNumRowEntries () const {
    return getCrsGraph ()->getNodeMaxNumRowEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getIndexBase () const {
    return getRowMap ()->getIndexBase ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getRowMap () const {
    return getCrsGraph ()->getRowMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getColMap () const {
    return getCrsGraph ()->getColMap ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<const Map<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getDomainMap () const {
    return getCrsGraph ()->getDomainMap ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<const Map<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getRangeMap () const {
    return getCrsGraph()->getRangeMap();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<const RowGraph<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGraph () const {
    if (staticGraph_ != Teuchos::null) {
      return staticGraph_;
    }
    return myGraph_;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<const CrsGraph<
                 LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getCrsGraph () const {
    if (staticGraph_ != Teuchos::null) {
      return staticGraph_;
    }
    return myGraph_;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isLowerTriangular () const {
    return getCrsGraph ()->isLowerTriangular ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isUpperTriangular () const {
    return getCrsGraph ()->isUpperTriangular ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isStaticGraph () const {
    return myGraph_.is_null ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  hasTransposeApply () const {
    return true;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  supportsRowViews () const {
    return true;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  allocateValues (ELocalGlobal lg, GraphAllocationStatus gas)
  {
#ifdef HAVE_TPETRA_DEBUG
    // If the graph indices are already allocated, then gas should be
    // GraphAlreadyAllocated.  Otherwise, gas should be
    // GraphNotYetAllocated.
    if ((gas == GraphAlreadyAllocated) != staticGraph_->indicesAreAllocated()) {
      const std::string err1 ("allocateValues: The caller has asserted that "
                              "the graph is ");
      const std::string err2 ("already allocated, but the static graph says "
                              "that its indices are ");
      const std::string err3 ("already allocated.  Please report this bug to "
                              "the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(gas == GraphAlreadyAllocated && ! staticGraph_->indicesAreAllocated(),
        std::logic_error, err1 << err2 << "not " << err3);
      TEUCHOS_TEST_FOR_EXCEPTION(gas != GraphAlreadyAllocated && staticGraph_->indicesAreAllocated(),
        std::logic_error, err1 << "not " << err2 << err3);
    }

    // If the graph is unallocated, then it had better be a
    // matrix-owned graph.  ("Matrix-owned graph" means that the
    // matrix gets to define the graph structure.  If the CrsMatrix
    // constructor that takes an RCP<const CrsGraph> was used, then
    // the matrix does _not_ own the graph.)
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! staticGraph_->indicesAreAllocated() && myGraph_.is_null(),
      std::logic_error,
      "allocateValues: The static graph says that its indices are not "
      "allocated, but the graph is not owned by the matrix.  Please report "
      "this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    if (gas == GraphNotYetAllocated) {
      myGraph_->allocateIndices (lg);
    }

    // Allocate matrix values.
    if (getProfileType () == StaticProfile) {
      // "Static profile" means that the number of matrix entries in
      // each row was fixed at the time the CrsMatrix constructor was
      // called.  This lets us use 1-D storage for the matrix's
      // values.  ("1-D storage" means the same as that used by the
      // three arrays in the classic compressed sparse row format.)

      const size_t lclNumRows = staticGraph_->getNodeNumRows ();
      typename Graph::local_graph_type::row_map_type k_ptrs =
        staticGraph_->k_rowPtrs_;
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_ptrs.dimension_0 () != lclNumRows+1, std::logic_error,
        "Tpetra::CrsMatrix::allocateValues: With StaticProfile, row offsets "
        "array has length " << k_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      // FIXME (mfh 08 Aug 2014) This assumes UVM.  We could fix this
      // either by storing the row offsets in the graph as a DualView,
      // or by making a device View of that entry, and copying it back
      // to host.
      const size_t lclTotalNumEntries = k_ptrs(lclNumRows);

      // Allocate array of (packed???) matrix values.
      typedef typename local_matrix_type::values_type values_type;
      k_values1D_ = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);
    }
    else {
      // "Dynamic profile" means the number of matrix entries in each
      // row is not fixed and may expand.  Thus, we store the matrix's
      // values in "2-D storage," meaning an array of arrays.  The
      // outer array has as many inner arrays as there are rows in the
      // matrix, and each inner array stores the values in that row.
      values2D_ = staticGraph_->template allocateValues2D<impl_scalar_type> ();
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers,
                Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices,
                Teuchos::ArrayRCP<const Scalar>& values) const
  {
    const char tfecfFuncName[] = "getAllValues: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      columnIndices.size () != values.size (), std::runtime_error,
      "Requires that columnIndices and values are the same size.");

    RCP<const crs_graph_type> relevantGraph = getCrsGraph ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      relevantGraph.is_null (), std::runtime_error,
      "Requires that getCrsGraph() is not null.");
    try {
      rowPointers = relevantGraph->getNodeRowPtrs ();
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error,
        "Caught exception while calling graph->getNodeRowPtrs(): "
        << e.what ());
    }
    try {
      columnIndices = relevantGraph->getNodePackedIndices ();
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error,
        "Caught exception while calling graph->getNodePackedIndices(): "
        << e.what ());
    }
    Teuchos::ArrayRCP<const impl_scalar_type> vals =
      Kokkos::Compat::persistingView (k_values1D_);
    values = Teuchos::arcp_reinterpret_cast<const Scalar> (vals);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  fillLocalGraphAndMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Kokkos::create_mirror_view;
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef ArrayRCP<size_t>::size_type size_type;
    typedef typename local_matrix_type::row_map_type row_map_type;
    typedef typename Graph::t_numRowEntries_ row_entries_type;
    typedef typename Graph::local_graph_type::entries_type::non_const_type lclinds_1d_type;
    typedef typename local_matrix_type::values_type values_type;

    // fillComplete() only calls fillLocalGraphAndMatrix() if the
    // matrix owns the graph, which means myGraph_ is not null.
    TEUCHOS_TEST_FOR_EXCEPTION(
      myGraph_.is_null (), std::logic_error, "Tpetra::CrsMatrix::"
      "fillLocalGraphAndMatrix (called from fillComplete or "
      "expertStaticFillComplete): The nonconst graph (myGraph_) is null.  This "
      "means that the matrix has a const (a.k.a. \"static\") graph.  This may "
      "mean that fillComplete or expertStaticFillComplete has a bug, since it "
      "should never call fillLocalGraphAndMatrix in that case.  "
      "Please report this bug to the Tpetra developers.");

    const size_t lclNumRows = this->getNodeNumRows ();

    // This method's goal is to fill in the three arrays (compressed
    // sparse row format) that define the sparse graph's and matrix's
    // structure, and the sparse matrix's values.
    //
    // Use the nonconst version of row_map_type for k_ptrs,
    // because row_map_type is const and we need to modify k_ptrs here.
    typename row_map_type::non_const_type k_ptrs;
    row_map_type k_ptrs_const;
    lclinds_1d_type k_inds;
    values_type k_vals;

    // Get references to the data in myGraph_, so we can modify them
    // as well.  Note that we only call fillLocalGraphAndMatrix() if
    // the matrix owns the graph, which means myGraph_ is not null.
    lclinds_1d_type k_lclInds1D_ = myGraph_->k_lclInds1D_;

    // The number of entries in each locally owned row.  This is a
    // DualView.  2-D storage lives on host and is currently not
    // thread-safe for parallel kernels even on host, so we have to
    // work sequentially with host storage in that case.
    row_entries_type k_numRowEnt = myGraph_->k_numRowEntries_;
    typename row_entries_type::t_host h_numRowEnt = k_numRowEnt.h_view;

    if (getProfileType () == DynamicProfile) {
      // Pack 2-D storage (DynamicProfile) into 1-D packed storage.
      //
      // DynamicProfile means that the matrix's column indices and
      // values are currently stored in a 2-D "unpacked" format, in
      // the arrays-of-arrays myGraph_->lclInds2D_ (for column
      // indices) and values2D_ (for values).  We allocate 1-D storage
      // (k_inds resp. k_vals), and then copy from 2-D storage
      // (lclInds2D_ resp. values2D_) into 1-D storage (k_inds
      // resp. k_vals).
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_numRowEnt.dimension_0 ()) != lclNumRows,
        std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix (called "
        "from fillComplete or expertStaticFillComplete): For the "
        "DynamicProfile branch, k_numRowEnt has the wrong length.  "
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
      typename row_map_type::non_const_type::HostMirror h_ptrs;
      {
        // Allocate the packed row offsets array.  We use a nonconst
        // temporary (packedRowOffsets) here, because k_ptrs is const.
        // We will assign packedRowOffsets to k_ptrs below.
        typename row_map_type::non_const_type packedRowOffsets ("Tpetra::CrsGraph::ptr",
                                                                lclNumRows+1);
        //
        // FIXME hack until we get parallel_scan in kokkos
        //
        h_ptrs = create_mirror_view (packedRowOffsets);
        h_ptrs(0) = 0;
        for (size_type i = 0; i < static_cast<size_type> (lclNumRows); ++i) {
          const size_t numEnt = h_numRowEnt(i);
          lclTotalNumEntries += numEnt;
          h_ptrs(i+1) = h_ptrs(i) + numEnt;
        }
        Kokkos::deep_copy (packedRowOffsets, h_ptrs);
        // packedRowOffsets is modifiable; k_ptrs isn't, so we have to
        // use packedRowOffsets in the loop above and assign here.
        k_ptrs = packedRowOffsets;
        k_ptrs_const = k_ptrs;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: In "
        "DynamicProfile branch, after packing k_ptrs, k_ptrs.dimension_0()"
        " = " << k_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (h_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: In "
        "DynamicProfile branch, after packing h_ptrs, h_ptrs.dimension_0()"
        " = " << h_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      // FIXME (mfh 08 Aug 2014) This assumes UVM.
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_ptrs(lclNumRows) != lclTotalNumEntries, std::logic_error,
        "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: In DynamicProfile branch, "
        "after packing k_ptrs, k_ptrs(lclNumRows = " << lclNumRows << ") = " <<
        k_ptrs(lclNumRows) << " != total number of entries on the calling "
        "process = " << lclTotalNumEntries << ".");

      // Allocate the arrays of packed column indices and values.
      k_inds = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);
      k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);

      // We need host views of the above, since 2-D storage lives on host.
      typename lclinds_1d_type::HostMirror h_inds = create_mirror_view (k_inds);
      typename values_type::HostMirror h_vals = create_mirror_view (k_vals);

      // Pack the column indices and values on the host.
      ArrayRCP<Array<LocalOrdinal> > lclInds2D = myGraph_->lclInds2D_;
      for (size_t row = 0; row < lclNumRows; ++row) {
        const size_t numEnt = h_numRowEnt(row);
        std::copy (lclInds2D[row].begin(),
                   lclInds2D[row].begin() + numEnt,
                   h_inds.ptr_on_device() + h_ptrs(row));
        std::copy (values2D_[row].begin(),
                   values2D_[row].begin() + numEnt,
                   h_vals.ptr_on_device() + h_ptrs(row));
      }
      // Copy the packed column indices and values to the device.
      Kokkos::deep_copy (k_inds, h_inds);
      Kokkos::deep_copy (k_vals, h_vals);

      // Sanity check of packed row offsets.
      if (k_ptrs.dimension_0 () != 0) {
        const size_t numOffsets = static_cast<size_t> (k_ptrs.dimension_0 ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs(numOffsets-1)) != k_vals.dimension_0 (),
          std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
          "In DynamicProfile branch, after packing, k_ptrs(" << (numOffsets-1)
          << ") = " << k_ptrs(numOffsets-1) << " != k_vals.dimension_0() = "
          << k_vals.dimension_0 () << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs(numOffsets-1)) != k_inds.dimension_0 (),
          std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
          "In DynamicProfile branch, after packing, k_ptrs(" << (numOffsets-1)
          << ") = " << k_ptrs(numOffsets-1) << " != k_inds.dimension_0() = "
          << k_inds.dimension_0 () << ".");
      }
    }
    else if (getProfileType () == StaticProfile) {
      // StaticProfile means that the matrix's column indices and
      // values are currently stored in a 1-D format, with row offsets
      // in k_rowPtrs_ and local column indices in k_lclInds1D_.

      // StaticProfile also means that the graph's array of row
      // offsets must already be allocated.
      typename Graph::local_graph_type::row_map_type curRowOffsets =
        myGraph_->k_rowPtrs_;
      TEUCHOS_TEST_FOR_EXCEPTION(
        curRowOffsets.dimension_0 () == 0, std::logic_error,
      "curRowOffsets has size zero, but shouldn't");
      TEUCHOS_TEST_FOR_EXCEPTION(
        curRowOffsets.dimension_0 () != lclNumRows + 1, std::logic_error,
        "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: curRowOffsets has size "
        << curRowOffsets.dimension_0 () << " != lclNumRows + 1 = "
        << (lclNumRows + 1) << ".")
      {
        const size_t numOffsets = curRowOffsets.dimension_0 ();
        // FIXME (mfh 06 Aug 2014) This relies on UVM.
        TEUCHOS_TEST_FOR_EXCEPTION(
          numOffsets != 0 &&
          myGraph_->k_lclInds1D_.dimension_0 () != curRowOffsets(numOffsets - 1),
          std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
          "numOffsets = " << numOffsets << " != 0 and "
          "myGraph_->k_lclInds1D_.dimension_0() = "
          << myGraph_->k_lclInds1D_.dimension_0 ()
          << " != curRowOffsets(" << numOffsets << ") = "
          << curRowOffsets(numOffsets - 1) << ".");
      }

      if (myGraph_->nodeNumEntries_ != myGraph_->nodeNumAllocated_) {
        // The matrix's current 1-D storage is "unpacked."  This means
        // the row offsets may differ from what the final row offsets
        // should be.  This could happen, for example, if the user
        // specified StaticProfile in the constructor and set an upper
        // bound on the number of entries per row, but didn't fill all
        // those entries.
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_numRowEnt.dimension_0 ()) != lclNumRows,
          std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix (called"
          " from fillComplete or expertStaticFillComplete): In StaticProfile "
          "unpacked branch, k_numRowEnt has the wrong length.  "
          "k_numRowEnt.dimension_0() = " << k_numRowEnt.dimension_0 ()
          << " != getNodeNumRows() = " << lclNumRows << ".");

        if (curRowOffsets.dimension_0 () != 0) {
          const size_t numOffsets =
            static_cast<size_t> (curRowOffsets.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            curRowOffsets(numOffsets-1) != static_cast<size_t> (k_values1D_.dimension_0 ()),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile branch, before allocating or packing, "
            "curRowOffsets(" << (numOffsets-1) << ") = "
            << curRowOffsets(numOffsets - 1)
            << " != k_values1D_.dimension_0() = "
            << k_values1D_.dimension_0 () << ".");
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (curRowOffsets(numOffsets - 1)) !=
            myGraph_->k_lclInds1D_.dimension_0 (),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile branch, before allocating or packing, "
            "curRowOffsets(" << (numOffsets-1) << ") = "
            << curRowOffsets(numOffsets - 1)
            << " != myGraph_->k_lclInds1D_.dimension_0() = "
            << myGraph_->k_lclInds1D_.dimension_0 () << ".");
        }

        // Pack the row offsets into k_ptrs, by doing a sum-scan of
        // the array of valid entry counts per row (h_numRowEnt).

        // Total number of entries in the matrix on the calling
        // process.  We will compute this in the loop below.  It's
        // cheap to compute and useful as a sanity check.
        size_t lclTotalNumEntries = 0;
        // This will be a host view of packed row offsets.
        typename row_map_type::non_const_type::HostMirror h_ptrs;
        {
          // Allocate the packed row offsets array.  We use a nonconst
          // temporary (packedRowOffsets) here, because k_ptrs is
          // const.  We will assign packedRowOffsets to k_ptrs below.
          typename row_map_type::non_const_type packedRowOffsets ("Tpetra::CrsGraph::ptr",
                                                                  lclNumRows+1);
          //
          // FIXME hack until we get parallel_scan in Kokkos
          //
          // Unlike in the 2-D storage case above, we don't need the
          // host view of the packed row offsets array after packing
          // the row offsets.
          h_ptrs = create_mirror_view (packedRowOffsets);
          h_ptrs(0) = 0;
          for (size_type i = 0; i < static_cast<size_type> (lclNumRows); ++i) {
            const size_t numEnt = h_numRowEnt(i);
            lclTotalNumEntries += numEnt;
            h_ptrs(i+1) = h_ptrs(i) + numEnt;
          }
          Kokkos::deep_copy (packedRowOffsets, h_ptrs);
          // packedRowOffsets is modifiable; k_ptrs isn't, so we have
          // to use packedRowOffsets in the loop above and assign here.
          k_ptrs = packedRowOffsets;
          k_ptrs_const = k_ptrs;
        }

        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs.dimension_0 ()) != lclNumRows + 1,
          std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: For "
          "the StaticProfile unpacked-but-pack branch, after packing k_ptrs, "
          "k_ptrs.dimension_0() = " << k_ptrs.dimension_0 () << " != "
          "lclNumRows+1 = " << (lclNumRows+1) << ".");
        // FIXME (mfh 06 Aug 2014) This assumes UVM.
        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs(lclNumRows) != lclTotalNumEntries, std::logic_error,
          "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: In StaticProfile "
          "unpacked-but-pack branch, after filling k_ptrs, k_ptrs(lclNumRows="
          << lclNumRows << ") = " << k_ptrs(lclNumRows) << " != total number "
          "of entries on the calling process = " << lclTotalNumEntries << ".");

        // Allocate the arrays of packed column indices and values.
        k_inds = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);
        k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);

        // curRowOffsets (myGraph_->k_rowPtrs_) (???), k_lclInds1D_,
        // and k_values1D_ are currently unpacked.  Pack them, using
        // the packed row offsets array k_ptrs that we created above.
        //
        // FIXME (mfh 06 Aug 2014) If "Optimize Storage" is false, we
        // need to keep around the unpacked row offsets, column
        // indices, and values arrays.

        // Pack the column indices from unpacked k_lclInds1D_ into
        // packed k_inds.  We will replace k_lclInds1D_ below.
        typedef pack_functor<typename Graph::local_graph_type::entries_type::non_const_type,
          typename Graph::local_graph_type::row_map_type>
          inds_packer_type;
        inds_packer_type indsPacker (k_inds, myGraph_->k_lclInds1D_,
                                     k_ptrs, curRowOffsets);
        Kokkos::parallel_for (lclNumRows, indsPacker);

        // Pack the values from unpacked k_values1D_ into packed
        // k_vals.  We will replace k_values1D_ below.
        typedef pack_functor<values_type, row_map_type> vals_packer_type;
        vals_packer_type valsPacker (k_vals, this->k_values1D_,
                                     k_ptrs, curRowOffsets);
        Kokkos::parallel_for (lclNumRows, valsPacker);

        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs.dimension_0 () == 0, std::logic_error, "Tpetra::CrsMatrix::"
          "fillLocalGraphAndMatrix: In StaticProfile \"Optimize Storage\" = "
          "true branch, after packing, k_ptrs.dimension_0() = 0.  This "
          "probably means that k_rowPtrs_ was never allocated.");
        if (k_ptrs.dimension_0 () != 0) {
          const size_t numOffsets = static_cast<size_t> (k_ptrs.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs(numOffsets - 1)) != k_vals.dimension_0 (),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile \"Optimize Storage\"=true branch, after packing, "
            "k_ptrs(" << (numOffsets-1) << ") = " << k_ptrs(numOffsets-1) <<
            " != k_vals.dimension_0() = " << k_vals.dimension_0 () << ".");
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs(numOffsets - 1)) != k_inds.dimension_0 (),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile \"Optimize Storage\"=true branch, after packing, "
            "k_ptrs(" << (numOffsets-1) << ") = " << k_ptrs(numOffsets-1) <<
            " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
        }
      }
      else { // We don't have to pack, so just set the pointers.
        k_ptrs_const = myGraph_->k_rowPtrs_;
        k_inds = myGraph_->k_lclInds1D_;
        k_vals = this->k_values1D_;

        TEUCHOS_TEST_FOR_EXCEPTION(
          k_ptrs_const.dimension_0 () == 0, std::logic_error, "Tpetra::CrsMatrix::"
          "fillLocalGraphAndMatrix: In StaticProfile \"Optimize Storage\" = "
          "false branch, k_ptrs_const.dimension_0() = 0.  This probably means that "
          "k_rowPtrs_ was never allocated.");
        if (k_ptrs_const.dimension_0 () != 0) {
          const size_t numOffsets = static_cast<size_t> (k_ptrs_const.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_vals.dimension_0 (),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile \"Optimize Storage\" = false branch, "
            "k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets - 1)
            << " != k_vals.dimension_0() = " << k_vals.dimension_0 () << ".");
          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_inds.dimension_0 (),
            std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: "
            "In StaticProfile \"Optimize Storage\" = false branch, "
            "k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets - 1)
            << " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
        }
      }
    }

    // Extra sanity checks.
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (k_ptrs_const.dimension_0 ()) != lclNumRows + 1,
      std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: After "
      "packing, k_ptrs_const.dimension_0() = " << k_ptrs_const.dimension_0 ()
      << " != lclNumRows+1 = " << (lclNumRows+1) << ".");
    if (k_ptrs_const.dimension_0 () != 0) {
      const size_t numOffsets = static_cast<size_t> (k_ptrs_const.dimension_0 ());
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_vals.dimension_0 (),
        std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: After "
        "packing, k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets-1)
        << " != k_vals.dimension_0() = " << k_vals.dimension_0 () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs_const(numOffsets - 1)) != k_inds.dimension_0 (),
        std::logic_error, "Tpetra::CrsMatrix::fillLocalGraphAndMatrix: After "
        "packing, k_ptrs_const(" << (numOffsets-1) << ") = " << k_ptrs_const(numOffsets-1)
        << " != k_inds.dimension_0() = " << k_inds.dimension_0 () << ".");
    }

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Optimize
    // storage if the graph is not static, or if the graph already has
    // optimized storage.
    const bool defaultOptStorage =
      ! isStaticGraph () || staticGraph_->isStorageOptimized ();
    const bool requestOptimizedStorage =
      (! params.is_null () && params->get ("Optimize Storage", defaultOptStorage)) ||
      (params.is_null () && defaultOptStorage);

    // The graph has optimized storage when indices are allocated,
    // myGraph_->k_numRowEntries_ is empty, and there are more than
    // zero rows on this process.  It's impossible for the graph to
    // have dynamic profile (getProfileType() == DynamicProfile) and
    // be optimized (isStorageOptimized()).
    if (requestOptimizedStorage) {
      // Free the old, unpacked, unoptimized allocations.
      // Change the graph from dynamic to static allocation profile

      // Free graph data structures that are only needed for 2-D or
      // unpacked 1-D storage.
      myGraph_->lclInds2D_ = null; // legacy KokkosClassic 2-D storage
      myGraph_->k_numRowEntries_ = row_entries_type ();

      // Free the matrix's 2-D storage.
      this->values2D_ = null;

      // Keep the new 1-D packed allocations.
      myGraph_->k_rowPtrs_ = k_ptrs_const;
      myGraph_->k_lclInds1D_ = k_inds;
      this->k_values1D_ = k_vals;

      // Storage is packed now, so the number of allocated entries is
      // the same as the actual number of entries.
      myGraph_->nodeNumAllocated_ = myGraph_->nodeNumEntries_;
      // The graph is definitely StaticProfile now, whether or not it
      // was before.
      myGraph_->pftype_ = StaticProfile;
      myGraph_->storageStatus_ = Details::STORAGE_1D_PACKED;
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }

    // Make the local graph, using the arrays of row offsets and
    // column indices that we built above.  The local graph should be
    // null, but we delete it first so that any memory can be freed
    // before we allocate the new one.
    //
    // FIXME (mfh 06,28 Aug 2014) It would make more sense for
    // Tpetra::CrsGraph to have a protected method that accepts k_inds
    // and k_ptrs, and creates the local graph lclGraph_.
    myGraph_->lclGraph_ =
      typename Graph::local_graph_type (k_inds, k_ptrs_const);

    // Make the local matrix, using the local graph and vals array.
    lclMatrix_ = local_matrix_type ("Tpetra::CrsMatrix::lclMatrix_",
                                    getNodeNumCols (), k_vals,
                                    myGraph_->lclGraph_);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  fillLocalMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Kokkos::create_mirror_view;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef LocalOrdinal LO;
    typedef typename Graph::t_numRowEntries_ row_entries_type;
    typedef typename Graph::local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::non_const_type non_const_row_map_type;
    typedef typename local_matrix_type::values_type values_type;

    const size_t lclNumRows = getNodeNumRows();
    const map_type& rowMap = * (getRowMap ());
    RCP<node_type> node = rowMap.getNode ();

    // The goals of this routine are first, to allocate and fill
    // packed 1-D storage (see below for an explanation) in the vals
    // array, and second, to give vals to the local matrix and
    // finalize the local matrix.  We only need k_ptrs, the packed 1-D
    // row offsets, within the scope of this routine, since we're only
    // filling the local matrix here (use fillLocalGraphAndMatrix() to
    // fill both the graph and the matrix at the same time).

    // get data from staticGraph_
    ArrayRCP<Array<LO> > lclInds2D = staticGraph_->lclInds2D_;
    size_t nodeNumEntries   = staticGraph_->nodeNumEntries_;
    size_t nodeNumAllocated = staticGraph_->nodeNumAllocated_;
    row_map_type k_rowPtrs_ = staticGraph_->lclGraph_.row_map;

    row_map_type k_ptrs; // "packed" row offsets array
    values_type k_vals; // "packed" values array

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Request
    // optimized storage by default.
    bool requestOptimizedStorage = true;
    const bool default_OptimizeStorage =
      ! isStaticGraph () || staticGraph_->isStorageOptimized ();
    if (! params.is_null () && ! params->get ("Optimize Storage", default_OptimizeStorage)) {
      requestOptimizedStorage = false;
    }
    // If we're not allowed to change a static graph, then we can't
    // change the storage of the matrix, either.  This means that if
    // the graph's storage isn't already optimized, we can't optimize
    // the matrix's storage either.  Check and give warning, as
    // appropriate.
    if (! staticGraph_->isStorageOptimized () && requestOptimizedStorage) {
      TPETRA_ABUSE_WARNING(true, std::runtime_error,
        "::fillLocalMatrix(): You requested optimized storage by setting the"
        "\"Optimize Storage\" flag to \"true\" in the parameter list, or by virtue"
        "of default behavior. However, the associated CrsGraph was filled separately"
        "and requested not to optimize storage. Therefore, the CrsMatrix cannot"
        "optimize storage.");
      requestOptimizedStorage = false;
    }

    // The number of entries in each locally owned row.  This is a
    // DualView.  2-D storage lives on host and is currently not
    // thread-safe for parallel kernels even on host, so we have to
    // work sequentially with host storage in that case.
    row_entries_type k_numRowEnt = staticGraph_->k_numRowEntries_;
    typename row_entries_type::t_host h_numRowEnt = k_numRowEnt.h_view;

    if (getProfileType() == DynamicProfile) {
      // Pack 2-D storage (DynamicProfile) into 1-D packed storage.
      //
      // DynamicProfile means that the matrix's values are currently
      // stored in a 2-D "unpacked" format, in the array-of-arrays
      // values2D_.  We allocate 1-D storage and then copy from 2-D
      // storage in values2D_ into 1-D storage in k_vals.  Since we're
      // only allocating the local matrix here, not the local graph,
      // we don't need to keep the row offsets array, but we do need
      // it here temporarily in order to convert to 1-D storage.  (The
      // allocStorage() function needs it.)  We'll free ptrs later in
      // this method.
      //
      // FIXME (mfh 08 Aug 2014) If we're in this method, then the
      // graph should already have packed 1-D storage.  Why can't we
      // just use the graph's current row offsets array?

      // Pack the row offsets into k_ptrs, by doing a sum-scan of
      // the array of valid entry counts per row (h_numRowEnt).
      //
      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      // This will be a host view of packed row offsets.
      typename non_const_row_map_type::HostMirror h_ptrs;
      {
        non_const_row_map_type packedRowOffsets ("Tpetra::CrsGraph::ptr",
                                                 lclNumRows+1);
        //
        // FIXME hack until we get parallel_scan in Kokkos
        //
        h_ptrs = create_mirror_view (packedRowOffsets);
        h_ptrs(0) = 0;
        for (size_t i = 0; i < lclNumRows; ++i) {
          const size_t numEnt = h_numRowEnt(i);
          lclTotalNumEntries += numEnt;
          h_ptrs(i+1) = h_ptrs(i) + numEnt;
        }
        Kokkos::deep_copy (packedRowOffsets, h_ptrs);
        k_ptrs = packedRowOffsets;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (k_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsMatrix::fillLocalMatrix: In "
        "DynamicProfile branch, after packing k_ptrs, k_ptrs.dimension_0()"
        " = " << k_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (h_ptrs.dimension_0 ()) != lclNumRows + 1,
        std::logic_error, "Tpetra::CrsMatrix::fillLocalMatrix: In "
        "DynamicProfile branch, after packing h_ptrs, h_ptrs.dimension_0()"
        " = " << h_ptrs.dimension_0 () << " != (lclNumRows+1) = "
        << (lclNumRows+1) << ".");
      // FIXME (mfh 08 Aug 2014) This assumes UVM.
      TEUCHOS_TEST_FOR_EXCEPTION(
        k_ptrs(lclNumRows) != lclTotalNumEntries, std::logic_error,
        "Tpetra::CrsMatrix::fillLocalMatrix: In DynamicProfile branch, "
        "after packing k_ptrs, k_ptrs(lclNumRows = " << lclNumRows << ") = " <<
        k_ptrs(lclNumRows) << " != total number of entries on the calling "
        "process = " << lclTotalNumEntries << ".");

      // Allocate the array of packed values.
      k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);
      // We need a host view of the above, since 2-D storage lives on host.
      typename values_type::HostMirror h_vals = create_mirror_view (k_vals);
      // Pack the values on the host.
      for (size_t lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const size_t numEnt = h_numRowEnt(lclRow);
        std::copy (values2D_[lclRow].begin(),
                   values2D_[lclRow].begin() + numEnt,
                   h_vals.ptr_on_device() + h_ptrs(lclRow));
      }
      // Copy the packed values to the device.
      Kokkos::deep_copy (k_vals, h_vals);

      // Sanity check of packed row offsets.
      if (k_ptrs.dimension_0 () != 0) {
        const size_t numOffsets = static_cast<size_t> (k_ptrs.dimension_0 ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          static_cast<size_t> (k_ptrs(numOffsets-1)) != k_vals.dimension_0 (),
          std::logic_error, "Tpetra::CrsMatrix::fillLocalMatrix: "
          "In DynamicProfile branch, after packing, k_ptrs(" << (numOffsets-1)
          << ") = " << k_ptrs(numOffsets-1) << " != k_vals.dimension_0() = "
          << k_vals.dimension_0 () << ".");
      }
    }
    else if (getProfileType () == StaticProfile) {
      // StaticProfile means that the matrix's values are currently
      // stored in a 1-D format.  However, this format is "unpacked";
      // it doesn't necessarily have the same row offsets as indicated
      // by the ptrs array returned by allocRowPtrs.  This could
      // happen, for example, if the user specified StaticProfile in
      // the constructor and fixed the number of matrix entries in
      // each row, but didn't fill all those entries.
      //
      // As above, we don't need to keep the "packed" row offsets
      // array ptrs here, but we do need it here temporarily, so we
      // have to allocate it.  We'll free ptrs later in this method.
      //
      // Note that this routine checks whether storage has already
      // been packed.  This is a common case for solution of nonlinear
      // PDEs using the finite element method, as long as the
      // structure of the sparse matrix does not change between linear
      // solves.
      if (nodeNumEntries != nodeNumAllocated) {
        // We have to pack the 1-D storage, since the user didn't fill
        // up all requested storage.
        non_const_row_map_type tmpk_ptrs ("Tpetra::CrsGraph::ptr",
                                          lclNumRows+1);
        // Total number of entries in the matrix on the calling
        // process.  We will compute this in the loop below.  It's
        // cheap to compute and useful as a sanity check.
        size_t lclTotalNumEntries = 0;
        k_ptrs = tmpk_ptrs;
        {
          //
          // FIXME hack until we get parallel_scan in Kokkos
          //
          typename non_const_row_map_type::HostMirror h_ptrs =
            create_mirror_view (tmpk_ptrs);
          h_ptrs(0) = 0;
          for (size_t i = 0; i < lclNumRows; ++i) {
            const size_t numEnt = h_numRowEnt(i);
            lclTotalNumEntries += numEnt;
            h_ptrs(i+1) = h_ptrs(i) + numEnt;
          }
          Kokkos::deep_copy (tmpk_ptrs, h_ptrs);
        }

        // Allocate the "packed" values array.
        // It has exactly the right number of entries.
        k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);

        // Pack k_values1D_ into k_vals.  We will replace k_values1D_ below.
        typedef pack_functor<values_type, row_map_type> packer_type;
        packer_type valsPacker (k_vals, k_values1D_, tmpk_ptrs, k_rowPtrs_);
        Kokkos::parallel_for (lclNumRows, valsPacker);
      }
      else { // We don't have to pack, so just set the pointer.
        k_vals = k_values1D_;
      }
    }

    // May we ditch the old allocations for the packed one?
    if (requestOptimizedStorage) {
      // The user requested optimized storage, so we can dump the
      // unpacked 2-D and 1-D storage, and keep the packed storage.
      values2D_ = null;
      k_values1D_ = k_vals;
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }

    // Build the local sparse matrix object.
    lclMatrix_ = local_matrix_type ("Tpetra::CrsMatrix::lclMatrix_",
                                    getDomainMap ()->getNodeNumElements (),
                                    k_vals,
                                    staticGraph_->getLocalGraph ());
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  insertLocalValues (const LocalOrdinal localRow,
                     const Teuchos::ArrayView<const LocalOrdinal>& indices,
                     const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::toString;
    using std::endl;
    const char tfecfFuncName[] = "insertLocalValues";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive (), std::runtime_error,
      ": Fill is not active.  After calling fillComplete, you must call "
      "resumeFill before you may insert entries into the matrix again.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isStaticGraph (),  std::runtime_error,
      " cannot insert indices with static graph; use replaceLocalValues() instead.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(myGraph_->isGloballyIndexed(),
      std::runtime_error, ": graph indices are global; use insertGlobalValues().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! hasColMap (), std::runtime_error,
      " cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
      std::runtime_error, ": values.size() must equal indices.size().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! getRowMap()->isNodeLocalElement(localRow), std::runtime_error,
      ": Local row index " << localRow << " does not belong to this process.");

    if (! myGraph_->indicesAreAllocated ()) {
      try {
        allocateValues (LocalIndices, GraphNotYetAllocated);
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Tpetra::CrsMatrix::insertLocalValues: "
          "allocateValues(LocalIndices,GraphNotYetAllocated) threw an "
          "exception: " << e.what ());
      }
    }

    const size_t numEntriesToAdd = static_cast<size_t> (indices.size ());
#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, if the matrix has a column Map, test whether
    // any of the given column indices are not in the column Map.
    // Keep track of the invalid column indices so we can tell the
    // user about them.
    if (hasColMap ()) {
      const map_type& colMap = * (getColMap ());
      Array<LocalOrdinal> badColInds;
      bool allInColMap = true;
      for (size_t k = 0; k < numEntriesToAdd; ++k) {
        if (! colMap.isNodeLocalElement (indices[k])) {
          allInColMap = false;
          badColInds.push_back (indices[k]);
        }
      }
      if (! allInColMap) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::insertLocalValues: You attempted to insert "
          "entries in owned row " << localRow << ", at the following column "
          "indices: " << toString (indices) << "." << endl;
        os << "Of those, the following indices are not in the column Map on "
          "this process: " << toString (badColInds) << "." << endl << "Since "
          "the matrix has a column Map already, it is invalid to insert "
          "entries at those locations.";
        TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
      }
    }
#endif // HAVE_TPETRA_DEBUG

#ifdef HAVE_TPETRA_DEBUG
    RowInfo rowInfo;
    try {
      rowInfo = myGraph_->getRowInfo (localRow);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix::insertLocalValues: "
        "myGraph_->getRowInfo threw an exception: " << e.what ());
    }
#else
    RowInfo rowInfo = myGraph_->getRowInfo (localRow);
#endif // HAVE_TPETRA_DEBUG

    const size_t curNumEntries = rowInfo.numEntries;
    const size_t newNumEntries = curNumEntries + numEntriesToAdd;
    if (newNumEntries > rowInfo.allocSize) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() == StaticProfile, std::runtime_error,
        ": new indices exceed statically allocated graph structure.");

      // Make space for the new matrix entries.
      try {
        rowInfo = myGraph_->template updateLocalAllocAndValues<impl_scalar_type> (rowInfo,
                                                                             newNumEntries,
                                                                             values2D_[localRow]);
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Tpetra::CrsMatrix::insertLocalValues: "
          "myGraph_->updateGlobalAllocAndValues threw an exception: "
          << e.what ());
      }
    }
    typename Graph::SLocalGlobalViews indsView;
    indsView.linds = indices;

#ifdef HAVE_TPETRA_DEBUG
    ArrayView<impl_scalar_type> valsView;
    try {
      valsView = this->getViewNonConst (rowInfo);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix::insertLocalValues: "
        "getViewNonConst threw an exception: " << e.what ());
    }
#else
    ArrayView<impl_scalar_type> valsView = this->getViewNonConst (rowInfo);
#endif // HAVE_TPETRA_DEBUG

    ArrayView<const impl_scalar_type> valsIn =
      av_reinterpret_cast<const impl_scalar_type> (values);
    try {
      myGraph_->template insertIndicesAndValues<impl_scalar_type> (rowInfo, indsView,
                                                              valsView, valsIn,
                                                              LocalIndices,
                                                              LocalIndices);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Tpetra::CrsMatrix::insertLocalValues: "
        "myGraph_->insertIndicesAndValues threw an exception: "
        << e.what ());
    }

#ifdef HAVE_TPETRA_DEBUG
    const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow (localRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      chkNewNumEntries != newNumEntries, std::logic_error,
      ": The row should have " << newNumEntries << " entries after insert, but "
      "instead has " << chkNewNumEntries << ".  Please report this bug to the "
      "Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isLocallyIndexed(), std::logic_error,
      ": At end of insertLocalValues(), this CrsMatrix is not locally indexed.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  insertLocalValuesFiltered (const LocalOrdinal localRow,
                             const Teuchos::ArrayView<const LocalOrdinal>& indices,
                             const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    const char tfecfFuncName[] = "insertLocalValues: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive (), std::runtime_error,
      "Requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isStaticGraph (),  std::runtime_error,
      "Cannot insert indices with static graph; use replaceLocalValues() instead.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(myGraph_->isGloballyIndexed(),
      std::runtime_error, "Graph indices are global; use insertGlobalValues().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error, "The matrix has no column Map yet, "
      "so you cannot insert local indices.  If you created the matrix without "
      "a column Map (or without a fill-complete graph), you must call "
      "fillComplete to create the column Map, before you may work with local "
      "indices.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      values.size () != indices.size (), std::runtime_error, "values.size() = "
      << values.size () << " != indices.size() = " << indices.size ()<< ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! getRowMap()->isNodeLocalElement (localRow), std::runtime_error,
      "Local row index " << localRow << " does not belong to this process.");
    if (! myGraph_->indicesAreAllocated ()) {
      allocateValues (LocalIndices, GraphNotYetAllocated);
    }
    // Use the graph to filter incoming entries whose column indices
    // aren't in the column Map.
    Array<LocalOrdinal> f_inds (indices);
    ArrayView<const impl_scalar_type> valsIn =
      av_reinterpret_cast<const impl_scalar_type> (values);
    Array<impl_scalar_type> f_vals (valsIn);
    const size_t numFilteredEntries =
      myGraph_->template filterLocalIndicesAndValues<impl_scalar_type> (f_inds (),
                                                                   f_vals ());
    if (numFilteredEntries > 0) {
      RowInfo rowInfo = myGraph_->getRowInfo (localRow);
      const size_t curNumEntries = rowInfo.numEntries;
      const size_t newNumEntries = curNumEntries + numFilteredEntries;
      if (newNumEntries > rowInfo.allocSize) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          getProfileType () == StaticProfile, std::runtime_error,
          ": new indices exceed statically allocated graph structure.  "
          "newNumEntries (" << newNumEntries << " > rowInfo.allocSize ("
          << rowInfo.allocSize << ").");
        // Make space for the new matrix entries.
        rowInfo =
          myGraph_->template updateLocalAllocAndValues<impl_scalar_type> (rowInfo,
                                                                     newNumEntries,
                                                                     values2D_[localRow]);
      }
      typename Graph::SLocalGlobalViews inds_view;
      inds_view.linds = f_inds (0, numFilteredEntries);
      myGraph_->template insertIndicesAndValues<impl_scalar_type> (rowInfo, inds_view,
                                                              this->getViewNonConst (rowInfo),
                                                              f_vals, LocalIndices,
                                                              LocalIndices);
#ifdef HAVE_TPETRA_DEBUG
      const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow (localRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries,
        std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isLocallyIndexed(), std::logic_error,
      ": At end of insertLocalValues(), this CrsMatrix is not locally indexed.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  insertGlobalValues (const GlobalOrdinal globalRow,
                      const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                      const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::toString;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;
    const char tfecfFuncName[] = "insertGlobalValues: ";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      values.size () != indices.size (), std::runtime_error,
      "values.size() = " << values.size() << " != indices.size() = "
      << indices.size() << ".");
#endif // HAVE_TPETRA_DEBUG

    const LO localRow = getRowMap ()->getLocalElement (globalRow);

    if (localRow == OTL::invalid ()) { // globalRow _not_ owned by calling process
      insertNonownedGlobalValues (globalRow, indices, values);
    }
    else { // globalRow _is_ owned by calling process
      if (this->isStaticGraph ()) {
        // Uh oh!  Not allowed to insert into owned rows in that case.
        std::ostringstream err;
        const int myRank = getRowMap ()->getComm ()->getRank ();
        const int numProcs = getRowMap ()->getComm ()->getSize ();

        err << "The matrix was constructed with a constant (\"static\") graph, "
          "yet the given global row index " << globalRow << " is in the row "
          "Map on the calling process (with rank " << myRank << ", of " <<
          numProcs << " process(es)).  In this case, you may not insert new "
          "entries into rows owned by the calling process.";

        if (! getRowMap ()->isNodeGlobalElement (globalRow)) {
          err << "  Furthermore, GID->LID conversion with the row Map claims that "
            "the global row index is owned on the calling process, yet "
            "getRowMap()->isNodeGlobalElement(globalRow) returns false.  That's"
            " weird!  This might indicate a Map bug.  Please report this to the"
            " Tpetra developers.";
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          this->isStaticGraph (), std::runtime_error, err.str ());
      }

      if (! myGraph_->indicesAreAllocated ()) {
        try {
          allocateValues (GlobalIndices, GraphNotYetAllocated);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::runtime_error, "Tpetra::CrsMatrix::insertGlobalValues: "
            "allocateValues(GlobalIndices,GraphNotYetAllocated) threw an "
            "exception: " << e.what ());
        }
      }

      const size_type numEntriesToInsert = indices.size ();
      // If the matrix has a column Map, check at this point whether
      // the column indices belong to the column Map.
      //
      // FIXME (mfh 16 May 2013) We may want to consider deferring the
      // test to the CrsGraph method, since it may have to do this
      // anyway.
      if (hasColMap ()) {
        const map_type& colMap = * (getColMap ());
        // In a debug build, keep track of the nonowned ("bad") column
        // indices, so that we can display them in the exception
        // message.  In a release build, just ditch the loop early if
        // we encounter a nonowned column index.
#ifdef HAVE_TPETRA_DEBUG
        Array<GO> badColInds;
#endif // HAVE_TPETRA_DEBUG
        bool allInColMap = true;
        for (size_type k = 0; k < numEntriesToInsert; ++k) {
          if (! colMap.isNodeGlobalElement (indices[k])) {
            allInColMap = false;
#ifdef HAVE_TPETRA_DEBUG
            badColInds.push_back (indices[k]);
#else
            break;
#endif // HAVE_TPETRA_DEBUG
          }
        }
        if (! allInColMap) {
          std::ostringstream os;
          os << "You attempted to insert entries in owned row " << globalRow
             << ", at the following column indices: " << toString (indices)
             << "." << endl;
#ifdef HAVE_TPETRA_DEBUG
          os << "Of those, the following indices are not in the column Map on "
            "this process: " << toString (badColInds) << "." << endl << "Since "
            "the matrix has a column Map already, it is invalid to insert "
            "entries at those locations.";
#else
          os << "At least one of those indices is not in the column Map on this "
            "process." << endl << "It is invalid to insert into columns not in "
            "the column Map on the process that owns the row.";
#endif // HAVE_TPETRA_DEBUG
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            ! allInColMap, std::invalid_argument, os.str ());
        }
      }

      typename Graph::SLocalGlobalViews inds_view;
      ArrayView<const impl_scalar_type> vals_view;

      inds_view.ginds = indices;
      vals_view       = av_reinterpret_cast<const impl_scalar_type> (values);

#ifdef HAVE_TPETRA_DEBUG
      RowInfo rowInfo;
      try {
        rowInfo = myGraph_->getRowInfo (localRow);
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "myGraph_->getRowInfo(localRow=" << localRow
          << ") threw an exception: " << e.what ());
      }
#else
      RowInfo rowInfo = myGraph_->getRowInfo (localRow);
#endif // HAVE_TPETRA_DEBUG

      const size_t curNumEntries = rowInfo.numEntries;
      const size_t newNumEntries =
        curNumEntries + static_cast<size_t> (numEntriesToInsert);
      if (newNumEntries > rowInfo.allocSize) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          getProfileType () == StaticProfile && newNumEntries > rowInfo.allocSize,
          std::runtime_error, "Tpetra::CrsMatrix::insertGlobalValues: new "
          "indices exceed statically allocated graph structure.  curNumEntries"
          " (" << curNumEntries << ") + numEntriesToInsert (" <<
          numEntriesToInsert << ") > allocSize (" << rowInfo.allocSize << ").");

        // Update allocation only as much as necessary
        try {
          rowInfo =
            myGraph_->template updateGlobalAllocAndValues<impl_scalar_type> (rowInfo,
                                                                        newNumEntries,
                                                                        values2D_[localRow]);
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::runtime_error, "myGraph_->updateGlobalAllocAndValues"
            "(...) threw an exception: " << e.what ());
        }
      }
      try {
        if (isGloballyIndexed ()) {
          // lg=GlobalIndices, I=GlobalIndices means the method calls
          // getGlobalViewNonConst() and does direct copying, which
          // should be reasonably fast.
          myGraph_->template insertIndicesAndValues<impl_scalar_type> (rowInfo, inds_view,
                                                                  this->getViewNonConst (rowInfo),
                                                                  vals_view,
                                                                  GlobalIndices, GlobalIndices);
        }
        else {
          // lg=GlobalIndices, I=LocalIndices means the method calls
          // the Map's getLocalElement() method once per entry to
          // insert.  This may be slow.
          myGraph_->template insertIndicesAndValues<impl_scalar_type> (rowInfo, inds_view,
                                                                  this->getViewNonConst (rowInfo),
                                                                  vals_view,
                                                                  GlobalIndices, LocalIndices);
        }
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "myGraph_->insertIndicesAndValues(...) "
          "threw an exception: " << e.what ());
      }

#ifdef HAVE_TPETRA_DEBUG
      const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow (localRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries,
        std::logic_error, ": There should be a total of " << newNumEntries
        << " entries in the row, but the graph now reports " << chkNewNumEntries
        << " entries.  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  insertGlobalValuesFiltered (const GlobalOrdinal globalRow,
                              const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                              const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;
    const char tfecfFuncName[] = "insertGlobalValuesFiltered: ";

    // mfh 14 Dec 2012: Defer test for static graph until we know that
    // globalRow is in the row Map.  If it's not in the row Map, it
    // doesn't matter whether or not the graph is static; the data
    // just get stashed for later use by globalAssemble().
    //
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   isStaticGraph(), std::runtime_error,
    //   ": matrix was constructed with static graph. Cannot insert new entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      values.size () != indices.size (), std::runtime_error,
      "values.size() = " << values.size() << " != indices.size() = "
      << indices.size() << ".");
#endif // HAVE_TPETRA_DEBUG

    ArrayView<const ST> valsIn = av_reinterpret_cast<const ST> (values);
    const LO lrow = getRowMap ()->getLocalElement (globalRow);

    if (lrow != Teuchos::OrdinalTraits<LO>::invalid ()) { // globalRow is in our row Map.
      // If the matrix has a static graph, this process is now allowed
      // to insert into rows it owns.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        this->isStaticGraph (), std::runtime_error,
        "The matrix was constructed with a static graph.  In that case, "
        "it is forbidden to insert new entries into rows owned by the "
        "calling process.");
      if (! myGraph_->indicesAreAllocated ()) {
        allocateValues (GlobalIndices, GraphNotYetAllocated);
      }
      typename Graph::SLocalGlobalViews inds_view;
      ArrayView<const ST> vals_view;

      // We have to declare these Arrays here rather than in the
      // hasColMap() if branch, so that views to them will remain
      // valid for the whole scope.
      Array<GO> filtered_indices;
      Array<ST> filtered_values;
      if (hasColMap ()) { // We have a column Map.
        // Use column Map to filter the indices and corresponding
        // values, so that we only insert entries into columns we own.
        filtered_indices.assign (indices.begin (), indices.end ());
        filtered_values.assign (valsIn.begin (), valsIn.end ());
        const size_t numFilteredEntries =
          myGraph_->template filterGlobalIndicesAndValues<ST> (filtered_indices (),
                                                               filtered_values ());
        inds_view.ginds = filtered_indices (0, numFilteredEntries);
        vals_view       = filtered_values (0, numFilteredEntries);
      }
      else { // we don't have a column Map.
        inds_view.ginds = indices;
        vals_view       = valsIn;
      }
      const size_t numFilteredEntries = vals_view.size ();
      // add the new indices and values
      if (numFilteredEntries > 0) {
        RowInfo rowInfo = myGraph_->getRowInfo (lrow);
        const size_t curNumEntries = rowInfo.numEntries;
        const size_t newNumEntries = curNumEntries + numFilteredEntries;
        if (newNumEntries > rowInfo.allocSize) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            getProfileType () == StaticProfile, std::runtime_error,
            "New indices exceed statically allocated graph structure.");

          // Update allocation only as much as necessary
          rowInfo = myGraph_->template updateGlobalAllocAndValues<ST> (rowInfo,
                                                                       newNumEntries,
                                                                       values2D_[lrow]);
        }
        if (isGloballyIndexed ()) {
          // lg=GlobalIndices, I=GlobalIndices means the method calls
          // getGlobalViewNonConst() and does direct copying, which
          // should be reasonably fast.
          myGraph_->template insertIndicesAndValues<ST> (rowInfo, inds_view,
                                                         this->getViewNonConst (rowInfo),
                                                         vals_view,
                                                         GlobalIndices, GlobalIndices);
        }
        else {
          // lg=GlobalIndices, I=LocalIndices means the method calls
          // the Map's getLocalElement() method once per entry to
          // insert.  This may be slow.
          myGraph_->template insertIndicesAndValues<ST> (rowInfo, inds_view,
                                                         this->getViewNonConst (rowInfo),
                                                         vals_view,
                                                         GlobalIndices, LocalIndices);
        }
#ifdef HAVE_TPETRA_DEBUG
        {
          const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow (lrow);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries,
            std::logic_error, ": There should be a total of " << newNumEntries
            << " entries in the row, but the graph now reports " << chkNewNumEntries
            << " entries.  Please report this bug to the Tpetra developers.");
        }
#endif // HAVE_TPETRA_DEBUG
      }
    }
    else { // The calling process doesn't own the given row.
      insertNonownedGlobalValues (globalRow, indices, values);
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  LocalOrdinal
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceLocalValues (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal> &indices,
                      const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;
    // project2nd is a binary function that returns its second
    // argument.  This replaces entries in the given row with their
    // corresponding entry of values.
    typedef Tpetra::project2nd<ST, ST> f_type;
    typedef typename ArrayView<GO>::size_type size_type;

    ArrayView<const ST> valsIn = av_reinterpret_cast<const ST> (values);
    if (! isFillActive ()) {
      // Fill must be active in order to call this method.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (! this->hasColMap ()) {
      // There is no such thing as local column indices without a column Map.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (values.size () != indices.size ()) {
      // The sizes of values and indices must match.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const bool isLocalRow = getRowMap ()->isNodeLocalElement (localRow);
    if (! isLocalRow) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      //
      // FIXME (mfh 02 Jan 2015) replaceGlobalValues returns invalid
      // in this case.
      return static_cast<LO> (0);
    }

    if (indices.size () == 0) {
      return static_cast<LO> (0);
    }
    else {
      RowInfo rowInfo = staticGraph_->getRowInfo (localRow);
      ArrayView<ST> curVals = this->getViewNonConst (rowInfo);
      if (isLocallyIndexed ()) {
        return staticGraph_->template transformLocalValues<ST, f_type> (rowInfo,
                                                                        curVals,
                                                                        indices,
                                                                        valsIn,
                                                                        f_type ());
      }
      else if (isGloballyIndexed ()) {
        // Convert the given local indices to global indices.
        //
        // FIXME (mfh 27 Jun 2014) Why can't we ask the graph to do
        // that?  It could do the conversions in place, so that we
        // wouldn't need temporary storage.
        const map_type& colMap = * (this->getColMap ());
        const size_type numInds = indices.size ();

        // mfh 27 Jun 2014: Some of the given local indices might be
        // invalid.  That's OK, though, since the graph ignores them
        // and their corresponding values in transformGlobalValues.
        // Thus, we don't have to count how many indices are valid.
        // We do so just as a sanity check.
        Array<GO> gblInds (numInds);
        size_type numValid = 0; // sanity check count of # valid indices
        for (size_type k = 0; k < numInds; ++k) {
          const GO gid = colMap.getGlobalElement (indices[k]);
          gblInds[k] = gid;
          if (gid != Teuchos::OrdinalTraits<GO>::invalid ()) {
            ++numValid; // sanity check count of # valid indices
          }
        }
        const LO numXformed =
          staticGraph_->template transformGlobalValues<ST, f_type> (rowInfo,
                                                                    curVals, // target
                                                                    gblInds,
                                                                    valsIn, // source
                                                                    f_type ());
        if (static_cast<size_type> (numXformed) != numValid) {
          return Teuchos::OrdinalTraits<LO>::invalid ();
        } else {
          return numXformed;
        }
      }
      // NOTE (mfh 26 Jun 2014) In the current version of CrsMatrix,
      // it's possible for a matrix (or graph) to be neither locally
      // nor globally indexed on a process.  This means that the graph
      // or matrix has no entries on that process.  Epetra also works
      // like this.  It's related to lazy allocation (on first
      // insertion, not at graph / matrix construction).  Lazy
      // allocation will go away because it is not thread scalable.
      return static_cast<LO> (0);
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  LocalOrdinal
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceGlobalValues (GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                       const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;
    // project2nd is a binary function that returns its second
    // argument.  This replaces entries in the given row with their
    // corresponding entry of values.
    typedef Tpetra::project2nd<ST, ST> f_type;
    typedef typename ArrayView<GO>::size_type size_type;

    ArrayView<const ST> valsIn = av_reinterpret_cast<const ST> (values);
    if (! isFillActive ()) {
      // Fill must be active in order to call this method.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (values.size () != indices.size ()) {
      // The sizes of values and indices must match.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }

    const LO lrow = this->getRowMap ()->getLocalElement (globalRow);
    if (lrow == Teuchos::OrdinalTraits<LO>::invalid ()) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      //
      // FIXME (mfh 02 Jan 2015) replaceLocalValues returns 0 in this case.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }

    if (staticGraph_.is_null ()) {
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = *staticGraph_;
    RowInfo rowInfo = graph.getRowInfo (lrow);
    if (indices.size () == 0) {
      return static_cast<LO> (0);
    }
    else {
      ArrayView<ST> curVals = this->getViewNonConst (rowInfo);
      if (isLocallyIndexed ()) {
        // Convert the given global indices to local indices.
        //
        // FIXME (mfh 08 Jul 2014) Why can't we ask the graph to do
        // that?  It could do the conversions in place, so that we
        // wouldn't need temporary storage.
        const map_type& colMap = * (this->getColMap ());
        const size_type numInds = indices.size ();
        Array<LO> lclInds (numInds);
        for (size_type k = 0; k < numInds; ++k) {
          // There is no need to filter out indices not in the
          // column Map.  Those that aren't will be mapped to
          // invalid(), which the graph's transformGlobalValues()
          // will filter out (but not count in its return value).
          lclInds[k] = colMap.getLocalElement (indices[k]);
        }
        return graph.template transformLocalValues<ST, f_type> (rowInfo,
                                                                curVals,
                                                                lclInds (),
                                                                valsIn,
                                                                f_type ());
      }
      else if (isGloballyIndexed ()) {
        return graph.template transformGlobalValues<ST, f_type> (rowInfo,
                                                                 curVals,
                                                                 indices,
                                                                 valsIn,
                                                                 f_type ());
      }
      else {
        // If the graph is neither locally nor globally indexed on
        // the calling process, that means that the calling process
        // can't possibly have any entries in the owned row.  Thus,
        // there are no entries to transform, so we return zero.
        return static_cast<LO> (0);
      }
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  LocalOrdinal
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoGlobalValues (const GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                       const Teuchos::ArrayView<const Scalar>& values)

  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;
    typedef std::plus<Scalar> f_type;
    typedef typename ArrayView<GO>::size_type size_type;

    ArrayView<const ST> valsIn = av_reinterpret_cast<const ST> (values);
    if (! isFillActive ()) {
      // Fill must be active in order to call this method.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (values.size () != indices.size ()) {
      // The sizes of values and indices must match.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }

    const LO lrow = this->getRowMap ()->getLocalElement (globalRow);
    if (lrow == Teuchos::OrdinalTraits<LO>::invalid ()) {
      // globalRow is not in the row Map, so stash the given entries
      // away in a separate data structure.  globalAssemble() (called
      // during fillComplete()) will exchange that data and sum it in
      // using sumIntoGlobalValues().
      this->insertNonownedGlobalValues (globalRow, indices, values);
      // FIXME (mfh 08 Jul 2014) It's not clear what to return here,
      // since we won't know whether the given indices were valid
      // until globalAssemble (called in fillComplete) is called.
      // That's why insertNonownedGlobalValues doesn't return
      // anything.  Just for consistency, I'll return the number of
      // entries that the user gave us.
      return static_cast<LO> (indices.size ());
    }

    if (staticGraph_.is_null ()) {
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = *staticGraph_;
    RowInfo rowInfo = graph.getRowInfo (lrow);
    if (indices.size () == 0) {
      return static_cast<LO> (0);
    }
    else {
      ArrayView<ST> curVals = this->getViewNonConst (rowInfo);
      if (isLocallyIndexed ()) {
        // Convert the given global indices to local indices.
        //
        // FIXME (mfh 08 Jul 2014) Why can't we ask the graph to do
        // that?  It could do the conversions in place, so that we
        // wouldn't need temporary storage.
        const map_type& colMap = * (this->getColMap ());
        const size_type numInds = indices.size ();
        Array<LO> lclInds (numInds);
        for (size_type k = 0; k < numInds; ++k) {
          // There is no need to filter out indices not in the
          // column Map.  Those that aren't will be mapped to
          // invalid(), which the graph's transformGlobalValues()
          // will filter out (but not count in its return value).
          lclInds[k] = colMap.getLocalElement (indices[k]);
        }
        return graph.template transformLocalValues<ST, f_type> (rowInfo,
                                                                curVals,
                                                                lclInds (),
                                                                valsIn,
                                                                f_type ());
      }
      else if (isGloballyIndexed ()) {
        return graph.template transformGlobalValues<ST, f_type> (rowInfo,
                                                                 curVals,
                                                                 indices,
                                                                 valsIn,
                                                                 f_type ());
      }
      else {
        // If the graph is neither locally nor globally indexed on
        // the calling process, that means that the calling process
        // can't possibly have any entries in the owned row.  Thus,
        // there are no entries to transform, so we return zero.
        return static_cast<LO> (0);
      }
    }
  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  LocalOrdinal
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoLocalValues (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal>& indices,
                      const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;
    typedef std::plus<Scalar> f_type;
    typedef typename ArrayView<GO>::size_type size_type;

    ArrayView<const ST> valsIn = av_reinterpret_cast<const ST> (values);
    if (! isFillActive ()) {
      // Fill must be active in order to call this method.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (! this->hasColMap ()) {
      // There is no such thing as local column indices without a column Map.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else if (values.size () != indices.size ()) {
      // The sizes of values and indices must match.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const bool isLocalRow = getRowMap ()->isNodeLocalElement (localRow);
    if (! isLocalRow) {
      // The calling process doesn't own the local row, so we can't
      // insert into it.
      return static_cast<LO> (0);
    }

    if (indices.size () == 0) {
      return static_cast<LO> (0);
    }
    else {
      RowInfo rowInfo = staticGraph_->getRowInfo (localRow);
      ArrayView<ST> curVals = this->getViewNonConst (rowInfo);
      if (isLocallyIndexed ()) {
        return staticGraph_->template transformLocalValues<ST, f_type> (rowInfo,
                                                                        curVals,
                                                                        indices,
                                                                        valsIn,
                                                                        f_type ());
      }
      else if (isGloballyIndexed ()) {
        // Convert the given local indices to global indices.
        //
        // FIXME (mfh 27 Jun 2014) Why can't we ask the graph to do
        // that?  It could do the conversions in place, so that we
        // wouldn't need temporary storage.
        const map_type& colMap = * (this->getColMap ());
        const size_type numInds = indices.size ();

        // mfh 27 Jun 2014: Some of the given local indices might be
        // invalid.  That's OK, though, since the graph ignores them
        // and their corresponding values in transformGlobalValues.
        // Thus, we don't have to count how many indices are valid.
        // We do so just as a sanity check.
        Array<GO> gblInds (numInds);
        size_type numValid = 0; // sanity check count of # valid indices
        for (size_type k = 0; k < numInds; ++k) {
          const GO gid = colMap.getGlobalElement (indices[k]);
          gblInds[k] = gid;
          if (gid != Teuchos::OrdinalTraits<GO>::invalid ()) {
            ++numValid; // sanity check count of # valid indices
          }
        }
        const LO numXformed =
          staticGraph_->template transformGlobalValues<ST, f_type> (rowInfo,
                                                                    curVals, // target
                                                                    gblInds,
                                                                    valsIn, // source
                                                                    f_type ());
        if (static_cast<size_type> (numXformed) != numValid) {
          return Teuchos::OrdinalTraits<LO>::invalid ();
        } else {
          return numXformed;
        }
      }
      // NOTE (mfh 26 Jun 2014) In the current version of CrsMatrix,
      // it's possible for a matrix (or graph) to be neither locally
      // nor globally indexed on a process.  This means that the graph
      // or matrix has no entries on that process.  Epetra also works
      // like this.  It's related to lazy allocation (on first
      // insertion, not at graph / matrix construction).  Lazy
      // allocation will go away because it is not thread scalable.
      return static_cast<LO> (0);
    }
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  Teuchos::ArrayView<const typename CrsMatrix<
                       Scalar, LocalOrdinal, GlobalOrdinal,
                       Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::impl_scalar_type>
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getView (RowInfo rowinfo) const
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;
    using Teuchos::ArrayView;
    typedef impl_scalar_type ST;
    typedef std::pair<size_t, size_t> range_type;

    if (k_values1D_.dimension_0 () != 0 && rowinfo.allocSize > 0) {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        rowinfo.offset1D + rowinfo.allocSize > k_values1D_.dimension_0 (),
        std::range_error, "Tpetra::CrsMatrix::getView: Invalid access "
        "to 1-D storage of values." << std::endl << "rowinfo.offset1D (" <<
        rowinfo.offset1D << ") + rowinfo.allocSize (" << rowinfo.allocSize <<
        ") > k_values1D_.dimension_0() (" << k_values1D_.dimension_0 () << ").");
#endif // HAVE_TPETRA_DEBUG
      range_type range (rowinfo.offset1D, rowinfo.offset1D + rowinfo.allocSize);
      typedef View<const ST*, execution_space, MemoryUnmanaged> subview_type;
      subview_type sv = Kokkos::subview (k_values1D_, range);

      const ST* const sv_raw = (rowinfo.allocSize == 0) ? NULL : sv.ptr_on_device ();
      return ArrayView<const ST> (sv_raw, rowinfo.allocSize);
    }
    else if (values2D_ != null) {
      return values2D_[rowinfo.localRow] ();
    }
    else {
      return ArrayView<impl_scalar_type> ();
    }
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  Teuchos::ArrayView<typename CrsMatrix<
                       Scalar, LocalOrdinal, GlobalOrdinal,
                       Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::impl_scalar_type>
  CrsMatrix<
    Scalar, LocalOrdinal,GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getViewNonConst (RowInfo rowinfo)
  {
    return Teuchos::av_const_cast<impl_scalar_type> (this->getView (rowinfo));
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalRowCopy (LocalOrdinal localRow,
                   const Teuchos::ArrayView<LocalOrdinal>& indices,
                   const Teuchos::ArrayView<Scalar>& values,
                   size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed () && ! hasColMap (), std::runtime_error,
      "Tpetra::CrsMatrix::getLocalRowCopy: The matrix is globally indexed and "
      "does not have a column Map yet.  That means we don't have local indices "
      "for columns yet, so it doesn't make sense to call this method.  If the "
      "matrix doesn't have a column Map yet, you should call fillComplete on "
      "it first.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! staticGraph_->hasRowInfo (), std::runtime_error,
      "Tpetra::CrsMatrix::getLocalRowCopy: The graph's row information was "
      "deleted at fillComplete().");

    if (! this->getRowMap ()->isNodeLocalElement (localRow)) {
      // The calling process owns no entries in this row.
      numEntries = 0;
      return;
    }

    const RowInfo rowinfo = staticGraph_->getRowInfo (localRow);
    const size_t theNumEntries = rowinfo.numEntries;

    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (indices.size ()) < theNumEntries ||
      static_cast<size_t> (values.size ()) < theNumEntries,
      std::runtime_error,
      "Tpetra::CrsMatrix::getLocalRowCopy: The given row " << localRow
      << " has " << theNumEntries << " entries.  One or both of the given "
      "ArrayViews are not long enough to store that many entries.  indices "
      "can store " << indices.size() << " entries and values can store "
      << values.size() << " entries.");

    numEntries = theNumEntries;

    if (staticGraph_->isLocallyIndexed ()) {
      ArrayView<const LO> indrowview = staticGraph_->getLocalView (rowinfo);
      ArrayView<const Scalar> valrowview =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
      std::copy (indrowview.begin (), indrowview.begin () + numEntries, indices.begin ());
      std::copy (valrowview.begin (), valrowview.begin () + numEntries, values.begin ());
    }
    else if (staticGraph_->isGloballyIndexed ()) {
      ArrayView<const GO> indrowview = staticGraph_->getGlobalView (rowinfo);
      ArrayView<const Scalar> valrowview =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
      std::copy (valrowview.begin (), valrowview.begin () + numEntries, values.begin ());

      const map_type& colMap = * (this->getColMap ());
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = colMap.getLocalElement (indrowview[j]);
      }
    }
    else {
      numEntries = 0;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalRowCopy (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<GlobalOrdinal>& indices,
                    const Teuchos::ArrayView<Scalar>& values,
                    size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const char tfecfFuncName[] = "getGlobalRowCopy: ";
    const LocalOrdinal lrow = getRowMap ()->getLocalElement (globalRow);
    if (lrow == OTL::invalid ()) {
      // The calling process owns no entries in this row.
      numEntries = 0;
      return;
    }

    const RowInfo rowinfo = staticGraph_->getRowInfo (lrow);
    const size_t theNumEntries = rowinfo.numEntries;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < theNumEntries ||
      static_cast<size_t> (values.size ()) < theNumEntries,
      std::runtime_error,
      "The given row " << globalRow << ", corresponding to local row " << lrow
      << ", has " << theNumEntries << " entries.  One or both of the given "
      "ArrayView input arguments are not long enough to store that many "
      "entries.  indices.size() = " << indices.size() << ", values.size() = "
      << values.size () << ", but the number of entries in the row is "
      << theNumEntries << ".");

    // Don't "commit" the value until we know that the input arrays are valid.
    numEntries = theNumEntries;

    if (staticGraph_->isGloballyIndexed ()) {
      ArrayView<const GO> indrowview = staticGraph_->getGlobalView (rowinfo);
      ArrayView<const Scalar> valrowview =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
      std::copy (indrowview.begin (), indrowview.begin () + numEntries, indices.begin ());
      std::copy (valrowview.begin (), valrowview.begin () + numEntries, values.begin ());
    }
    else if (staticGraph_->isLocallyIndexed ()) {
      ArrayView<const LO> indrowview = staticGraph_->getLocalView(rowinfo);
      ArrayView<const Scalar> valrowview =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
      std::copy (valrowview.begin (), valrowview.begin () + numEntries, values.begin ());
      for (size_t j = 0; j < numEntries; ++j) {
        indices[j] = getColMap ()->getGlobalElement (indrowview[j]);
      }
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above if indices are allocated
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        staticGraph_->indicesAreAllocated (), std::logic_error,
        "Internal logic error. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
      numEntries = 0;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalRowView (LocalOrdinal localRow,
                   Teuchos::ArrayView<const LocalOrdinal>& indices,
                   Teuchos::ArrayView<const Scalar>& values) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;

    const char tfecfFuncName[] = "getLocalRowView: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "The matrix currently stores "
      "its indices as global indices, so you cannot get a view with local "
      "column indices.  If the matrix has a column Map, you may call "
      "getLocalRowCopy() to get local column indices; otherwise, you may get "
      "a view with global column indices by calling getGlobalRowCopy().");
    indices = Teuchos::null;
    values = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
    size_t numEntries = 0;
#endif // HAVE_TPETRA_DEBUG
    if (getRowMap ()->isNodeLocalElement (localRow)) {
      const RowInfo rowinfo = staticGraph_->getRowInfo (localRow);
#ifdef HAVE_TPETRA_DEBUG
      numEntries = rowinfo.numEntries;
#endif // HAVE_TPETRA_DEBUG
      if (rowinfo.numEntries > 0) {
        ArrayView<const LO> indTmp = staticGraph_->getLocalView (rowinfo);
        ArrayView<const Scalar> valTmp =
          av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
        indices = indTmp (0, rowinfo.numEntries);
        values = valTmp (0, rowinfo.numEntries);
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    const char suffix[] = ".  This should never happen.  Please report this "
      "bug to the Tpetra developers.";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        static_cast<size_t>(indices.size ()) != static_cast<size_t>(values.size ()), std::logic_error,
      "At the end of this method, for local row " << localRow << ", "
      "indices.size() = " << indices.size () << " != values.size () = "
      << values.size () << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        static_cast<size_t>(indices.size ()) != static_cast<size_t>(numEntries), std::logic_error,
      "At the end of this method, for local row " << localRow << ", "
      "indices.size() = " << indices.size () << " != numEntries = "
      << numEntries << suffix);
    const size_t expectedNumEntries = this->getNumEntriesInLocalRow (localRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numEntries != expectedNumEntries, std::logic_error,
      "At the end of this method, for local row " << localRow << ", numEntries"
      " = " << numEntries << " != getNumEntriesInLocalRow(localRow)"
      " = "<< expectedNumEntries << suffix);
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalRowView (GlobalOrdinal globalRow,
                    Teuchos::ArrayView<const GlobalOrdinal>& indices,
                    Teuchos::ArrayView<const Scalar>& values) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "getGlobalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      "The matrix is locally indexed, so we cannot return a view of the row "
      "with global column indices.  Use getGlobalRowCopy() instead.");
    indices = Teuchos::null;
    values  = Teuchos::null;
    const LocalOrdinal lrow = getRowMap ()->getLocalElement (globalRow);
    if (lrow != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      // getRowInfo() requires a local row index, whether or not
      // storage has been optimized.
      const RowInfo rowinfo = staticGraph_->getRowInfo(lrow);
      if (rowinfo.numEntries > 0) {
        ArrayView<const GO> indTmp = staticGraph_->getGlobalView (rowinfo);
        ArrayView<const Scalar> valTmp =
          av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
        indices = indTmp (0, rowinfo.numEntries);
        values = valTmp (0, rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) != this->getNumEntriesInGlobalRow (globalRow) ||
      indices.size () != values.size (),
      std::logic_error,
      "Violated stated post-conditions. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  scale (const Scalar& alpha)
  {
    typedef LocalOrdinal LO;
    typedef Kokkos::SparseRowView<local_matrix_type> row_view_type;
    typedef typename Teuchos::Array<Scalar>::size_type size_type;
    const char tfecfFuncName[] = "scale: ";
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error,
      "Fill must be active before you may call this method.  "
      "Please call resumeFill() to make fill active.");

    const size_t nlrs = staticGraph_->getNodeNumRows ();
    const size_t numAlloc = staticGraph_->getNodeAllocationSize ();
    const size_t numEntries = staticGraph_->getNodeNumEntries ();
    if (! staticGraph_->indicesAreAllocated () || nlrs == 0 ||
        numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      if (staticGraph_->getProfileType () == StaticProfile) {
        const LO lclNumRows = lclMatrix_.numRows ();
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          row_view_type row_i = lclMatrix_.template row<typename row_view_type::size_type> (lclRow);
          for (LO k = 0; k < row_i.length; ++k) {
            // FIXME (mfh 02 Jan 2015) This assumes CUDA UVM.
            row_i.value (k) *= theAlpha;
          }
        }
      }
      else if (staticGraph_->getProfileType () == DynamicProfile) {
        for (size_t row = 0; row < nlrs; ++row) {
          const size_type numEnt = getNumEntriesInLocalRow (row);
          Teuchos::ArrayView<impl_scalar_type> rowVals = values2D_[row] ();
          for (size_type k = 0; k < numEnt; ++k) {
            rowVals[k] *= theAlpha;
          }
        }
      }
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  setAllToScalar (const Scalar& alpha)
  {
    const char tfecfFuncName[] = "setAllToScalar: ";
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error,
      "Fill must be active before you may call this method.  "
      "Please call resumeFill() to make fill active.");

    // replace all values in the matrix
    // it is easiest to replace all allocated values, instead of replacing only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (no entry have been inserted so far)
    const size_t     nlrs = staticGraph_->getNodeNumRows(),
                 numAlloc = staticGraph_->getNodeAllocationSize(),
               numEntries = staticGraph_->getNodeNumEntries();
    if (! staticGraph_->indicesAreAllocated () || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      const ProfileType profType = staticGraph_->getProfileType ();
      if (profType == StaticProfile) {
        // FIXME (mfh 24 Dec 2014) Once CrsMatrix implements DualView
        // semantics, this would be the place to mark memory as
        // modified.
        typedef typename local_matrix_type::values_type values_type;
        Kokkos::Impl::ViewFill<values_type> (k_values1D_, theAlpha);
      }
      else if (profType == DynamicProfile) {
        for (size_t row = 0; row < nlrs; ++row) {
          std::fill (values2D_[row].begin (), values2D_[row].end (), theAlpha);
        }
      }
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  setAllValues (const typename local_matrix_type::row_map_type& rowPointers,
                const typename local_graph_type::entries_type::non_const_type& columnIndices,
                const typename local_matrix_type::values_type& values)
  {
    const char tfecfFuncName[] = "setAllValues";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      columnIndices.size () != values.size (), std::runtime_error,
      ": columnIndices and values must have the same size.  columnIndices.size() = "
      << columnIndices.size () << " != values.size() = " << values.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error, ": myGraph_ must not be null.");

    try {
      myGraph_->setAllIndices (rowPointers, columnIndices);
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error, ": Caught exception while calling myGraph_->"
        "setAllIndices(): " << e.what ());
    }
    k_values1D_ = values;
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  setAllValues (const Teuchos::ArrayRCP<size_t>& rowPointers,
                const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
                const Teuchos::ArrayRCP<Scalar>& values)
  {
    const char tfecfFuncName[] = "setAllValues: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      columnIndices.size () != values.size (), std::runtime_error,
      "columnIndices.size() = " << columnIndices.size () << " != "
      "values.size() = " << values.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error, "myGraph_ must not be null.");

    try {
      myGraph_->setAllIndices (rowPointers, columnIndices);
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error, "Caught exception while calling myGraph_->"
        "setAllIndices(): " << e.what ());
    }
    Teuchos::ArrayRCP<impl_scalar_type> vals =
      Teuchos::arcp_reinterpret_cast<impl_scalar_type> (values);
    k_values1D_ = Kokkos::Compat::getKokkosViewDeepCopy<DeviceType> (vals ());
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "getLocalDiagOffsets";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error,
      ": This method requires that the matrix have a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      staticGraph_.is_null (), std::runtime_error,
      ": This method requires that the matrix have a graph.");

    const map_type& rowMap = * (this->getRowMap ());
    const map_type& colMap = * (this->getColMap ());

    const size_t myNumRows = getNodeNumRows ();
    if (static_cast<size_t> (offsets.size ()) != myNumRows) {
      offsets.resize (static_cast<size_t> (myNumRows));
    }

#ifdef HAVE_TPETRA_DEBUG
    bool allRowMapDiagEntriesInColMap = true;
    bool allDiagEntriesFound = true;
#endif // HAVE_TPETRA_DEBUG

    for (size_t r = 0; r < myNumRows; ++r) {
      const GlobalOrdinal rgid = rowMap.getGlobalElement (r);
      const LocalOrdinal rlid = colMap.getLocalElement (rgid);

#ifdef HAVE_TPETRA_DEBUG
      if (rlid == Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
        allRowMapDiagEntriesInColMap = false;
      }
#endif // HAVE_TPETRA_DEBUG

      if (rlid != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
        RowInfo rowinfo = staticGraph_->getRowInfo (r);
        if (rowinfo.numEntries > 0) {
          offsets[r] = staticGraph_->findLocalIndex (rowinfo, rlid);
        }
        else {
          offsets[r] = Teuchos::OrdinalTraits<size_t>::invalid ();
#ifdef HAVE_TPETRA_DEBUG
          allDiagEntriesFound = false;
#endif // HAVE_TPETRA_DEBUG
        }
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    using Teuchos::reduceAll;
    using std::endl;

    const bool localSuccess =
      allRowMapDiagEntriesInColMap && allDiagEntriesFound;
    int localResults[3];
    localResults[0] = allRowMapDiagEntriesInColMap ? 1 : 0;
    localResults[1] = allDiagEntriesFound ? 1 : 0;
    // min-all-reduce will compute least rank of all the processes
    // that didn't succeed.
    localResults[2] =
      ! localSuccess ? getComm ()->getRank () : getComm ()->getSize ();
    int globalResults[3];
    globalResults[0] = 0;
    globalResults[1] = 0;
    globalResults[2] = 0;
    reduceAll<int, int> (* (getComm ()), Teuchos::REDUCE_MIN,
                         3, localResults, globalResults);
    if (globalResults[0] == 0 || globalResults[1] == 0) {
      std::ostringstream os; // build error message
      const bool both =
        globalResults[0] == 0 && globalResults[1] == 0;
      os << ": At least one process (including Process " << globalResults[2]
         << ") had the following issue" << (both ? "s" : "") << ":" << endl;
      if (globalResults[0] == 0) {
        os << "  - The column Map does not contain at least one diagonal entry "
          "of the matrix." << endl;
      }
      if (globalResults[1] == 0) {
        os << "  - There is a row on that / those process(es) that does not "
          "contain a diagonal entry." << endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error, os.str());
    }
#endif // HAVE_TPETRA_DEBUG
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& dvec) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    const char tfecfFuncName[] = "getLocalDiagCopy: ";
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> vec_type;
    typedef typename vec_type::dual_view_type dual_view_type;
    typedef typename dual_view_type::host_mirror_space::execution_space host_execution_space;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error,
      "This method requires that the matrix have a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      staticGraph_.is_null (), std::runtime_error,
      "This method requires that the matrix have a graph.");
    const map_type& rowMap = * (this->getRowMap ());
    const map_type& colMap = * (this->getColMap ());

#ifdef HAVE_TPETRA_DEBUG
    // isCompatible() requires an all-reduce, and thus this check
    // should only be done in debug mode.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! dvec.getMap ()->isCompatible (rowMap), std::runtime_error,
      ": The input Vector's Map must be compatible with the CrsMatrix's row "
      "Map.  You may check this by using Map's isCompatible method: "
      "dvec.getMap ()->isCompatible (A.getRowMap ());");
#endif // HAVE_TPETRA_DEBUG

    // For now, we fill the Vector on the host and sync to device.
    // Later, we may write a parallel kernel that works entirely on
    // device.
    dual_view_type lclVec = dvec.getDualView ();
    lclVec.template modify<host_execution_space> ();
    typedef typename dual_view_type::t_host host_view_type;
    host_view_type lclVecHost = lclVec.h_view;

    // 1-D subview of lclVecHost.  All the "typename" stuff ensures
    // that we get the same layout and memory traits as the original
    // 2-D view.
    typedef typename Kokkos::View<impl_scalar_type*,
      typename host_view_type::array_layout,
      typename host_view_type::device_type,
      typename host_view_type::memory_traits>
      host_view_1d_type;
    host_view_1d_type lclVecHost1d =
      Kokkos::subview (lclVecHost, Kokkos::ALL (), 0);

    // Find the diagonal entries and put them in lclVecHost1d.
    const size_t myNumRows = getNodeNumRows ();
    for (size_t r = 0; r < myNumRows; ++r) {
      lclVecHost1d(r) = STS::zero (); // default value if no diag entry
      const GlobalOrdinal rgid = rowMap.getGlobalElement (r);
      const LocalOrdinal rlid = colMap.getLocalElement (rgid);

      if (rlid != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
        RowInfo rowinfo = staticGraph_->getRowInfo (r);
        if (rowinfo.numEntries > 0) {
          const size_t j = staticGraph_->findLocalIndex (rowinfo, rlid);
          if (j != Teuchos::OrdinalTraits<size_t>::invalid ()) {
            // NOTE (mfh 02 Jan 2015) This technically does not assume
            // UVM, since getView and getViewNonConst are supposed to
            // return views of host data.
            ArrayView<const impl_scalar_type> view = this->getView (rowinfo);
            lclVecHost1d(r) = view[j];
          }
        }
      }
    }
    lclVec.template sync<execution_space> (); // sync changes back to device
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> vec_type;
    typedef typename vec_type::dual_view_type dual_view_type;
    typedef typename dual_view_type::host_mirror_space::execution_space host_execution_space;

#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getLocalDiagCopy: ";
    const map_type& rowMap = * (this->getRowMap ());
    // isCompatible() requires an all-reduce, and thus this check
    // should only be done in debug mode.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! diag.getMap ()->isCompatible (rowMap), std::runtime_error,
      "The input Vector's Map must be compatible with (in the sense of Map::"
      "isCompatible) the CrsMatrix's row Map.");
#endif // HAVE_TPETRA_DEBUG

    // For now, we fill the Vector on the host and sync to device.
    // Later, we may write a parallel kernel that works entirely on
    // device.
    dual_view_type lclVec = diag.getDualView ();
    lclVec.template modify<host_execution_space> ();
    typedef typename dual_view_type::t_host host_view_type;
    host_view_type lclVecHost = lclVec.h_view;

    // 1-D subview of lclVecHost.  All the "typename" stuff ensures
    // that we get the same layout and memory traits as the original
    // 2-D view.
    typedef typename Kokkos::View<impl_scalar_type*,
      typename host_view_type::array_layout,
      typename host_view_type::device_type,
      typename host_view_type::memory_traits>
      host_view_1d_type;
    host_view_1d_type lclVecHost1d =
      Kokkos::subview (lclVecHost, Kokkos::ALL (), 0);

    // Find the diagonal entries and put them in lclVecHost1d.
    const size_t myNumRows = getNodeNumRows ();
    for (size_t i = 0; i < myNumRows; ++i) {
      lclVecHost1d(i) = STS::zero (); // default value if no diag entry
      if (offsets[i] != Teuchos::OrdinalTraits<size_t>::invalid ()) {
        ArrayView<const LocalOrdinal> ind;
        ArrayView<const Scalar> val;
        // NOTE (mfh 02 Jan 2015) This technically does not assume
        // UVM, since the get{Global,Local}RowView methods are
        // supposed to return views of host data.
        this->getLocalRowView (i, ind, val);
        lclVecHost1d(i) = static_cast<impl_scalar_type> (val[offsets[i]]);
      }
    }
    lclVec.template sync<execution_space> (); // sync changes back to device
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& x)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> vec_type;
    const char tfecfFuncName[] = "leftScale";

    // FIXME (mfh 06 Aug 2014) This doesn't make sense.  The matrix
    // should only be modified when it is not fill complete.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete (), std::runtime_error,
      ": matrix must be fill complete.");
    RCP<const vec_type> xp;

    if (getRangeMap ()->isSameAs (* (x.getMap ()))){
      // Take from Epetra: If we have a non-trivial exporter, we must
      // import elements that are permuted or are on other processors.
      // (We will use the exporter to perform the import ("reverse
      // mode").)
      if (getCrsGraph ()->getExporter () != null) {
        RCP<vec_type> tempVec = rcp (new vec_type (getRowMap ()));
        tempVec->doImport (x, * (getCrsGraph ()->getExporter ()), INSERT);
        xp = tempVec;
      }
      else {
        xp = rcpFromRef (x);
      }
    }
    else if (getRowMap ()->isSameAs (* (x.getMap ()))) {
      xp = rcpFromRef (x);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::invalid_argument, ": The "
        "input scaling vector x's Map must be the same as either the row Map or "
        "the range Map of the CrsMatrix.");
    }
    ArrayRCP<const Scalar> vectorVals = xp->getData (0);
    ArrayView<impl_scalar_type> rowValues = null;

    const size_t lclNumRows = this->getNodeNumRows ();
    for (size_t i = 0; i < lclNumRows; ++i) {
      const RowInfo rowinfo = staticGraph_->getRowInfo (static_cast<LocalOrdinal> (i));
      rowValues = this->getViewNonConst (rowinfo);
      const impl_scalar_type scaleValue = static_cast<impl_scalar_type> (vectorVals[i]);
      for (size_t j = 0; j < rowinfo.numEntries; ++j) {
        rowValues[j] *= scaleValue;
      }
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& x)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> vec_type;
    const char tfecfFuncName[] = "rightScale: ";

    // FIXME (mfh 06 Aug 2014) This doesn't make sense.  The matrix
    // should only be modified when it is not fill complete.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete (), std::runtime_error, "Matrix must be fill complete.");
    RCP<const vec_type> xp;
    if (getDomainMap ()->isSameAs (* (x.getMap ()))) {
      // Take from Epetra: If we have a non-trivial exporter, we must
      // import elements that are permuted or are on other processors.
      // (We will use the exporter to perform the import.)
      if (getCrsGraph ()->getImporter () != null) {
        RCP<vec_type> tempVec = rcp (new vec_type (getColMap ()));
        tempVec->doImport (x, * (getCrsGraph ()->getImporter ()), INSERT);
        xp = tempVec;
      }
      else {
        xp = rcpFromRef (x);
      }
    }
    else if (getRowMap ()->isSameAs (* (x.getMap ()))) {
      xp = rcpFromRef (x);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error, "The vector x must have the same Map as "
        "either the row Map or the range Map.");
    }

    ArrayRCP<const Scalar> vectorVals = xp->getData (0);
    ArrayView<impl_scalar_type> rowValues = null;

    const size_t lclNumRows = this->getNodeNumRows ();
    for (size_t i = 0; i < lclNumRows; ++i) {
      const RowInfo rowinfo = staticGraph_->getRowInfo (static_cast<LocalOrdinal> (i));
      rowValues = this->getViewNonConst (rowinfo);
      ArrayView<const LocalOrdinal> colInds;
      getCrsGraph ()->getLocalRowView (static_cast<LocalOrdinal> (i), colInds);
      for (size_t j = 0; j < rowinfo.numEntries; ++j) {
        rowValues[j] *= static_cast<impl_scalar_type> (vectorVals[colInds[j]]);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::mag_type
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getFrobeniusNorm () const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    typedef typename Teuchos::ArrayRCP<const impl_scalar_type>::size_type size_type;

    // FIXME (mfh 05 Aug 2014) Write a thread-parallel kernel for the
    // local part of this computation.  It could make sense to put
    // this operation in the Kokkos::CrsMatrix.

    // check the cache first
    mag_type frobNorm = frobNorm_;
    if (frobNorm == -STM::one ()) {
      mag_type mySum = STM::zero ();
      if (getNodeNumEntries() > 0) {
        if (isStorageOptimized ()) {
          // "Optimized" storage is packed storage.  That means we can
          // iterate in one pass through the 1-D values array.
          const size_type numEntries =
            static_cast<size_type> (getNodeNumEntries ());
          for (size_type k = 0; k < numEntries; ++k) {
            // FIXME (mfh 05 Aug 2014) This assumes UVM.
            const impl_scalar_type val = k_values1D_(k);
            // Note (etp 06 Jan 2015) We need abs() here for composite types
            // (in general, if mag_type is on the left-hand-side, we need
            // abs() on the right-hand-side)
            const mag_type val_abs = STS::abs (val);
            mySum += val_abs * val_abs;
          }
        }
        else {
          const size_t numRows = getNodeNumRows ();
          for (size_t r = 0; r < numRows; ++r) {
            RowInfo rowInfo = myGraph_->getRowInfo (r);
            const size_type numEntries =
              static_cast<size_type> (rowInfo.numEntries);
            ArrayView<const impl_scalar_type> A_r =
              this->getView (rowInfo).view (0, numEntries);
            for (size_type k = 0; k < numEntries; ++k) {
              const impl_scalar_type val = A_r[k];
              const mag_type val_abs = STS::abs (val);
              mySum += val_abs * val_abs;
            }
          }
        }
      }
      mag_type totalSum = STM::zero ();
      reduceAll<int, mag_type> (* (getComm ()), REDUCE_SUM,
                                mySum, outArg (totalSum));
      frobNorm = STM::sqrt (totalSum);
    }
    if (isFillComplete ()) {
      // Only cache the result if the matrix is fill complete.
      // Otherwise, the values might still change.  resumeFill clears
      // the cache.
      frobNorm_ = frobNorm;
    }
    return frobNorm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceColMap (const Teuchos::RCP<const map_type>& newColMap)
  {
    const char tfecfFuncName[] = "replaceColMap: ";
    // FIXME (mfh 06 Aug 2014) What if the graph is locally indexed?
    // Then replacing the column Map might mean that we need to
    // reindex the column indices.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but "
      "this method necessarily must modify the graph, since the graph owns "
      "the matrix's column Map.");
    myGraph_->replaceColMap (newColMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reindexColumns (crs_graph_type* const graph,
                  const Teuchos::RCP<const map_type>& newColMap,
                  const Teuchos::RCP<const import_type>& newImport,
                  const bool sortEachRow)
  {
    const char tfecfFuncName[] = "reindexColumns: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      graph == NULL && myGraph_.is_null (), std::invalid_argument,
      "The input graph is NULL, but the matrix does not own its graph.");

    crs_graph_type& theGraph = (graph == NULL) ? *myGraph_ : *graph;
    const bool sortGraph = false; // we'll sort graph & matrix together below
    theGraph.reindexColumns (newColMap, newImport, sortGraph);
    if (sortEachRow && theGraph.isLocallyIndexed () && ! theGraph.isSorted ()) {
      // We can't just call sortEntries() here, because that fails if
      // the matrix has a const graph.  We want to use the given graph
      // in that case.
      const size_t lclNumRows = theGraph.getNodeNumRows ();
      for (size_t row = 0; row < lclNumRows; ++row) {
        RowInfo rowInfo = theGraph.getRowInfo (row);
        Teuchos::ArrayView<impl_scalar_type> rv = this->getViewNonConst (rowInfo);
        theGraph.template sortRowIndicesAndValues<impl_scalar_type> (rowInfo, rv);
      }
      theGraph.indicesAreSorted_ = true;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                               Teuchos::RCP<const import_type>& newImporter)
  {
    const char tfecfFuncName[] = "replaceDomainMapAndImporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but this"
      " method necessarily must modify the graph, since the graph owns the "
      "matrix's domain Map and Import objects.");
    myGraph_->replaceDomainMapAndImporter (newDomainMap, newImporter);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<Scalar,
            LocalOrdinal,
            GlobalOrdinal,
            Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  insertNonownedGlobalValues (const GlobalOrdinal globalRow,
                              const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                              const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    typedef GlobalOrdinal GO;
    typedef typename Array<GO>::size_type size_type;

    const size_type numToInsert = indices.size ();
    // Add the new data to the list of nonlocals.
    // This creates the arrays if they don't exist yet.
    std::pair<Array<GO>, Array<Scalar> >& curRow = nonlocals_[globalRow];
    Array<GO>& curRowInds = curRow.first;
    Array<Scalar>& curRowVals = curRow.second;
    const size_type newCapacity = curRowInds.size () + numToInsert;
    curRowInds.reserve (newCapacity);
    curRowVals.reserve (newCapacity);
    for (size_type k = 0; k < numToInsert; ++k) {
      curRowInds.push_back (indices[k]);
      curRowVals.push_back (values[k]);
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  globalAssemble ()
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::CommRequest;
    using Teuchos::gatherAll;
    using Teuchos::isend;
    using Teuchos::ireceive;
    using Teuchos::null;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::tuple;
    using Teuchos::waitAll;
    using std::make_pair;
    using std::pair;
    typedef GlobalOrdinal GO;
    typedef typename Array<GO>::size_type size_type;
    // nonlocals_ contains the entries stored by previous calls to
    // insertGlobalValues() for nonowned rows.
    typedef std::map<GO, pair<Array<GO>, Array<Scalar> > > nonlocals_map_type;
    typedef typename nonlocals_map_type::const_iterator nonlocals_iter_type;

    const char tfecfFuncName[] = "globalAssemble";
    const Teuchos::Comm<int>& comm = * (getComm ());
    const int numImages = comm.getSize ();
    const int myImageID = comm.getRank ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error, ": requires that fill is active.");

    // Determine (via a global all-reduce) if any processes have
    // nonlocal entries to share.  This is necessary even if the
    // matrix has a static graph, because insertGlobalValues allows
    // nonlocal entries in that case.
    size_t MyNonlocals = static_cast<size_t> (nonlocals_.size ());
    size_t MaxGlobalNonlocals = 0;
    reduceAll<int, size_t> (comm, REDUCE_MAX, MyNonlocals,
                            outArg (MaxGlobalNonlocals));
    if (MaxGlobalNonlocals == 0) {
      return;  // no entries to share
    }

    // FIXME (mfh 14 Dec 2012) The code below reimplements an Export
    // operation.  It would be better just to use an Export.  See
    // Comment #34 in discussion of Bug 5782.
    //
    // mfh 24 Feb 2014: On the other hand, this is not technically an
    // Export, since the row Map might not necessarily be one-to-one.

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images:
    //                  globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // Construct the set of all nonowned rows encountered by this
      // process in insertGlobalValues() or sumIntoGlobalValues().
      std::set<GlobalOrdinal> setOfRows;
      for (nonlocals_iter_type iter = nonlocals_.begin ();
           iter != nonlocals_.end (); ++iter) {
        setOfRows.insert (iter->first);
      }
      // Copy the resulting set of nonowned rows into an Array.
      Array<GlobalOrdinal> NLRs (setOfRows.size ());
      std::copy (setOfRows.begin (), setOfRows.end (), NLRs.begin ());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds (NLRs.size ());
      {
        const LookupStatus stat =
          getRowMap ()->getRemoteIndexList (NLRs (), NLRIds ());
        const int lclerr = (stat == IDNotPresent ? 1 : 0);
        int gblerr;
        reduceAll<int, int> (comm, REDUCE_MAX, lclerr, outArg (gblerr));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          gblerr, std::runtime_error, ": non-local entries correspond to "
          "invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve (NLRs.size ());
      Array<char> localNeighbors (numImages, 0);
      typename Array<GO>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin (), id = NLRIds.begin ();
           nlr != NLRs.end (); ++nlr, ++id) {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back (make_pair (*id, *nlr));
      }
      for (int j = 0; j < numImages; ++j) {
        if (localNeighbors[j]) {
          sendIDs.push_back (j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort (IdsAndRows.begin (), IdsAndRows.end ());
      // gather from other nodes to form the full graph
      //
      // FIXME (mfh 24 Feb 2014) Ugh, this is awful!!!  It's making a
      // P x P matrix which is the full graph of process connectivity.
      // Neither Export nor Import does this!  It would probably be
      // more efficient to do the following:
      //
      //   1. Form the one-to-one version of the row Map, tgtMap
      //   2. Form the (possibly overlapping) Map srcMap, with the
      //      global row indices which are the keys of nonlocals_ on
      //      each process
      //   3. Construct an Export from srcMap to tgtMap
      //   4. Execute the Export with Tpetra::ADD
      globalNeighbors.shapeUninitialized (numImages, numImages);
      gatherAll (comm, numImages, localNeighbors.getRawPtr (),
                 numImages*numImages, globalNeighbors.values ());
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
      if (globalNeighbors (myImageID, j)) {
        recvIDs.push_back (j);
      }
    }
    const size_t numRecvs = recvIDs.size ();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<Details::CrsIJV<GlobalOrdinal, Scalar> > IJVSendBuffer;
    Array<size_t> sendSizes (sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int, GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow)
    {
      const int id = IdAndRow->first;
      const GO row = IdAndRow->second;

      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          sendIDs[numSends] != id, std::logic_error,
          ": internal logic error. Contact Tpetra team.");
      }

      // copy data for row into contiguous storage
      pair<Array<GO>, Array<Scalar> >& nonlocalsRow = nonlocals_[row];
      ArrayView<const GO> nonlocalsRow_colInds = nonlocalsRow.first ();
      ArrayView<const Scalar> nonlocalsRow_values = nonlocalsRow.second ();
      const size_type numNonlocalsRow = nonlocalsRow_colInds.size ();

      for (size_type k = 0; k < numNonlocalsRow; ++k) {
        const Scalar val = nonlocalsRow_values[k];
        const GO col = nonlocalsRow_colInds[k];
        IJVSendBuffer.push_back (Details::CrsIJV<GO, Scalar> (row, col, val));
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size () > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_type> (numSends) != sendIDs.size(),
      std::logic_error, ": internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    // clear it before we start allocating a bunch of new memory
    nonlocals_.clear ();

    //////////////////////////////////////////////////////////////////////////////////////
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    //////////////////////////////////////////////////////////////////////////////////////
    // perform non-blocking sends: send sizes to our recipients
    Array<RCP<CommRequest<int> > > sendRequests;
    for (size_t s = 0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back (isend<int, size_t> (comm, rcpFromRef (sendSizes[s]), sendIDs[s]));
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<RCP<CommRequest<int> > > recvRequests;
    Array<size_t> recvSizes (numRecvs);
    for (size_t r = 0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication
      // will be local to this method and the scope of our data
      recvRequests.push_back (ireceive<int, size_t> (comm, rcpFromRef (recvSizes[r]), recvIDs[r]));
    }
    // wait on all
    if (! sendRequests.empty ()) {
      waitAll (comm, sendRequests ());
    }
    if (! recvRequests.empty ()) {
      waitAll (comm, recvRequests ());
    }
    comm.barrier ();
    sendRequests.clear ();
    recvRequests.clear ();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJVSendBuffer
    Array<ArrayView<Details::CrsIJV<GO, Scalar> > > sendBuffers (numSends, null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer (cur, sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s = 0; s < numSends; ++s) {
      // we'll fake the memory management, because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<Details::CrsIJV<GO, Scalar> > tmparcp =
        arcp (sendBuffers[s].getRawPtr (), 0, sendBuffers[s].size (), false);
      sendRequests.push_back (isend<int, Details::CrsIJV<GlobalOrdinal,Scalar> > (comm, tmparcp, sendIDs[s]));
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate (recvSizes.begin (), recvSizes.end (), 0);
    Array<Details::CrsIJV<GO, Scalar> > IJVRecvBuffer (totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<Details::CrsIJV<GO, Scalar> > > recvBuffers (numRecvs, null);
    {
      size_t cur = 0;
      for (size_t r = 0; r < numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer (cur, recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r = 0; r < numRecvs ; ++r) {
      // we'll fake the memory management, because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<Details::CrsIJV<GO, Scalar> > tmparcp =
        arcp (recvBuffers[r].getRawPtr (), 0, recvBuffers[r].size (), false);
      recvRequests.push_back (ireceive (comm, tmparcp, recvIDs[r]));
    }
    // perform waits
    if (! sendRequests.empty ()) {
      waitAll (comm, sendRequests ());
    }
    if (! recvRequests.empty ()) {
      waitAll (comm, recvRequests ());
    }
    comm.barrier ();
    sendRequests.clear ();
    recvRequests.clear ();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.

    typedef typename Array<Details::CrsIJV<GO, Scalar> >::const_iterator ijv_iter_type;
    if (this->isStaticGraph ()) {
      for (ijv_iter_type ijv = IJVRecvBuffer.begin ();
           ijv != IJVRecvBuffer.end (); ++ijv) {
        sumIntoGlobalValues (ijv->i, tuple (ijv->j), tuple (ijv->v));
      }
    }
    else { // Dynamic graph; can use insertGlobalValues ()
      for (ijv_iter_type ijv = IJVRecvBuffer.begin ();
           ijv != IJVRecvBuffer.end (); ++ijv) {
        try {
          insertGlobalValues (ijv->i, tuple (ijv->j), tuple (ijv->v));
        }
        catch (std::runtime_error &e) {
          std::ostringstream outmsg;
          outmsg << e.what() << std::endl
                 << "caught in globalAssemble() in " << __FILE__ << ":" << __LINE__
                 << std::endl ;
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, outmsg.str());
        }
      }
    }

    // WHEW! THAT WAS TIRING!
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    if (! isStaticGraph ()) { // Don't resume fill of a nonowned graph.
      myGraph_->resumeFill (params);
    }
    clearGlobalConstants ();
    fillComplete_ = false;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  computeGlobalConstants ()
  {
    // This method doesn't do anything.  The analogous method in
    // CrsGraph does actually compute something.
    //
    // Oddly enough, clearGlobalConstants() clears frobNorm_ (by
    // setting it to -1), but computeGlobalConstants() does _not_
    // compute the Frobenius norm; this is done on demand in
    // getFrobeniusNorm(), and the result is cached there.
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  clearGlobalConstants () {
    // We use -1 to indicate that the Frobenius norm needs to be
    // recomputed, since the values might change between now and the
    // next fillComplete call.
    //
    // Oddly enough, clearGlobalConstants() clears frobNorm_, but
    // computeGlobalConstants() does _not_ compute the Frobenius norm;
    // this is done on demand in getFrobeniusNorm(), and the result is
    // cached there.
    frobNorm_ = -STM::one ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  fillComplete (const RCP<ParameterList>& params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      getCrsGraph ().is_null (), std::logic_error, "Tpetra::CrsMatrix::"
      "fillComplete(params): getCrsGraph() returns null.  "
      "This should not happen at this point.  "
      "Please report this bug to the Tpetra developers.");

    if (isStaticGraph () && getCrsGraph ()->isFillComplete ()) {
      fillComplete (getCrsGraph ()->getDomainMap (),
                    getCrsGraph ()->getRangeMap (), params);
    } else {
      fillComplete (getRowMap (), getRowMap (), params);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                const Teuchos::RCP<const map_type>& rangeMap,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    const char tfecfFuncName[] = "fillComplete";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive () || isFillComplete (),
      std::runtime_error, ": Matrix fill state must be active (isFillActive() "
      "must be true) before you may call fillComplete().");
    const int numProcs = getComm ()->getSize ();

    //
    // Read parameters from the input ParameterList.
    //

    // If true, the caller promises that no process did nonlocal
    // changes since the last call to fillComplete.
    bool assertNoNonlocalInserts = false;
    // If true, makeColMap sorts remote GIDs (within each remote
    // process' group).
    bool sortGhosts = true;

    if (! params.is_null ()) {
      assertNoNonlocalInserts = params->get ("No Nonlocal Changes",
                                             assertNoNonlocalInserts);
      if (params->isParameter ("sort column map ghost gids")) {
        sortGhosts = params->get ("sort column map ghost gids", sortGhosts);
      }
      else if (params->isParameter ("Sort column Map ghost GIDs")) {
        sortGhosts = params->get ("Sort column Map ghost GIDs", sortGhosts);
      }
    }
    // We also don't need to do global assembly if there is only one
    // process in the communicator.
    const bool needGlobalAssemble = ! assertNoNonlocalInserts && numProcs > 1;
    // This parameter only matters if this matrix owns its graph.
    if (! myGraph_.is_null ()) {
      myGraph_->sortGhostsAssociatedWithEachProcessor_ = sortGhosts;
    }

    if (! getCrsGraph()->indicesAreAllocated()) {
      if (hasColMap ()) {
        // We have a column Map, so use local indices.
        allocateValues (LocalIndices, GraphNotYetAllocated);
      } else {
        // We don't have a column Map, so use global indices.
        allocateValues (GlobalIndices, GraphNotYetAllocated);
      }
    }
    // Global assemble, if we need to.  This call only costs a single
    // all-reduce if we didn't need global assembly after all.
    if (needGlobalAssemble) {
      globalAssemble ();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numProcs == 1 && nonlocals_.size() > 0,
        std::runtime_error, ": cannot have nonlocal entries on a serial run.  "
        "An invalid entry (i.e., with row index not in the row Map) must have "
        "been submitted to the CrsMatrix.");
    }

    if (isStaticGraph ()) {
      // FIXME (mfh 18 Jun 2014) This check for correctness of the
      // input Maps incurs a penalty of two all-reduces for the
      // otherwise optimal const graph case.
      //
      // We could turn these (max) 2 all-reduces into (max) 1, by
      // fusing them.  We could do this by adding a "locallySameAs"
      // method to Map, which would return one of four states:
      //
      //   a. Certainly globally the same
      //   b. Certainly globally not the same
      //   c. Locally the same
      //   d. Locally not the same
      //
      // The first two states don't require further communication.
      // The latter two states require an all-reduce to communicate
      // globally, but we only need one all-reduce, since we only need
      // to check whether at least one of the Maps is wrong.
      const bool domainMapsMatch = staticGraph_->getDomainMap ()->isSameAs (*domainMap);
      const bool rangeMapsMatch = staticGraph_->getRangeMap ()->isSameAs (*rangeMap);

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! domainMapsMatch, std::runtime_error,
        ": The CrsMatrix's domain Map does not match the graph's domain Map.  "
        "The graph cannot be changed because it was given to the CrsMatrix "
        "constructor as const.  You can fix this by passing in the graph's "
        "domain Map and range Map to the matrix's fillComplete call.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! rangeMapsMatch, std::runtime_error,
        ": The CrsMatrix's range Map does not match the graph's range Map.  "
        "The graph cannot be changed because it was given to the CrsMatrix "
        "constructor as const.  You can fix this by passing in the graph's "
        "domain Map and range Map to the matrix's fillComplete call.");
    }
    else {
      // Set the graph's domain and range Maps.  This will clear the
      // Import if the domain Map has changed (is a different
      // pointer), and the Export if the range Map has changed (is a
      // different pointer).
      myGraph_->setDomainRangeMaps (domainMap, rangeMap);

      // Make the graph's column Map, if necessary.
      if (! myGraph_->hasColMap ()) {
        myGraph_->makeColMap ();
      }

      // Make indices local, if necessary.  The method won't do
      // anything if the graph is already locally indexed.
      myGraph_->makeIndicesLocal ();

      if (! myGraph_->isSorted ()) {
        sortEntries ();
      }
      if (! myGraph_->isMerged ()) {
        mergeRedundantEntries ();
      }
      // Make the Import and Export, if they haven't been made already.
      myGraph_->makeImportExport ();
      myGraph_->computeGlobalConstants ();
      myGraph_->fillComplete_ = true;
      myGraph_->checkInternalState ();
    }
    computeGlobalConstants ();
    // fill local objects; will fill and finalize local graph if appropriate
    if (myGraph_.is_null ()) {
      // The matrix does _not_ own the graph, and the graph's
      // structure is already fixed, so just fill the local matrix.
      fillLocalMatrix (params);
    } else {
      // The matrix _does_ own the graph, so fill the local graph at
      // the same time as the local matrix.
      fillLocalGraphAndMatrix (params);
    }

    // Once we've initialized the sparse kernels, we're done with the
    // local objects.  We may now release them and their memory, since
    // they will persist in the local sparse ops if necessary.  We
    // keep the local graph if the parameters tell us to do so.

    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    fillComplete_ = true; // Now we're fill complete!
    checkInternalState ();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  expertStaticFillComplete (const Teuchos::RCP<const map_type> & domainMap,
                            const Teuchos::RCP<const map_type> & rangeMap,
                            const Teuchos::RCP<const import_type>& importer,
                            const Teuchos::RCP<const export_type>& exporter,
                            const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    const char tfecfFuncName[] = "expertStaticFillComplete: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, "Matrix fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::logic_error, "myGraph_ is null.  This is not allowed.");

    // We will presume globalAssemble is not needed, so we do the ESFC on the graph
    myGraph_->expertStaticFillComplete (domainMap, rangeMap, importer, exporter);

    computeGlobalConstants ();

    // Fill the local graph and matrix
    fillLocalGraphAndMatrix (params);

    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    // Now we're fill complete!
    fillComplete_ = true;

    // Sanity checks at the end.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillActive(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillComplete(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    checkInternalState();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sortEntries ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      isStaticGraph (), std::runtime_error, "Tpetra::CrsMatrix::sortEntries: "
      "Cannot sort with static graph.");
    if (! myGraph_->isSorted ()) {
      const size_t lclNumRows = this->getNodeNumRows ();
      for (size_t row = 0; row < lclNumRows; ++row) {
        RowInfo rowInfo = myGraph_->getRowInfo (row);
        Teuchos::ArrayView<impl_scalar_type> rv = this->getViewNonConst (rowInfo);
        myGraph_->template sortRowIndicesAndValues<impl_scalar_type> (rowInfo, rv);
      }
      // we just sorted every row
      myGraph_->indicesAreSorted_ = true;
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  mergeRedundantEntries ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      isStaticGraph (), std::runtime_error, "Tpetra::CrsMatrix::"
      "mergeRedundantEntries: Cannot merge with static graph.");
    if (! myGraph_->isMerged ()) {
      const size_t lclNumRows = this->getNodeNumRows ();
      for (size_t row = 0; row < lclNumRows; ++row) {
        RowInfo rowInfo = myGraph_->getRowInfo (row);
        Teuchos::ArrayView<impl_scalar_type> rv = this->getViewNonConst (rowInfo);
        myGraph_->template mergeRowIndicesAndValues<impl_scalar_type> (rowInfo, rv);
      }
      myGraph_->noRedundancies_ = true; // we just merged every row
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  applyNonTranspose (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,node_type> & X_in,
                     MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,node_type> & Y_in,
                     Scalar alpha,
                     Scalar beta) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();
    const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one ();

    // mfh 05 Jun 2014: Special case for alpha == 0.  I added this to
    // fix an Ifpack2 test (RILUKSingleProcessUnitTests), which was
    // failing only for the Kokkos refactor version of Tpetra.  It's a
    // good idea regardless to have the bypass.
    if (alpha == ZERO) {
      if (beta == ZERO) {
        Y_in.putScalar (ZERO);
      } else if (beta != ONE) {
        Y_in.scale (beta);
      }
      return;
    }

    // It's possible that X is a view of Y or vice versa.  We don't
    // allow this (apply() requires that X and Y not alias one
    // another), but it's helpful to detect and work around this case.
    // We don't try to to detect the more subtle cases (e.g., one is a
    // subview of the other, but their initial pointers differ).  We
    // only need to do this if this matrix's Import is trivial;
    // otherwise, we don't actually apply the operator from X into Y.

    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();

    // If beta == 0, then the output MV will be overwritten; none of
    // its entries should be read.  (Sparse BLAS semantics say that we
    // must ignore any Inf or NaN entries in Y_in, if beta is zero.)
    // This matters if we need to do an Export operation; see below.
    const bool Y_is_overwritten = (beta == ZERO);

    // We treat the case of a replicated MV output specially.
    const bool Y_is_replicated = ! Y_in.isDistributed ();

    // This is part of the special case for replicated MV output.
    // We'll let each process do its thing, but do an all-reduce at
    // the end to sum up the results.  Setting beta=0 on all processes
    // but Proc 0 makes the math work out for the all-reduce.  (This
    // assumes that the replicated data is correctly replicated, so
    // that the data are the same on all processes.)
    if (Y_is_replicated && this->getComm ()->getRank () > 0) {
      beta = ZERO;
    }

    // Temporary MV for Import operation.  After the block of code
    // below, this will be an (Imported if necessary) column Map MV
    // ready to give to localMultiply().
    RCP<const MV> X_colMap;
    if (importer.is_null ()) {
      if (! X_in.isConstantStride ()) {
        // Not all sparse mat-vec kernels can handle an input MV with
        // nonconstant stride correctly, so we have to copy it in that
        // case into a constant stride MV.  To make a constant stride
        // copy of X_in, we force creation of the column (== domain)
        // Map MV (if it hasn't already been created, else fetch the
        // cached copy).  This avoids creating a new MV each time.
        RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in, true);
        Tpetra::deep_copy (*X_colMapNonConst, X_in);
        X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
      }
      else {
        // The domain and column Maps are the same, so do the local
        // multiply using the domain Map input MV X_in.
        X_colMap = rcpFromRef (X_in);
      }
    }
    else {
      // We're doing an Import anyway, which will copy the relevant
      // elements of the domain Map MV X_in into a separate column Map
      // MV.  Thus, we don't have to worry whether X_in is constant
      // stride.
      RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in);

      // Import from the domain Map MV to the column Map MV.
      X_colMapNonConst->doImport (X_in, *importer, INSERT);
      X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
    }

    // Temporary MV for Export operation, or for copying a nonconstant
    // stride output MV into a constant stride MV.
    RCP<MV> Y_rowMap = getRowMapMultiVector (Y_in);

    // If we have a nontrivial Export object, we must perform an
    // Export.  In that case, the local multiply result will go into
    // the row Map multivector.  We don't have to make a
    // constant-stride version of Y_in in this case, because we had to
    // make a constant stride Y_rowMap MV and do an Export anyway.
    if (! exporter.is_null ()) {
      this->template localMultiply<Scalar, Scalar> (*X_colMap, *Y_rowMap,
                                                    Teuchos::NO_TRANS,
                                                    alpha, ZERO);
      // If we're overwriting the output MV Y_in completely (beta ==
      // 0), then make sure that it is filled with zeros before we do
      // the Export.  Otherwise, the ADD combine mode will use data in
      // Y_in, which is supposed to be zero.
      if (Y_is_overwritten) {
        Y_in.putScalar (ZERO);
      }
      else {
        // Scale the output MV by beta, so that the Export sums in the
        // mat-vec contribution: Y_in = beta*Y_in + alpha*A*X_in.
        Y_in.scale (beta);
      }
      // Do the Export operation.
      Y_in.doExport (*Y_rowMap, *exporter, ADD);
    }
    else { // Don't do an Export: row Map and range Map are the same.
      //
      // If Y_in does not have constant stride, or if the column Map
      // MV aliases Y_in, then we can't let the kernel write directly
      // to Y_in.  Instead, we have to use the cached row (== range)
      // Map MV as temporary storage.
      //
      // FIXME (mfh 05 Jun 2014) This test for aliasing only tests if
      // the user passed in the same MultiVector for both X and Y.  It
      // won't detect whether one MultiVector views the other.  We
      // should also check the MultiVectors' raw data pointers.
      if (! Y_in.isConstantStride () || X_colMap.getRawPtr () == &Y_in) {
        // Force creating the MV if it hasn't been created already.
        // This will reuse a previously created cached MV.
        Y_rowMap = getRowMapMultiVector (Y_in, true);

        // If beta == 0, we don't need to copy Y_in into Y_rowMap,
        // since we're overwriting it anyway.
        if (beta != ZERO) {
          Tpetra::deep_copy (*Y_rowMap, Y_in);
        }
        this->template localMultiply<Scalar, Scalar> (*X_colMap,
                                                      *Y_rowMap,
                                                      Teuchos::NO_TRANS,
                                                      alpha, beta);
        Tpetra::deep_copy (Y_in, *Y_rowMap);
      }
      else {
        this->template localMultiply<Scalar, Scalar> (*X_colMap, Y_in,
                                                      Teuchos::NO_TRANS,
                                                      alpha, beta);
      }
    }

    // If the range Map is a locally replicated Map, sum up
    // contributions from each process.  We set beta = 0 on all
    // processes but Proc 0 initially, so this will handle the scaling
    // factor beta correctly.
    if (Y_is_replicated) {
      Y_in.reduce ();
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal, class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  applyTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X_in,
                  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& Y_in,
                  const Teuchos::ETransp mode,
                  Scalar alpha,
                  Scalar beta) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    // Take shortcuts for alpha == 0.
    if (alpha == ZERO) {
      // Follow the Sparse BLAS convention by ignoring both the matrix
      // and X_in, in this case.
      if (beta == ZERO) {
        // Follow the Sparse BLAS convention by overwriting any Inf or
        // NaN values in Y_in, in this case.
        Y_in.putScalar (ZERO);
      }
      else {
        Y_in.scale (beta);
      }
      return;
    }

    const size_t numVectors = X_in.getNumVectors ();

    // We don't allow X_in and Y_in to alias one another.  It's hard
    // to check this, because advanced users could create views from
    // raw pointers.  However, if X_in and Y_in reference the same
    // object, we will do the user a favor by copying X into new
    // storage (with a warning).  We only need to do this if we have
    // trivial importers; otherwise, we don't actually apply the
    // operator from X into Y.
    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    // access X indirectly, in case we need to create temporary storage
    RCP<const MV> X;

    // some parameters for below
    const bool Y_is_replicated = ! Y_in.isDistributed ();
    const bool Y_is_overwritten = (beta == ZERO);
    if (Y_is_replicated && this->getComm ()->getRank () > 0) {
      beta = ZERO;
    }

    // The kernels do not allow input or output with nonconstant stride.
    if (! X_in.isConstantStride () && importer.is_null ()) {
      X = rcp (new MV (X_in)); // Constant-stride copy of X_in
    } else {
      X = rcpFromRef (X_in); // Reference to X_in
    }

    // Set up temporary multivectors for Import and/or Export.
    if (importer != null) {
      if (importMV_ != null && importMV_->getNumVectors() != numVectors) {
        importMV_ = null;
      }
      if (importMV_ == null) {
        importMV_ = rcp (new MV (this->getColMap (), numVectors));
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) {
        exportMV_ = null;
      }
      if (exportMV_ == null) {
        exportMV_ = rcp (new MV (this->getRowMap (), numVectors));
      }
    }

    // If we have a non-trivial exporter, we must import elements that
    // are permuted or are on other processors.
    if (! exporter.is_null ()) {
      exportMV_->doImport (X_in, *exporter, INSERT);
      X = exportMV_; // multiply out of exportMV_
    }

    // If we have a non-trivial importer, we must export elements that
    // are permuted or belong to other processors.  We will compute
    // solution into the to-be-exported MV; get a view.
    if (importer != null) {
      // Do the local computation.
      this->template localMultiply<Scalar, Scalar> (*X, *importMV_, mode,
                                                    alpha, ZERO);
      if (Y_is_overwritten) {
        Y_in.putScalar (ZERO);
      } else {
        Y_in.scale (beta);
      }
      Y_in.doExport (*importMV_, *importer, ADD);
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't multiply into non-strided multivector
      //
      // FIXME (mfh 05 Jun 2014) This test for aliasing only tests if
      // the user passed in the same MultiVector for both X and Y.  It
      // won't detect whether one MultiVector views the other.  We
      // should also check the MultiVectors' raw data pointers.
      if (! Y_in.isConstantStride () || X.getRawPtr () == &Y_in) {
        // Make a deep copy of Y_in, into which to write the multiply result.
        MV Y (Y_in, Teuchos::Copy);
        this->template localMultiply<Scalar, Scalar> (*X, Y, mode, alpha, beta);
        Tpetra::deep_copy (Y_in, Y);
      } else {
        this->template localMultiply<Scalar, Scalar> (*X, Y_in, mode, alpha, beta);
      }
    }

    // If the range Map is a locally replicated map, sum the
    // contributions from each process.  (That's why we set beta=0
    // above for all processes but Proc 0.)
    if (Y_is_replicated) {
      Y_in.reduce ();
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &X,
         MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isFillComplete (), std::runtime_error,
      "Tpetra::CrsMatrix::apply(): Cannot call apply() until fillComplete() "
      "has been called.");

    if (mode == Teuchos::NO_TRANS) {
      applyNonTranspose (X, Y, alpha, beta);
    } else {
      //Thyra was implicitly assuming that Y gets set to zero / or is overwritten
      //when bets==0. This was not the case with transpose in a multithreaded
      //environment where a multiplication with subsequent atomic_adds is used
      //since 0 is effectively not special cased. Doing the explicit set to zero here
      //This catches cases where Y is nan or inf.
      const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();
      if(beta == ZERO)
        Y.putScalar (ZERO);
      applyTranspose (X, Y, mode, alpha, beta);
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
               const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
               const Scalar& dampingFactor,
               const ESweepDirection direction,
               const int numSweeps) const
  {
    reorderedGaussSeidel (B, X, D, Teuchos::null, dampingFactor, direction, numSweeps);
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reorderedGaussSeidel (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& B,
                        MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& X,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& D,
                        const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                        const Scalar& dampingFactor,
                        const ESweepDirection direction,
                        const int numSweeps) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    typedef Scalar ST;

    TEUCHOS_TEST_FOR_EXCEPTION(
      isFillComplete() == false, std::runtime_error,
      "Tpetra::CrsMatrix::gaussSeidel: cannot call this method until "
      "fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0,
      std::invalid_argument,
      "Tpetra::CrsMatrix::gaussSeidel: The number of sweeps must be , "
      "nonnegative but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    KokkosClassic::ESweepDirection localDirection;
    if (direction == Forward) {
      localDirection = KokkosClassic::Forward;
    }
    else if (direction == Backward) {
      localDirection = KokkosClassic::Backward;
    }
    else if (direction == Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = KokkosClassic::Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Tpetra::CrsMatrix::gaussSeidel: The 'direction' enum does not have "
        "any of its valid values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return; // Nothing to do.
    }

    // We don't need the Export object because this method assumes
    // that the row, domain, and range Maps are the same.  We do need
    // the Import object, if there is one, though.
    RCP<const import_type> importer = this->getGraph()->getImporter();
    RCP<const export_type> exporter = this->getGraph()->getExporter();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! exporter.is_null (), std::runtime_error,
      "Tpetra's gaussSeidel implementation requires that the row, domain, "
      "and range Maps be the same.  This cannot be the case, because the "
      "matrix has a nontrivial Export object.");

    RCP<const map_type> domainMap = this->getDomainMap ();
    RCP<const map_type> rangeMap = this->getRangeMap ();
    RCP<const map_type> rowMap = this->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = this->getGraph ()->getColMap ();

#ifdef HAVE_TEUCHOS_DEBUG
    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }
#else
    // Forestall any compiler warnings for unused variables.
    (void) rangeMap;
    (void) rowMap;
#endif // HAVE_TEUCHOS_DEBUG

    // If B is not constant stride, copy it into a constant stride
    // multivector.  We'l handle the right-hand side B first and deal
    // with X right before the sweeps, to improve locality of the
    // first sweep.  (If the problem is small enough, then that will
    // hopefully keep more of the entries of X in cache.  This
    // optimizes for the typical case of a small number of sweeps.)
    RCP<const MV> B_in;
    if (B.isConstantStride()) {
      B_in = rcpFromRef (B);
    }
    else {
      // The range Map and row Map are the same in this case, so we
      // can use the (possibly cached) row Map multivector to store a
      // constant stride copy of B.  We don't have to copy back, since
      // Gauss-Seidel won't modify B.
      RCP<MV> B_in_nonconst = getRowMapMultiVector (B, true);
      deep_copy (*B_in_nonconst, B); // Copy from B into B_in(_nonconst).
      B_in = rcp_const_cast<const MV> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
        "requires that X and B both have constant stride.  Since B does not "
        "have constant stride, we had to make a copy.  This is a limitation of "
        "the current implementation and not your fault, but we still report it "
        "as an efficiency warning for your information.");
    }

    // If X is not constant stride, copy it into a constant stride
    // multivector.  Also, make the column Map multivector X_colMap,
    // and its domain Map view X_domainMap.  (X actually must be a
    // domain Map view of a column Map multivector; exploit this, if X
    // has constant stride.)

    RCP<MV> X_domainMap;
    RCP<MV> X_colMap;
    bool copiedInput = false;

    if (importer.is_null ()) { // Domain and column Maps are the same.
      if (X.isConstantStride ()) {
        X_domainMap = rcpFromRef (X);
        X_colMap = X_domainMap;
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector, make a domain Map
        // view of it, and copy X into the domain Map view.  We have
        // to copy here because we won't be doing Import operations.
        X_colMap = getColumnMapMultiVector (X, true);
        X_domainMap = X_colMap; // Domain and column Maps are the same.
        deep_copy (*X_domainMap, X); // Copy X into the domain Map view.
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant "
          "stride.  Since X does not have constant stride, we had to make a "
          "copy.  This is a limitation of the current implementation and not "
          "your fault, but we still report it as an efficiency warning for "
          "your information.");
      }
    }
    else { // We will be doing Import operations in the sweeps.
      if (X.isConstantStride ()) {
        X_domainMap = rcpFromRef (X);
        // This kernel assumes that X is a domain Map view of a column
        // Map multivector.  We will only check if this is valid if
        // the CMake configure Teuchos_ENABLE_DEBUG is ON.
        X_colMap = X_domainMap->offsetViewNonConst (colMap, 0);

        // FIXME (mfh 19 Mar 2013) Do we need to fill the remote
        // entries of X_colMap with zeros?  Do we need to fill all of
        // X_domainMap initially with zeros?  Ifpack
        // (Ifpack_PointRelaxation.cpp, line 906) creates an entirely
        // new MultiVector each time.

        // Do the first Import for the first sweep.  This simplifies
        // the logic in the sweeps.
        X_colMap->doImport (X, *importer, INSERT);
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector X_colMap, and make a
        // domain Map view X_domainMap of it.  Instead of copying, we
        // do an Import from X into X_domainMap.  This saves us a
        // copy, since the Import has to copy the data anyway.
        X_colMap = getColumnMapMultiVector (X, true);
        X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);
        X_colMap->doImport (X, *importer, INSERT);
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant stride.  "
          "Since X does not have constant stride, we had to make a copy.  "
          "This is a limitation of the current implementation and not your fault, "
          "but we still report it as an efficiency warning for your information.");
      }
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep.
        X_colMap->doImport (*X_domainMap, *importer, INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Symmetric) {
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   localDirection);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            localDirection);
        }
      }
      else { // direction == Symmetric
        const bool doImportBetweenDirections = false;
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   KokkosClassic::Forward);
          // mfh 18 Mar 2013: Aztec's implementation of "symmetric
          // Gauss-Seidel" does _not_ do an Import between the forward
          // and backward sweeps.  This makes sense, because Aztec
          // considers "symmetric Gauss-Seidel" a subdomain solver.
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, INSERT);
            }
          }
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   KokkosClassic::Backward);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            KokkosClassic::Forward);
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, INSERT);
            }
          }
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            KokkosClassic::Backward);
        }
      }
    }

    if (copiedInput) {
      deep_copy (X, *X_domainMap); // Copy back from X_domainMap to X.
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  gaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                   const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                   const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps,
                   const bool zeroInitialGuess) const
  {
    reorderedGaussSeidelCopy (X, B, D, Teuchos::null, dampingFactor, direction,
                              numSweeps, zeroInitialGuess);
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reorderedGaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                            const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                            const Scalar& dampingFactor,
                            const ESweepDirection direction,
                            const int numSweeps,
                            const bool zeroInitialGuess) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::rcp_const_cast;
    typedef Scalar ST;
    const char prefix[] = "Tpetra::CrsMatrix::(reordered)gaussSeidelCopy: ";
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isFillComplete (), std::runtime_error,
      prefix << "The matrix is not fill complete.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0, std::invalid_argument,
      prefix << "The number of sweeps must be nonnegative, "
      "but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    KokkosClassic::ESweepDirection localDirection;
    if (direction == Forward) {
      localDirection = KokkosClassic::Forward;
    }
    else if (direction == Backward) {
      localDirection = KokkosClassic::Backward;
    }
    else if (direction == Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = KokkosClassic::Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        prefix << "The 'direction' enum does not have any of its valid "
        "values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return;
    }

    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! exporter.is_null (), std::runtime_error,
      "This method's implementation currently requires that the matrix's row, "
      "domain, and range Maps be the same.  This cannot be the case, because "
      "the matrix has a nontrivial Export object.");

    RCP<const map_type> domainMap = this->getDomainMap ();
    RCP<const map_type> rangeMap = this->getRangeMap ();
    RCP<const map_type> rowMap = this->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = this->getGraph ()->getColMap ();

#ifdef HAVE_TEUCHOS_DEBUG
    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }
#else
    // Forestall any compiler warnings for unused variables.
    (void) rangeMap;
    (void) rowMap;
#endif // HAVE_TEUCHOS_DEBUG

    // Fetch a (possibly cached) temporary column Map multivector
    // X_colMap, and a domain Map view X_domainMap of it.  Both have
    // constant stride by construction.  We know that the domain Map
    // must include the column Map, because our Gauss-Seidel kernel
    // requires that the row Map, domain Map, and range Map are all
    // the same, and that each process owns all of its own diagonal
    // entries of the matrix.

    RCP<MV> X_colMap;
    RCP<MV> X_domainMap;
    bool copyBackOutput = false;
    if (importer.is_null ()) {
      if (X.isConstantStride ()) {
        X_colMap = rcpFromRef (X);
        X_domainMap = rcpFromRef (X);
        // Column Map and domain Map are the same, so there are no
        // remote entries.  Thus, if we are not setting the initial
        // guess to zero, we don't have to worry about setting remote
        // entries to zero, even though we are not doing an Import in
        // this case.
        if (zeroInitialGuess) {
          X_colMap->putScalar (ZERO);
        }
        // No need to copy back to X at end.
      }
      else { // We must copy X into a constant stride multivector.
        // Just use the cached column Map multivector for that.
        // force=true means fill with zeros, so no need to fill
        // remote entries (not in domain Map) with zeros.
        X_colMap = getColumnMapMultiVector (X, true);
        // X_domainMap is always a domain Map view of the column Map
        // multivector.  In this case, the domain and column Maps are
        // the same, so X_domainMap _is_ X_colMap.
        X_domainMap = X_colMap;
        if (! zeroInitialGuess) { // Don't copy if zero initial guess
          try {
            deep_copy (*X_domainMap , X); // Copy X into constant stride MV
          } catch (std::exception& e) {
            std::ostringstream os;
            os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
              "deep_copy(*X_domainMap, X) threw an exception: "
               << e.what () << ".";
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
          }
        }
        copyBackOutput = true; // Don't forget to copy back at end.
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (),
          std::runtime_error,
          "gaussSeidelCopy: The current implementation of the Gauss-Seidel "
          "kernel requires that X and B both have constant stride.  Since X "
          "does not have constant stride, we had to make a copy.  This is a "
          "limitation of the current implementation and not your fault, but we "
          "still report it as an efficiency warning for your information.");
      }
    }
    else { // Column Map and domain Map are _not_ the same.
      X_colMap = getColumnMapMultiVector (X);
      X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);

#ifdef HAVE_TPETRA_DEBUG
      typename MV::dual_view_type X_colMap_view = X_colMap->getDualView ();
      typename MV::dual_view_type X_domainMap_view = X_domainMap->getDualView ();

      if (X_colMap->getLocalLength () != 0 && X_domainMap->getLocalLength ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          X_colMap_view.h_view.ptr_on_device () != X_domainMap_view.h_view.ptr_on_device (),
          std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
          "Pointer to start of column Map view of X is not equal to pointer to "
          "start of (domain Map view of) X.  This may mean that "
          "Tpetra::MultiVector::offsetViewNonConst is broken.  "
          "Please report this bug to the Tpetra developers.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap_view.dimension_0 () < X_domainMap_view.dimension_0 () ||
        X_colMap->getLocalLength () < X_domainMap->getLocalLength (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has fewer local rows than X_domainMap.  "
        "X_colMap_view.dimension_0() = " << X_colMap_view.dimension_0 ()
        << ", X_domainMap_view.dimension_0() = "
        << X_domainMap_view.dimension_0 ()
        << ", X_colMap->getLocalLength() = " << X_colMap->getLocalLength ()
        << ", and X_domainMap->getLocalLength() = "
        << X_domainMap->getLocalLength ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap->getNumVectors () != X_domainMap->getNumVectors (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has a different number of columns than X_domainMap.  "
        "X_colMap->getNumVectors() = " << X_colMap->getNumVectors ()
        << " != X_domainMap->getNumVectors() = "
        << X_domainMap->getNumVectors ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");

      // TEUCHOS_TEST_FOR_EXCEPTION(
      //   X_colMap->getLocalMV ().getStride () !=
      //   X_domainMap->getLocalMV ().getStride (),
      //   std::logic_error,
      //   "Tpetra::CrsMatrix::gaussSeidelCopy: "
      //   "X_colMap has local stride " << X_colMap->getLocalMV ().getStride ()
      //   << ", which does not equal the local stride "
      //   << X_domainMap->getLocalMV ().getStride () << " of X_domainMap.  "
      //   "This means that Tpetra::MultiVector::offsetViewNonConst is broken.  "
      //   "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      if (zeroInitialGuess) {
        // No need for an Import, since we're filling with zeros.
        X_colMap->putScalar (ZERO);
      } else {
        // We could just copy X into X_domainMap.  However, that
        // wastes a copy, because the Import also does a copy (plus
        // communication).  Since the typical use case for
        // Gauss-Seidel is a small number of sweeps (2 is typical), we
        // don't want to waste that copy.  Thus, we do the Import
        // here, and skip the first Import in the first sweep.
        // Importing directly from X effects the copy into X_domainMap
        // (which is a view of X_colMap).
        X_colMap->doImport (X, *importer, INSERT);
      }
      copyBackOutput = true; // Don't forget to copy back at end.
    } // if column and domain Maps are (not) the same

    // The Gauss-Seidel / SOR kernel expects multivectors of constant
    // stride.  X_colMap is by construction, but B might not be.  If
    // it's not, we have to make a copy.
    RCP<const MV> B_in;
    if (B.isConstantStride ()) {
      B_in = rcpFromRef (B);
    }
    else {
      // Range Map and row Map are the same in this case, so we can
      // use the cached row Map multivector to store a constant stride
      // copy of B.
      RCP<MV> B_in_nonconst = getRowMapMultiVector (B, true);
      try {
        deep_copy (*B_in_nonconst, B);
      } catch (std::exception& e) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
          "deep_copy(*B_in_nonconst, B) threw an exception: "
           << e.what () << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
      }
      B_in = rcp_const_cast<const MV> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidelCopy: The current implementation requires that B have "
        "constant stride.  Since B does not have constant stride, we had to "
        "copy it into a separate constant-stride multivector.  This is a "
        "limitation of the current implementation and not your fault, but we "
        "still report it as an efficiency warning for your information.");
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep above,
        // if it was necessary.
        X_colMap->doImport (*X_domainMap, *importer, INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Symmetric) {
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   localDirection);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            localDirection);
        }
      }
      else { // direction == Symmetric
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   KokkosClassic::Forward);
          // mfh 18 Mar 2013: Aztec's implementation of "symmetric
          // Gauss-Seidel" does _not_ do an Import between the forward
          // and backward sweeps.  This makes symmetric Gauss-Seidel a
          // symmetric preconditioner if the matrix A is symmetric.  We
          // imitate Aztec's behavior here.
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   KokkosClassic::Backward);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            KokkosClassic::Forward);
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            KokkosClassic::Backward);

        }
      }
    }

    if (copyBackOutput) {
      try {
        deep_copy (X , *X_domainMap); // Copy result back into X.
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, prefix << "deep_copy(X, *X_domainMap) "
          "threw an exception: " << e.what ());
      }
    }
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  template <class DomainScalar, class RangeScalar>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  localMultiply (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                 MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& Y,
                 Teuchos::ETransp mode,
                 RangeScalar alpha,
                 RangeScalar beta) const
  {
    using Teuchos::NO_TRANS;
    // Just like Scalar and impl_scalar_type may differ in CrsMatrix,
    // RangeScalar and its corresponding impl_scalar_type may differ in
    // MultiVector.
    typedef typename MultiVector<RangeScalar, LocalOrdinal, GlobalOrdinal,
                                 node_type>::impl_scalar_type range_impl_scalar_type;
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "localMultiply: ";
#endif // HAVE_TPETRA_DEBUG

    const range_impl_scalar_type theAlpha = static_cast<range_impl_scalar_type> (alpha);
    const range_impl_scalar_type theBeta = static_cast<range_impl_scalar_type> (beta);
    const bool conjugate = (mode == Teuchos::CONJ_TRANS);
    const bool transpose = (mode != Teuchos::NO_TRANS);

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
      "X.getNumVectors() = " << X.getNumVectors () << " != Y.getNumVectors() = "
      << Y.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! transpose && X.getLocalLength () != getColMap ()->getNodeNumElements (),
      std::runtime_error, "NO_TRANS case: X has the wrong number of local rows.  "
      "X.getLocalLength() = " << X.getLocalLength () << " != getColMap()->"
      "getNodeNumElements() = " << getColMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! transpose && Y.getLocalLength () != getRowMap ()->getNodeNumElements (),
      std::runtime_error, "NO_TRANS case: Y has the wrong number of local rows.  "
      "Y.getLocalLength() = " << Y.getLocalLength () << " != getRowMap()->"
      "getNodeNumElements() = " << getRowMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      transpose && X.getLocalLength () != getRowMap ()->getNodeNumElements (),
      std::runtime_error, "TRANS or CONJ_TRANS case: X has the wrong number of "
      "local rows.  X.getLocalLength() = " << X.getLocalLength () << " != "
      "getRowMap()->getNodeNumElements() = "
      << getRowMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      transpose && Y.getLocalLength () != getColMap ()->getNodeNumElements (),
      std::runtime_error, "TRANS or CONJ_TRANS case: X has the wrong number of "
      "local rows.  Y.getLocalLength() = " << Y.getLocalLength () << " != "
      "getColMap()->getNodeNumElements() = "
      << getColMap ()->getNodeNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete (), std::runtime_error, "The matrix is not fill "
      "complete.  You must call fillComplete() (possibly with domain and range "
      "Map arguments) without an intervening resumeFill() call before you may "
      "call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! X.isConstantStride () || ! Y.isConstantStride (), std::runtime_error,
      "X and Y must be constant stride.");
    // If the two pointers are NULL, then they don't alias one
    // another, even though they are equal.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getDualView ().d_view.ptr_on_device () == Y.getDualView ().d_view.ptr_on_device () &&
      X.getDualView ().d_view.ptr_on_device () != NULL,
      std::runtime_error, "X and Y may not alias one another.");
#endif // HAVE_TPETRA_DEBUG

      // Y = alpha*op(M) + beta*Y
      if (transpose) {
        KokkosSparse::spmv ( conjugate?KokkosSparse::ConjugateTranspose:KokkosSparse::Transpose,
                             theAlpha,
                             lclMatrix_,
                             X.template getLocalView<DeviceType> (),
                             theBeta,
                             Y.template getLocalView<DeviceType> ()
                            );
      }
      else {
        KokkosSparse::spmv ( KokkosSparse::NoTranspose,
                             theAlpha,
                             lclMatrix_,
                             X.template getLocalView<DeviceType> (),
                             theBeta,
                             Y.template getLocalView<DeviceType> ()
                            );
      }

  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  template <class DomainScalar, class RangeScalar>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  localGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                    MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                    const RangeScalar& dampingFactor,
                    const KokkosClassic::ESweepDirection direction) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::MultiVector<DomainScalar, LO, GO, node_type> DMV;
    typedef Tpetra::MultiVector<RangeScalar, LO, GO, node_type> RMV;
    typedef Tpetra::MultiVector<Scalar, LO, GO, node_type> MMV;
    typedef typename DMV::dual_view_type::host_mirror_space HMDT ;
    typedef typename Graph::local_graph_type k_local_graph_type;
    typedef typename k_local_graph_type::size_type offset_type;
    const char prefix[] = "Tpetra::CrsMatrix::localGaussSeidel: ";

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->isFillComplete (), std::runtime_error,
      prefix << "The matrix is not fill complete.");
    const size_t lclNumRows = this->getNodeNumRows ();
    const size_t numVecs = B.getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors () != numVecs, std::invalid_argument,
      prefix << "B.getNumVectors() = " << numVecs << " != "
      "X.getNumVectors() = " << X.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      B.getLocalLength () != lclNumRows, std::invalid_argument,
      prefix << "B.getLocalLength() = " << B.getLocalLength ()
      << " != this->getNodeNumRows() = " << lclNumRows << ".");

    typename DMV::dual_view_type::t_host B_lcl = B.template getLocalView<HMDT> ();
    typename RMV::dual_view_type::t_host X_lcl = X.template getLocalView<HMDT> ();
    typename MMV::dual_view_type::t_host D_lcl = D.template getLocalView<HMDT> ();

    offset_type B_stride[8], X_stride[8], D_stride[8];
    B_lcl.stride (B_stride);
    X_lcl.stride (X_stride);
    D_lcl.stride (D_stride);

    local_matrix_type lclMatrix = this->getLocalMatrix ();
    k_local_graph_type lclGraph = lclMatrix.graph;
    typename local_matrix_type::row_map_type ptr = lclGraph.row_map;
    typename local_matrix_type::index_type ind = lclGraph.entries;
    typename local_matrix_type::values_type val = lclMatrix.values;
    const offset_type* const ptrRaw = ptr.ptr_on_device ();
    const LO* const indRaw = ind.ptr_on_device ();
    const impl_scalar_type* const valRaw = val.ptr_on_device ();

    Kokkos::Sequential::gaussSeidel (static_cast<LO> (lclNumRows),
                                     static_cast<LO> (numVecs),
                                     ptrRaw, indRaw, valRaw,
                                     B_lcl.ptr_on_device (), B_stride[1],
                                     X_lcl.ptr_on_device (), X_stride[1],
                                     D_lcl.ptr_on_device (),
                                     static_cast<impl_scalar_type> (dampingFactor),
                                     direction);
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  template<class DomainScalar,
           class RangeScalar>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reorderedLocalGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                             MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                             const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                             const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                             const RangeScalar& dampingFactor,
                             const KokkosClassic::ESweepDirection direction) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::MultiVector<DomainScalar, LO, GO, node_type> DMV;
    typedef Tpetra::MultiVector<RangeScalar, LO, GO, node_type> RMV;
    typedef Tpetra::MultiVector<Scalar, LO, GO, node_type> MMV;
    typedef typename DMV::dual_view_type::host_mirror_space HMDT ;
    typedef typename Graph::local_graph_type k_local_graph_type;
    typedef typename k_local_graph_type::size_type offset_type;
    const char prefix[] = "Tpetra::CrsMatrix::reorderedLocalGaussSeidel: ";

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->isFillComplete (), std::runtime_error,
      prefix << "The matrix is not fill complete.");
    const size_t lclNumRows = this->getNodeNumRows ();
    const size_t numVecs = B.getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors () != numVecs, std::invalid_argument,
      prefix << "B.getNumVectors() = " << numVecs << " != "
      "X.getNumVectors() = " << X.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      B.getLocalLength () != lclNumRows, std::invalid_argument,
      prefix << "B.getLocalLength() = " << B.getLocalLength ()
      << " != this->getNodeNumRows() = " << lclNumRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (rowIndices.size ()) < lclNumRows,
      std::invalid_argument, prefix << "rowIndices.size() = "
      << rowIndices.size () << " < this->getNodeNumRows() = "
      << lclNumRows << ".");

    typename DMV::dual_view_type::t_host B_lcl = B.template getLocalView<HMDT> ();
    typename RMV::dual_view_type::t_host X_lcl = X.template getLocalView<HMDT> ();
    typename MMV::dual_view_type::t_host D_lcl = D.template getLocalView<HMDT> ();

    offset_type B_stride[8], X_stride[8], D_stride[8];
    B_lcl.stride (B_stride);
    X_lcl.stride (X_stride);
    D_lcl.stride (D_stride);

    local_matrix_type lclMatrix = this->getLocalMatrix ();
    typename Graph::local_graph_type lclGraph = lclMatrix.graph;
    typename local_matrix_type::index_type ind = lclGraph.entries;
    typename local_matrix_type::row_map_type ptr = lclGraph.row_map;
    typename local_matrix_type::values_type val = lclMatrix.values;
    const offset_type* const ptrRaw = ptr.ptr_on_device ();
    const LO* const indRaw = ind.ptr_on_device ();
    const impl_scalar_type* const valRaw = val.ptr_on_device ();

    Kokkos::Sequential::reorderedGaussSeidel (static_cast<LO> (lclNumRows),
                                              static_cast<LO> (numVecs),
                                              ptrRaw, indRaw, valRaw,
                                              B_lcl.ptr_on_device (),
                                              B_stride[1],
                                              X_lcl.ptr_on_device (),
                                              X_stride[1],
                                              D_lcl.ptr_on_device (),
                                              rowIndices.getRawPtr (),
                                              static_cast<LO> (lclNumRows),
                                              static_cast<impl_scalar_type> (dampingFactor),
                                              direction);
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  template<class DomainScalar,
           class RangeScalar>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  localSolve (const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& Y,
              MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
              Teuchos::ETransp mode) const
  {
    using Kokkos::Sequential::triSolveKokkos;
    using Teuchos::CONJ_TRANS;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::MultiVector<DomainScalar, LO, GO, node_type> DMV;
    typedef Tpetra::MultiVector<RangeScalar, LO, GO, node_type> RMV;
    typedef typename DMV::dual_view_type::host_mirror_space HMDT ;

    const char tfecfFuncName[] = "localSolve: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete (), std::runtime_error,
      "The matrix is not fill complete.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! X.isConstantStride () || ! Y.isConstantStride (), std::invalid_argument,
      "X and Y must be constant stride.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isUpperTriangular () && ! isLowerTriangular (), std::runtime_error,
      "The matrix is neither upper triangular or lower triangular.  "
      "You may only call this method if the matrix is triangular.  "
      "Remember that this is a local (per MPI process) property, and that "
      "Tpetra only knows how to do a local (per process) triangular solve.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      STS::isComplex && mode == TRANS, std::logic_error, "This method does "
      "not currently support non-conjugated transposed solve (mode == "
      "Teuchos::TRANS) for complex scalar types.");

    // FIXME (mfh 27 Aug 2014) Tpetra has always made the odd decision
    // that if _some_ diagonal entries are missing locally, then it
    // assumes that the matrix has an implicitly stored unit diagonal.
    // Whether the matrix has an implicit unit diagonal or not should
    // be up to the user to decide.  What if the graph has no diagonal
    // entries, and the user wants it that way?  The only reason this
    // matters, though, is for the triangular solve, and in that case,
    // missing diagonal entries will cause trouble anyway.  However,
    // it would make sense to warn the user if they ask for a
    // triangular solve with an incomplete diagonal.  Furthermore,
    // this code should only assume an implicitly stored unit diagonal
    // if the matrix has _no_ explicitly stored diagonal entries.
    const Teuchos::EDiag diag = getNodeNumDiags () < getNodeNumRows () ?
      Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG;
    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    if (isUpperTriangular ()) {
      uplo = Teuchos::UPPER_TRI;
    } else if (isLowerTriangular ()) {
      uplo = Teuchos::LOWER_TRI;
    }

    local_matrix_type A_lcl = this->getLocalMatrix ();
    typename DMV::dual_view_type::t_host X_lcl = X.template getLocalView<HMDT> ();
    typename RMV::dual_view_type::t_host Y_lcl = Y.template getLocalView<HMDT> ();
    triSolveKokkos (X_lcl, A_lcl, Y_lcl, uplo, diag, mode);
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  template<class T>
  Teuchos::RCP<CrsMatrix<
                 T, LocalOrdinal, GlobalOrdinal,
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  convert () const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef CrsMatrix<T, LocalOrdinal, GlobalOrdinal, node_type> out_mat_type;
    typedef typename out_mat_type::local_matrix_type out_lcl_mat_type;
    typedef typename out_lcl_mat_type::values_type out_vals_type;
    typedef ArrayRCP<size_t>::size_type size_type;
    const char tfecfFuncName[] = "convert";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete (), std::runtime_error, "This matrix (the source of "
      "the conversion) is not fill complete.  You must first call "
      "fillComplete() (possibly with the domain and range Map) without an "
      "intervening call to resumeFill(), before you may call this method.");

    // mfh 27 Feb 2014: It seems reasonable that if this matrix has a
    // const graph, then the returned matrix should also.  However, if
    // this matrix does not have a const graph, then neither should
    // the returned matrix.  The code below implements this strategy.

    RCP<out_mat_type> newmat; // the matrix to return

    if (this->isStaticGraph ()) {
      // This matrix has a const graph, so the returned matrix should too.
      newmat = rcp (new out_mat_type (this->getCrsGraph ()));

      // Convert the values from Scalar to T, and stuff them directly
      // into the matrix to return.
      const size_type numVals =
        static_cast<size_type> (this->lclMatrix_.values.dimension_0 ());

      // FIXME (mfh 05 Aug 2014) Write a copy kernel (impl_scalar_type and
      // T differ, so we can't use Kokkos::deep_copy).
      //
      // FIXME (mfh 05 Aug 2014) This assumes UVM.
      out_vals_type newVals1D ("Tpetra::CrsMatrix::val", numVals);
      for (size_type k = 0; k < numVals; ++k) {
        newVals1D(k) = static_cast<T> (this->k_values1D_(k));
      }
      newmat->lclMatrix_ =
        out_lcl_mat_type ("Tpetra::CrsMatrix::lclMatrix_",
                          this->lclMatrix_.numCols (), newVals1D,
                          this->lclMatrix_.graph);
      newmat->k_values1D_ = newVals1D;
      // Since newmat has a static (const) graph, the graph already
      // has a column Map, and Import and Export objects already exist
      // (if applicable).  Thus, calling fillComplete is cheap.
      newmat->fillComplete (this->getDomainMap (), this->getRangeMap ());
    }
    else {
      // This matrix has a nonconst graph, so the returned matrix
      // should also have a nonconst graph.  However, it's fine for
      // the returned matrix to have static profile.  This will
      // certainly speed up its fillComplete.

      //
      // FIXME (mfh 05 Aug 2014) Instead of the slow stuff below, we
      // should copy the values and existing graph into a new local
      // matrix (lclMatrix), and then use the Tpetra::CrsMatrix
      // constructor that takes (rowMap, colMap, lclMatrix, params).
      //

      // Get this matrix's local data.
      ArrayRCP<const size_t> ptr;
      ArrayRCP<const LocalOrdinal> ind;
      ArrayRCP<const Scalar> oldVal;
      this->getAllValues (ptr, ind, oldVal);

      RCP<const map_type> rowMap = this->getRowMap ();
      RCP<const map_type> colMap = this->getColMap ();

      // Get an array of the number of entries in each (locally owned)
      // row, so that we can make the new matrix with static profile.
      const size_type numLocalRows =
        static_cast<size_type> (rowMap->getNodeNumElements ());
      ArrayRCP<size_t> numEntriesPerRow (numLocalRows);
      for (size_type localRow = 0; localRow < numLocalRows; ++localRow) {
        numEntriesPerRow[localRow] =
          static_cast<size_type> (getNumEntriesInLocalRow (localRow));
      }

      newmat = rcp (new out_mat_type (rowMap, colMap, numEntriesPerRow,
                                      StaticProfile));

      // Convert this matrix's values from Scalar to T.
      const size_type numVals = this->lclMatrix_.values.dimension_0 ();
      ArrayRCP<T> newVals1D (numVals);
      // FIXME (mfh 05 Aug 2014) This assumes UVM.
      for (size_type k = 0; k < numVals; ++k) {
        newVals1D[k] = static_cast<T> (this->k_values1D_(k));
      }

      // Give this matrix all of its local data.  We can all this
      // method because newmat was _not_ created with a const graph.
      // The data must be passed in as nonconst, so we have to copy it
      // first.
      ArrayRCP<size_t> newPtr (ptr.size ());
      std::copy (ptr.begin (), ptr.end (), newPtr.begin ());
      ArrayRCP<LocalOrdinal> newInd (ind.size ());
      std::copy (ind.begin (), ind.end (), newInd.begin ());
      newmat->setAllValues (newPtr, newInd, newVals1D);

      // We already have the Import and Export (if applicable) objects
      // from the graph, so we can save a lot of time by passing them
      // in to expertStaticFillComplete.
      RCP<const map_type> domainMap = this->getDomainMap ();
      RCP<const map_type> rangeMap = this->getRangeMap ();
      RCP<const import_type> importer = this->getCrsGraph ()->getImporter ();
      RCP<const export_type> exporter = this->getCrsGraph ()->getExporter ();
      newmat->expertStaticFillComplete (domainMap, rangeMap, importer, exporter);
    }

    return newmat;
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  checkInternalState () const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "checkInternalState: ";
    const char err[] = "Internal state is not consistent.  "
      "Please report this bug to the Tpetra developers.";

    // This version of the graph (RCP<const crs_graph_type>) must
    // always be nonnull.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      staticGraph_.is_null (),
      std::logic_error, err);
    // myGraph == null means that the matrix has a const ("static")
    // graph.  Otherwise, the matrix has a dynamic graph (it owns its
    // graph).
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! myGraph_.is_null () && myGraph_ != staticGraph_,
      std::logic_error, err);
    // if matrix is fill complete, then graph must be fill complete
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete () && ! staticGraph_->isFillComplete (),
      std::logic_error, err << "  Specifically, the matrix is fill complete, "
      "but its graph is NOT fill complete.");
    // if matrix is storage optimized, it should have a 1D allocation
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized () && ! values2D_.is_null (),
      std::logic_error, err);
    // if matrix/graph are static profile, then 2D allocation should not be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getProfileType() == StaticProfile && values2D_ != null,
      std::logic_error, err);
    // if matrix/graph are dynamic profile, then 1D allocation should not be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getProfileType() == DynamicProfile && k_values1D_.dimension_0 () > 0,
      std::logic_error, err);
    // if values are allocated and they are non-zero in number, then
    // one of the allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      staticGraph_->indicesAreAllocated () &&
      staticGraph_->getNodeAllocationSize() > 0 &&
      staticGraph_->getNodeNumRows() > 0
      && values2D_.is_null () &&
      k_values1D_.dimension_0 () == 0,
      std::logic_error, err);
    // we cannot have both a 1D and 2D allocation
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_values1D_.dimension_0 () > 0 && values2D_ != null,
      std::logic_error, err << "  Specifically, k_values1D_ is allocated (has "
      "size " << k_values1D_.dimension_0 () << " > 0) and values2D_ is also "
      "allocated.  CrsMatrix is not suppose to have both a 1-D and a 2-D "
      "allocation at the same time.");
#endif
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  std::string
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  description () const
  {
    std::ostringstream os;

    os << "Tpetra::CrsMatrix (Kokkos refactor): {";
    if (this->getObjectLabel () != "") {
      os << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    if (isFillComplete ()) {
      os << "isFillComplete: true"
         << ", global dimensions: [" << getGlobalNumRows () << ", "
         << getGlobalNumCols () << "]"
         << ", global number of entries: " << getGlobalNumEntries ()
         << "}";
    }
    else {
      os << "isFillComplete: false"
         << ", global dimensions: [" << getGlobalNumRows () << ", "
         << getGlobalNumCols () << "]}";
    }
    return os.str ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // Don't print anything at all
    }
    // By convention, describe() always begins with a tab.
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = this->getComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, static_cast<size_t> (11)) + 2;

    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per process
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (myRank == 0) {
      out << "Tpetra::CrsMatrix (Kokkos refactor):" << endl;
    }
    Teuchos::OSTab tab1 (out);

    if (myRank == 0) {
      if (this->getObjectLabel () != "") {
        out << "Label: \"" << this->getObjectLabel () << "\", ";
      }
      {
        out << "Template parameters:" << endl;
        Teuchos::OSTab tab2 (out);
        out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
            << "LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name () << endl
            << "GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name () << endl
            << "Node: " << Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>::name () << endl;
      }
      if (isFillComplete()) {
        out << "isFillComplete: true" << endl
            << "Global dimensions: [" << getGlobalNumRows () << ", "
            << getGlobalNumCols () << "]" << endl
            << "Global number of entries: " << getGlobalNumEntries () << endl
            << "Global number of diagonal entries: " << getGlobalNumDiags ()
            << endl << "Global max number of entries in a row: "
            << getGlobalMaxNumRowEntries () << endl;
      }
      else {
        out << "isFillComplete: false" << endl
            << "Global dimensions: [" << getGlobalNumRows () << ", "
            << getGlobalNumCols () << "]" << endl;
      }
    }

    if (vl < VERB_MEDIUM) {
      return; // all done!
    }

    // Describe the Map Map.
    if (myRank == 0) {
      out << endl << "Row Map:" << endl;
    }
    getRowMap ()->describe (out, vl);

    // Describe the column Map.
    if (myRank == 0) {
      out << "Column Map: ";
    }
    if (getColMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getColMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getColMap ()->describe (out, vl);
    }

    // Describe the domain Map.
    if (myRank == 0) {
      out << "Domain Map: ";
    }
    if (getDomainMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getDomainMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else if (getDomainMap () == getColMap ()) {
      if (myRank == 0) {
        out << "same as column Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getColMap ()->describe (out, vl);
    }

    // Describe the range Map.
    if (myRank == 0) {
      out << "Range Map: ";
    }
    if (getRangeMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getRangeMap () == getDomainMap ()) {
      if (myRank == 0) {
        out << "same as domain Map" << endl;
      }
    } else if (getRangeMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getColMap ()->describe (out, vl);
    }

    // O(P) data
    for (int curRank = 0; curRank < numProcs; ++curRank) {
      if (myRank == curRank) {
        out << "Process rank: " << curRank << endl;
        Teuchos::OSTab tab2 (out);
        if (! staticGraph_->indicesAreAllocated ()) {
          out << "Graph indices not allocated" << endl;
        }
        else {
          out << "Number of allocated entries: "
              << staticGraph_->getNodeAllocationSize () << endl;
        }
        out << "Number of entries: " << getNodeNumEntries () << endl;
        if (isFillComplete ()) {
          out << "Number of diagonal entries: " << getNodeNumDiags () << endl;
        }
        out << "Max number of entries per row: " << getNodeMaxNumRowEntries ()
            << endl;
      }
      // Give output time to complete by executing some barriers.
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    }

    if (vl < VERB_HIGH) {
      return; // all done!
    }

    // O(N) and O(NNZ) data
    for (int curRank = 0; curRank < numProcs; ++curRank) {
      if (myRank == curRank) {
        out << std::setw(width) << "Proc Rank"
            << std::setw(width) << "Global Row"
            << std::setw(width) << "Num Entries";
        if (vl == VERB_EXTREME) {
          out << std::setw(width) << "(Index,Value)";
        }
        out << endl;
        for (size_t r = 0; r < getNodeNumRows (); ++r) {
          const size_t nE = getNumEntriesInLocalRow(r);
          GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
          out << std::setw(width) << myRank
              << std::setw(width) << gid
              << std::setw(width) << nE;
          if (vl == VERB_EXTREME) {
            if (isGloballyIndexed()) {
              ArrayView<const GlobalOrdinal> rowinds;
              ArrayView<const Scalar> rowvals;
              getGlobalRowView (gid, rowinds, rowvals);
              for (size_t j = 0; j < nE; ++j) {
                out << " (" << rowinds[j]
                    << ", " << rowvals[j]
                    << ") ";
              }
            }
            else if (isLocallyIndexed()) {
              ArrayView<const LocalOrdinal> rowinds;
              ArrayView<const Scalar> rowvals;
              getLocalRowView (r, rowinds, rowvals);
              for (size_t j=0; j < nE; ++j) {
                out << " (" << getColMap()->getGlobalElement(rowinds[j])
                    << ", " << rowvals[j]
                    << ") ";
              }
            } // globally or locally indexed
          } // vl == VERB_EXTREME
          out << endl;
        } // for each row r on this process
      } // if (myRank == curRank)

      // Give output time to complete
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    } // for each process p
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  checkSizes (const SrcDistObject& source)
  {
    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.

    // Currently, the source object must be a RowMatrix with the same
    // four template parameters as the target CrsMatrix.  We might
    // relax this requirement later.
    typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> row_matrix_type;
    const row_matrix_type* srcRowMat =
      dynamic_cast<const row_matrix_type*> (&source);
    return (srcRowMat != NULL);
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  copyAndPermute (const SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                  const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef node_type NT;
    // Method name string for TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC.
    const char tfecfFuncName[] = "copyAndPermute: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(),
      std::invalid_argument, "permuteToLIDs.size() = " << permuteToLIDs.size()
      << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");

    // This dynamic cast should succeed, because we've already tested
    // it in checkSizes().
    typedef RowMatrix<Scalar, LO, GO, NT> row_matrix_type;
    const row_matrix_type& srcMat = dynamic_cast<const row_matrix_type&> (source);

    const bool sourceIsLocallyIndexed = srcMat.isLocallyIndexed ();
    //
    // Copy the first numSame row from source to target (this matrix).
    // This involves copying rows corresponding to LIDs [0, numSame-1].
    //
    const map_type& srcRowMap = * (srcMat.getRowMap ());
    Array<GO> rowInds;
    Array<Scalar> rowVals;
    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);
    for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; ++sourceLID) {
      // Global ID for the current row index in the source matrix.
      // The first numSameIDs GIDs in the two input lists are the
      // same, so sourceGID == targetGID in this case.
      const GO sourceGID = srcRowMap.getGlobalElement (sourceLID);
      const GO targetGID = sourceGID;

      // Input views for the combineGlobalValues() call below.
      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.size())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy (sourceGID, rowIndsView, rowValsView, checkRowLength);

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowLength != checkRowLength,
          std::logic_error, "For global row index " << sourceGID << ", the source"
          " matrix's getNumEntriesInGlobalRow() method returns a row length of "
          << rowLength << ", but the getGlobalRowCopy() method reports that "
          "the row length is " << checkRowLength << ".  Please report this bug "
          "to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else { // source matrix is globally indexed.
        srcMat.getGlobalRowView (sourceGID, rowIndsConstView, rowValsConstView);
      }

      // Combine the data into the target matrix.
      if (isStaticGraph()) {
        // Applying a permutation to a matrix with a static graph
        // means REPLACE-ing entries.
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, REPLACE);
      }
      else {
        // Applying a permutation to a matrix with a dynamic graph
        // means INSERT-ing entries.  This has the same effect as
        // ADD, if the target graph already has an entry there.
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, INSERT);
      }
    } // For each of the consecutive source and target IDs that are the same

    //
    // Permute the remaining rows.
    //
    const map_type& tgtRowMap = * (this->getRowMap ());
    const size_t numPermuteToLIDs = static_cast<size_t> (permuteToLIDs.size ());
    for (size_t p = 0; p < numPermuteToLIDs; ++p) {
      const GO sourceGID = srcRowMap.getGlobalElement (permuteFromLIDs[p]);
      const GO targetGID = tgtRowMap.getGlobalElement (permuteToLIDs[p]);

      // Input views for the combineGlobalValues() call below.
      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.size ())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy (sourceGID, rowIndsView, rowValsView, checkRowLength);

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowLength != checkRowLength,
          std::logic_error, "For the source matrix's global row index "
          << sourceGID << ", the source matrix's getNumEntriesInGlobalRow() "
          "method returns a row length of " << rowLength << ", but the "
          "getGlobalRowCopy() method reports that the row length is "
          << checkRowLength << ".  Please report this bug to the Tpetra "
          "developers.");
#endif // HAVE_TPETRA_DEBUG

        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else {
        srcMat.getGlobalRowView (sourceGID, rowIndsConstView, rowValsConstView);
      }

      // Combine the data into the target matrix.
      if (isStaticGraph()) {
        this->combineGlobalValues (targetGID, rowIndsConstView,
                                   rowValsConstView, REPLACE);
      }
      else {
        this->combineGlobalValues (targetGID, rowIndsConstView,
                                   rowValsConstView, INSERT);
      }
    } // For each ID to permute
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  packAndPrepare (const SrcDistObject& source,
                  const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                  Teuchos::Array<char>& exports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor& distor)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "packAndPrepare: ";

    // Attempt to cast the source object to RowMatrix.  If the cast
    // succeeds, use the source object's pack method to pack its data
    // for communication.  If the source object is really a CrsMatrix,
    // this will pick up the CrsMatrix's more efficient override.  If
    // the RowMatrix cast fails, then the source object doesn't have
    // the right type.
    //
    // FIXME (mfh 30 Jun 2013) We don't even need the RowMatrix to
    // have the same Node type.  Unfortunately, we don't have a way to
    // ask if the RowMatrix is "a RowMatrix with any Node type," since
    // RowMatrix doesn't have a base class.  A hypothetical
    // RowMatrixBase<Scalar, LO, GO> class, which does not currently
    // exist, would satisfy this requirement.
    //
    // Why RowMatrixBase<Scalar, LO, GO>?  The source object's Scalar
    // type doesn't technically need to match the target object's
    // Scalar type, so we could just have RowMatrixBase<LO, GO>.  LO
    // and GO need not be the same, as long as there is no overflow of
    // the indices.  However, checking for index overflow is global
    // and therefore undesirable.
    typedef RowMatrix<Scalar, LO, GO, node_type> row_matrix_type;
    const row_matrix_type* srcRowMat =
      dynamic_cast<const row_matrix_type*> (&source);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      srcRowMat == NULL, std::invalid_argument,
      "The source object of the Import or Export operation is neither a "
      "CrsMatrix (with the same template parameters as the target object), "
      "nor a RowMatrix (with the same first four template parameters as the "
      "target object).");
#ifdef HAVE_TPETRA_DEBUG
    {
      using Teuchos::reduceAll;
      std::ostringstream msg;
      int lclBad = 0;
      try {
        srcRowMat->pack (exportLIDs, exports, numPacketsPerLID,
                         constantNumPackets, distor);
      } catch (std::exception& e) {
        lclBad = 1;
        msg << e.what ();
      }
      int gblBad = 0;
      const Teuchos::Comm<int>& comm = * (this->getComm ());
      reduceAll<int, int> (comm, Teuchos::REDUCE_MAX,
                           lclBad, Teuchos::outArg (gblBad));
      if (gblBad != 0) {
        const int myRank = comm.getRank ();
        const int numProcs = comm.getSize ();
        for (int r = 0; r < numProcs; ++r) {
          if (r == myRank && lclBad != 0) {
            std::ostringstream os;
            os << "Proc " << myRank << ": " << msg.str () << std::endl;
            std::cerr << os.str ();
          }
          comm.barrier ();
          comm.barrier ();
          comm.barrier ();
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "pack() threw an exception on one or "
          "more participating processes.");
      }
    }
#else
    srcRowMat->pack (exportLIDs, exports, numPacketsPerLID,
                     constantNumPackets, distor);
#endif // HAVE_TPETRA_DEBUG
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  packRow (char* const numEntOut,
           char* const valOut,
           char* const indOut,
           const size_t numEnt,
           const LocalOrdinal lclRow) const
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const LO numEntLO = static_cast<LO> (numEnt);
    memcpy (numEntOut, &numEntLO, sizeof (LO));
    if (this->isLocallyIndexed ()) {
      // If the matrix is locally indexed on the calling process, we
      // have to use its column Map (which it _must_ have in this
      // case) to convert to global indices.
      ArrayView<const LO> indIn;
      ArrayView<const Scalar> valIn;
      this->getLocalRowView (lclRow, indIn, valIn);
      const map_type& colMap = * (this->getColMap ());
      // Copy column indices one at a time, so that we don't need
      // temporary storage.
      for (size_t k = 0; k < numEnt; ++k) {
        const GO gblIndIn = colMap.getGlobalElement (indIn[k]);
        memcpy (indOut + k * sizeof (GO), &gblIndIn, sizeof (GO));
      }
      memcpy (valOut, valIn.getRawPtr (), numEnt * sizeof (Scalar));
    }
    else if (this->isGloballyIndexed ()) {
      // If the matrix is globally indexed on the calling process,
      // then we can use the column indices directly.  However, we
      // have to get the global row index.  The calling process must
      // have a row Map, since otherwise it shouldn't be participating
      // in packing operations.
      ArrayView<const GO> indIn;
      ArrayView<const Scalar> valIn;
      const map_type& rowMap = * (this->getRowMap ());
      const GO gblRow = rowMap.getGlobalElement (lclRow);
      this->getGlobalRowView (gblRow, indIn, valIn);
      memcpy (indOut, indIn.getRawPtr (), numEnt * sizeof (GO));
      memcpy (valOut, valIn.getRawPtr (), numEnt * sizeof (Scalar));
    }
    else {
      if (numEnt != 0) {
        return false;
      }
    }
    return true;
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  bool
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  unpackRow (Scalar* const valInTmp,
             GlobalOrdinal* const indInTmp,
             const size_t tmpSize,
             const char* const valIn,
             const char* const indIn,
             const size_t numEnt,
             const LocalOrdinal lclRow,
             const Tpetra::CombineMode combineMode)
  {
    if (tmpSize < numEnt || (numEnt != 0 && (valInTmp == NULL || indInTmp == NULL))) {
      return false;
    }
    memcpy (valInTmp, valIn, numEnt * sizeof (Scalar));
    memcpy (indInTmp, indIn, numEnt * sizeof (GlobalOrdinal));
    const GlobalOrdinal gblRow = this->getRowMap ()->getGlobalElement (lclRow);
    Teuchos::ArrayView<Scalar> val ((numEnt == 0) ? NULL : valInTmp, numEnt);
    Teuchos::ArrayView<GlobalOrdinal> ind ((numEnt == 0) ? NULL : indInTmp, numEnt);
    this->combineGlobalValues (gblRow, ind, val, combineMode);
    return true;
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  allocatePackSpace (Teuchos::Array<char>& exports,
                     size_t& totalNumEntries,
                     const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
    //const char tfecfFuncName[] = "allocatePackSpace: ";
    const size_type numExportLIDs = exportLIDs.size ();

    // Count the total number of entries to send.
    totalNumEntries = 0;
    for (size_type i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs[i];
      size_t curNumEntries = this->getNumEntriesInLocalRow (lclRow);
      // FIXME (mfh 25 Jan 2015) We should actually report invalid row
      // indices as an error.  Just consider them nonowned for now.
      if (curNumEntries == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        curNumEntries = 0;
      }
      totalNumEntries += curNumEntries;
    }

    // FIXME (mfh 24 Feb 2013) This code is only correct if
    // sizeof(Scalar) is a meaningful representation of the amount of
    // data in a Scalar instance.  (LO and GO are always built-in
    // integer types.)
    //
    // Allocate the exports array.  It does NOT need padding for
    // alignment, since we use memcpy to write to / read from send /
    // receive buffers.
    const size_t allocSize =
      static_cast<size_t> (numExportLIDs) * sizeof (LO) +
      totalNumEntries * (sizeof (Scalar) + sizeof (GO));
    if (static_cast<size_t> (exports.size ()) < allocSize) {
      exports.resize (allocSize);
    }
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<char>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Distributor& distor) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::RCP;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const LO>::size_type size_type;
    const char tfecfFuncName[] = "pack: ";

    const size_type numExportLIDs = exportLIDs.size ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numExportLIDs != numPacketsPerLID.size (), std::invalid_argument,
      "exportLIDs.size() = " << numExportLIDs << " != numPacketsPerLID.size()"
      " = " << numPacketsPerLID.size () << ".");

    // Setting this to zero tells the caller to expect a possibly
    // different ("nonconstant") number of packets per local index
    // (i.e., a possibly different number of entries per row).
    constantNumPackets = 0;

    // The pack buffer 'exports' enters this method possibly
    // unallocated.  Do the first two parts of "Count, allocate, fill,
    // compute."
    size_t totalNumEntries = 0;
    allocatePackSpace (exports, totalNumEntries, exportLIDs);
    const size_t bufSize = static_cast<size_t> (exports.size ());

    // Compute the number of "packets" (in this case, bytes) per
    // export LID (in this case, local index of the row to send), and
    // actually pack the data.
    //
    // FIXME (mfh 24 Feb 2013, 25 Jan 2015) This code is only correct
    // if sizeof(Scalar) is a meaningful representation of the amount
    // of data in a Scalar instance.  (LO and GO are always built-in
    // integer types.)

    // Variables for error reporting in the loop.
    size_type firstBadIndex = 0; // only valid if outOfBounds == true.
    size_t firstBadOffset = 0;   // only valid if outOfBounds == true.
    size_t firstBadNumBytes = 0; // only valid if outOfBounds == true.
    bool outOfBounds = false;
    bool packErr = false;

    char* const exportsRawPtr = exports.getRawPtr ();
    size_t offset = 0; // current index into 'exports' array.
    for (size_type i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs[i];
      const size_t numEnt = this->getNumEntriesInLocalRow (lclRow);

      // Only pad this row if it has a nonzero number of entries.
      if (numEnt == 0) {
        numPacketsPerLID[i] = 0;
      }
      else {
        char* const numEntBeg = exportsRawPtr + offset;
        char* const numEntEnd = numEntBeg + sizeof (LO);
        char* const valBeg = numEntEnd;
        char* const valEnd = valBeg + numEnt * sizeof (Scalar);
        char* const indBeg = valEnd;
        const size_t numBytes = sizeof (LO) +
          numEnt * (sizeof (Scalar) + sizeof (GO));
        if (offset > bufSize || offset + numBytes > bufSize) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadNumBytes = numBytes;
          outOfBounds = true;
          break;
        }
        packErr = ! packRow (numEntBeg, valBeg, indBeg, numEnt, lclRow);
        if (packErr) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadNumBytes = numBytes;
          break;
        }
        // numPacketsPerLID[i] is the number of "packets" in the
        // current local row i.  Packet=char (really "byte") so use
        // the number of bytes of the packed data for that row.
        numPacketsPerLID[i] = numBytes;
        offset += numBytes;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      outOfBounds, std::logic_error, "First invalid offset into 'exports' "
      "pack buffer at index i = " << firstBadIndex << ".  exportLIDs[i]: "
      << exportLIDs[firstBadIndex] << ", bufSize: " << bufSize << ", offset: "
      << firstBadOffset << ", numBytes: " << firstBadNumBytes << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      packErr, std::logic_error, "First error in packRow() at index i = "
      << firstBadIndex << ".  exportLIDs[i]: " << exportLIDs[firstBadIndex]
      << ", bufSize: " << bufSize << ", offset: " << firstBadOffset
      << ", numBytes: " << firstBadNumBytes << ".");
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  combineGlobalValues (const GlobalOrdinal globalRowIndex,
                       const Teuchos::ArrayView<const GlobalOrdinal>& columnIndices,
                       const Teuchos::ArrayView<const Scalar>& values,
                       const Tpetra::CombineMode combineMode)
  {
    const char tfecfFuncName[] = "combineGlobalValues: ";

    if (isStaticGraph ()) {
      // INSERT doesn't make sense for a static graph, since you
      // aren't allowed to change the structure of the graph.
      // However, all the other combine modes work.
      if (combineMode == ADD) {
        sumIntoGlobalValues (globalRowIndex, columnIndices, values);
      }
      else if (combineMode == REPLACE) {
        replaceGlobalValues (globalRowIndex, columnIndices, values);
      }
      else if (combineMode == ABSMAX) {
        using Details::AbsMax;
        AbsMax<Scalar> f;
        this->template transformGlobalValues<AbsMax<Scalar> > (globalRowIndex,
                                                               columnIndices,
                                                               values, f);
      }
      else if (combineMode == INSERT) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          isStaticGraph () && combineMode == INSERT, std::invalid_argument,
          "INSERT combine mode is not allowed if the matrix has a static graph "
          "(i.e., was constructed with the CrsMatrix constructor that takes a "
          "const CrsGraph pointer).");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "Invalid combine mode; should never get "
          "here!  Please report this bug to the Tpetra developers.");
      }
    }
    else { // The matrix has a dynamic graph.
      if (combineMode == ADD || combineMode == INSERT) {
        // For a dynamic graph, all incoming column indices are
        // inserted into the target graph.  Duplicate indices will
        // have their values summed.  In this context, ADD and INSERT
        // are equivalent.  We need to call insertGlobalValues()
        // anyway if the column indices don't yet exist in this row,
        // so we just call insertGlobalValues() for both cases.
        insertGlobalValuesFiltered (globalRowIndex, columnIndices, values);
      }
      // FIXME (mfh 14 Mar 2012):
      //
      // Implementing ABSMAX or REPLACE for a dynamic graph would
      // require modifying assembly to attach a possibly different
      // combine mode to each inserted (i, j, A_ij) entry.  For
      // example, consider two different Export operations to the same
      // target CrsMatrix, the first with ABSMAX combine mode and the
      // second with REPLACE.  This isn't a common use case, so we
      // won't mess with it for now.
      else if (combineMode == ABSMAX) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          ! isStaticGraph () && combineMode == ABSMAX, std::logic_error,
          "ABSMAX combine mode when the matrix has a dynamic graph is not yet "
          "implemented.");
      }
      else if (combineMode == REPLACE) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          ! isStaticGraph () && combineMode == REPLACE, std::logic_error,
          "REPLACE combine mode when the matrix has a dynamic graph is not yet "
          "implemented.");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "Should never get here!  Please report this "
          "bug to the Tpetra developers.");
      }
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
                    const Teuchos::ArrayView<const char>& imports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor& distor,
                    CombineMode combineMode)
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "unpackAndCombine: ";
    const CombineMode validModes[4] = {ADD, REPLACE, ABSMAX, INSERT};
    const char* validModeNames[4] = {"ADD", "REPLACE", "ABSMAX", "INSERT"};
    const int numValidModes = 4;

    if (std::find (validModes, validModes+numValidModes, combineMode) ==
        validModes+numValidModes) {
      std::ostringstream os;
      os << "Invalid combine mode.  Valid modes are {";
      for (int k = 0; k < numValidModes; ++k) {
        os << validModeNames[k];
        if (k < numValidModes - 1) {
          os << ", ";
        }
      }
      os << "}.";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::invalid_argument, os.str ());
    }

    {
      using Teuchos::reduceAll;
      std::ostringstream msg;
      int lclBad = 0;
      try {
        this->unpackAndCombineImpl (importLIDs, imports, numPacketsPerLID,
                                    constantNumPackets, distor, combineMode);
      } catch (std::exception& e) {
        lclBad = 1;
        msg << e.what ();
      }
      int gblBad = 0;
      const Teuchos::Comm<int>& comm = * (this->getComm ());
      reduceAll<int, int> (comm, Teuchos::REDUCE_MAX,
                           lclBad, Teuchos::outArg (gblBad));
      if (gblBad != 0) {
        const int myRank = comm.getRank ();
        const int numProcs = comm.getSize ();
        for (int r = 0; r < numProcs; ++r) {
          if (r == myRank && lclBad != 0) {
            std::ostringstream os;
            os << "Proc " << myRank << ": " << msg.str () << std::endl;
            std::cerr << os.str ();
          }
          comm.barrier ();
          comm.barrier ();
          comm.barrier ();
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "unpackAndCombineImpl() threw an "
          "exception on one or more participating processes.");
      }
    }
#else
    this->unpackAndCombineImpl (importLIDs, imports, numPacketsPerLID,
                                constantNumPackets, distor, combineMode);
#endif // HAVE_TPETRA_DEBUG
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  unpackAndCombineImpl (const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
                        const Teuchos::ArrayView<const char>& imports,
                        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                        size_t constantNumPackets,
                        Distributor & /* distor */,
                        CombineMode combineMode)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
    const char tfecfFuncName[] = "unpackAndCombine: ";

#ifdef HAVE_TPETRA_DEBUG
    const CombineMode validModes[4] = {ADD, REPLACE, ABSMAX, INSERT};
    const char* validModeNames[4] = {"ADD", "REPLACE", "ABSMAX", "INSERT"};
    const int numValidModes = 4;

    if (std::find (validModes, validModes+numValidModes, combineMode) ==
        validModes+numValidModes) {
      std::ostringstream os;
      os << "Invalid combine mode.  Valid modes are {";
      for (int k = 0; k < numValidModes; ++k) {
        os << validModeNames[k];
        if (k < numValidModes - 1) {
          os << ", ";
        }
      }
      os << "}.";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::invalid_argument, os.str ());
    }
#endif // HAVE_TPETRA_DEBUG

    const size_type numImportLIDs = importLIDs.size ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numImportLIDs != numPacketsPerLID.size (), std::invalid_argument,
      "importLIDs.size() = " << numImportLIDs << "  != numPacketsPerLID.size()"
      << " = " << numPacketsPerLID.size () << ".");

    // If a sanity check fails, keep track of some state at the
    // "first" place where it fails.  After the first failure, "run
    // through the motions" until the end of this method, then raise
    // an error with an informative message.
    size_type firstBadIndex = 0;
    size_t firstBadOffset = 0;
    size_t firstBadExpectedNumBytes = 0;
    size_t firstBadNumBytes = 0;
    LO firstBadNumEnt = 0;
    // We have sanity checks for three kinds of errors:
    //
    //   1. Offset into array of all the incoming data (for all rows)
    //      is out of bounds
    //   2. Too few bytes of incoming data for a row, given the
    //      reported number of entries in those incoming data
    //   3. Error in unpacking the row's incoming data
    //
    bool outOfBounds = false;
    bool wrongNumBytes = false;
    bool unpackErr = false;

    const size_t bufSize = static_cast<size_t> (imports.size ());
    const char* const importsRawPtr = imports.getRawPtr ();
    size_t offset = 0;

    // Temporary storage for incoming values and indices.  We need
    // this because the receive buffer does not align storage; it's
    // just contiguous bytes.  In order to avoid violating ANSI
    // aliasing rules, we memcpy each incoming row's data into these
    // temporary arrays.  We double their size every time we run out
    // of storage.
    Array<Scalar> valInTmp;
    Array<GO> indInTmp;
    for (size_type i = 0; i < numImportLIDs; ++i) {
      const LO lclRow = importLIDs[i];
      const size_t numBytes = numPacketsPerLID[i];

      if (numBytes > 0) { // there is actually something in the row
        const char* const numEntBeg = importsRawPtr + offset;
        const char* const numEntEnd = numEntBeg + sizeof (LO);

        // Now we know how many entries to expect in the received data
        // for this row.
        LO numEnt = 0;
        memcpy (&numEnt, numEntBeg, sizeof (LO));

        const char* const valBeg = numEntEnd;
        const char* const valEnd =
          valBeg + static_cast<size_t> (numEnt) * sizeof (Scalar);
        const char* const indBeg = valEnd;
        const size_t expectedNumBytes = sizeof (LO) +
          static_cast<size_t> (numEnt) * (sizeof (Scalar) + sizeof (GO));

        if (expectedNumBytes > numBytes) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadExpectedNumBytes = expectedNumBytes;
          firstBadNumBytes = numBytes;
          firstBadNumEnt = numEnt;
          wrongNumBytes = true;
          break;
        }
        if (offset > bufSize || offset + numBytes > bufSize) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadExpectedNumBytes = expectedNumBytes;
          firstBadNumBytes = numBytes;
          firstBadNumEnt = numEnt;
          outOfBounds = true;
          break;
        }
        size_t tmpNumEnt = static_cast<size_t> (valInTmp.size ());
        if (tmpNumEnt < static_cast<size_t> (numEnt) ||
            static_cast<size_t> (indInTmp.size ()) < static_cast<size_t> (numEnt)) {
          // Double the size of the temporary arrays for incoming data.
          tmpNumEnt = std::max (static_cast<size_t> (numEnt), tmpNumEnt * 2);
          valInTmp.resize (tmpNumEnt);
          indInTmp.resize (tmpNumEnt);
        }
        unpackErr =
          ! unpackRow (valInTmp.getRawPtr (), indInTmp.getRawPtr (), tmpNumEnt,
                       valBeg, indBeg, numEnt, lclRow, combineMode);
        if (unpackErr) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadExpectedNumBytes = expectedNumBytes;
          firstBadNumBytes = numBytes;
          firstBadNumEnt = numEnt;
          break;
        }
        offset += numBytes;
      }
    }

    if (wrongNumBytes || outOfBounds || unpackErr) {
      std::ostringstream os;
      os << "  importLIDs[i]: " << importLIDs[firstBadIndex]
         << ", bufSize: " << bufSize
         << ", offset: " << firstBadOffset
         << ", numBytes: " << firstBadNumBytes
         << ", expectedNumBytes: " << firstBadExpectedNumBytes
         << ", numEnt: " << firstBadNumEnt;
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        wrongNumBytes, std::logic_error, "At index i = " << firstBadIndex
        << ", expectedNumBytes > numBytes." << os.str ());
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        outOfBounds, std::logic_error, "First invalid offset into 'imports' "
        "unpack buffer at index i = " << firstBadIndex << "." << os.str ());
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        unpackErr, std::logic_error, "First error in unpackRow() at index i = "
        << firstBadIndex << "." << os.str ());
    }
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
                           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getColumnMapMultiVector (const MV& X_domainMap,
                           const bool force) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->hasColMap (), std::runtime_error, "Tpetra::CrsMatrix::getColumn"
      "MapMultiVector: You may only call this method if the matrix has a "
      "column Map.  If the matrix does not yet have a column Map, you should "
      "first call fillComplete (with domain and range Map if necessary).");

    // If the graph is not fill complete, then the Import object (if
    // one should exist) hasn't been constructed yet.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::"
      "CrsMatrix::getColumnMapMultiVector: You may only call this method if "
      "this matrix's graph is fill complete.");

    const size_t numVecs = X_domainMap.getNumVectors ();
    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const map_type> colMap = this->getColMap ();

    RCP<MV> X_colMap; // null by default

    // If the Import object is trivial (null), then we don't need a
    // separate column Map multivector.  Just return null in that
    // case.  The caller is responsible for knowing not to use the
    // returned null pointer.
    //
    // If the Import is nontrivial, then we do need a separate
    // column Map multivector for the Import operation.  Check in
    // that case if we have to (re)create the column Map
    // multivector.
    if (! importer.is_null () || force) {
      if (importMV_.is_null () || importMV_->getNumVectors () != numVecs) {
        X_colMap = rcp (new MV (colMap, numVecs));

        // Cache the newly created multivector for later reuse.
        importMV_ = X_colMap;
      }
      else { // Yay, we can reuse the cached multivector!
        X_colMap = importMV_;
        // mfh 09 Jan 2013: We don't have to fill with zeros first,
        // because the Import uses INSERT combine mode, which overwrites
        // existing entries.
        //
        //X_colMap->putScalar (ZERO);
      }
    }
    return X_colMap;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
                           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getRowMapMultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& Y_rangeMap,
                        const bool force) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // If the graph is not fill complete, then the Export object (if
    // one should exist) hasn't been constructed yet.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::"
      "CrsMatrix::getRowMapMultiVector: You may only call this method if this "
      "matrix's graph is fill complete.");

    const size_t numVecs = Y_rangeMap.getNumVectors ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    // Every version of the constructor takes either a row Map, or a
    // graph (all of whose constructors take a row Map).  Thus, the
    // matrix always has a row Map.
    RCP<const map_type> rowMap = this->getRowMap ();

    RCP<MV> Y_rowMap; // null by default

    // If the Export object is trivial (null), then we don't need a
    // separate row Map multivector.  Just return null in that case.
    // The caller is responsible for knowing not to use the returned
    // null pointer.
    //
    // If the Export is nontrivial, then we do need a separate row
    // Map multivector for the Export operation.  Check in that case
    // if we have to (re)create the row Map multivector.
    if (! exporter.is_null () || force) {
      if (exportMV_.is_null () || exportMV_->getNumVectors () != numVecs) {
        Y_rowMap = rcp (new MV (rowMap, numVecs));
        exportMV_ = Y_rowMap; // Cache the newly created MV for later reuse.
      }
      else { // Yay, we can reuse the cached multivector!
        Y_rowMap = exportMV_;
      }
    }
    return Y_rowMap;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      myGraph_.is_null (), std::logic_error, "Tpetra::CrsMatrix::"
      "removeEmptyProcessesInPlace: This method does not work when the matrix "
      "was created with a constant graph (that is, when it was created using "
      "the version of its constructor that takes an RCP<const CrsGraph>).  "
      "This is because the matrix is not allowed to modify the graph in that "
      "case, but removing empty processes requires modifying the graph.");
    myGraph_->removeEmptyProcessesInPlace (newMap);
    // Even though CrsMatrix's row Map (as returned by getRowMap())
    // comes from its CrsGraph, CrsMatrix still implements DistObject,
    // so we also have to change the DistObject's Map.
    this->map_ = this->getRowMap ();
    // In the nonconst graph case, staticGraph_ is just a const
    // pointer to myGraph_.  This assignment is probably redundant,
    // but it doesn't hurt.
    staticGraph_ = Teuchos::rcp_const_cast<const Graph> (myGraph_);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  add (const Scalar& alpha,
       const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
       const Scalar& beta,
       const Teuchos::RCP<const map_type>& domainMap,
       const Teuchos::RCP<const map_type>& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    using Teuchos::sublist;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> row_matrix_type;
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> crs_matrix_type;

    const crs_matrix_type& B = *this; // a convenient abbreviation
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();
    const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one ();

    // If the user didn't supply a domain or range Map, then try to
    // get one from B first (if it has them), then from A (if it has
    // them).  If we don't have any domain or range Maps, scold the
    // user.
    RCP<const map_type> A_domainMap = A.getDomainMap ();
    RCP<const map_type> A_rangeMap = A.getRangeMap ();
    RCP<const map_type> B_domainMap = B.getDomainMap ();
    RCP<const map_type> B_rangeMap = B.getRangeMap ();

    RCP<const map_type> theDomainMap = domainMap;
    RCP<const map_type> theRangeMap = rangeMap;

    if (domainMap.is_null ()) {
      if (B_domainMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_domainMap.is_null (), std::invalid_argument,
          "Tpetra::CrsMatrix::add: If neither A nor B have a domain Map, "
          "then you must supply a nonnull domain Map to this method.");
        theDomainMap = A_domainMap;
      } else {
        theDomainMap = B_domainMap;
      }
    }
    if (rangeMap.is_null ()) {
      if (B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_rangeMap.is_null (), std::invalid_argument,
          "Tpetra::CrsMatrix::add: If neither A nor B have a range Map, "
          "then you must supply a nonnull range Map to this method.");
        theRangeMap = A_rangeMap;
      } else {
        theRangeMap = B_rangeMap;
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, check that A and B have matching domain and
    // range Maps, if they have domain and range Maps at all.  (If
    // they aren't fill complete, then they may not yet have them.)
    if (! A_domainMap.is_null () && ! A_rangeMap.is_null ()) {
      if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_domainMap->isSameAs (*A_domainMap), std::invalid_argument,
          "Tpetra::CrsMatrix::add: The input RowMatrix A must have a domain Map "
          "which is the same as (isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_rangeMap->isSameAs (*A_rangeMap), std::invalid_argument,
          "Tpetra::CrsMatrix::add: The input RowMatrix A must have a range Map "
          "which is the same as (isSameAs) this RowMatrix's range Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
          std::invalid_argument,
          "Tpetra::CrsMatrix::add: The input domain Map must be the same as "
          "(isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
          std::invalid_argument,
          "Tpetra::CrsMatrix::add: The input range Map must be the same as "
          "(isSameAs) this RowMatrix's range Map.");
      }
    }
    else if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
        std::invalid_argument,
        "Tpetra::CrsMatrix::add: The input domain Map must be the same as "
        "(isSameAs) this RowMatrix's domain Map.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
        std::invalid_argument,
        "Tpetra::CrsMatrix::add: The input range Map must be the same as "
        "(isSameAs) this RowMatrix's range Map.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        domainMap.is_null () || rangeMap.is_null (), std::invalid_argument,
        "Tpetra::CrsMatrix::add: If neither A nor B have a domain and range "
        "Map, then you must supply a nonnull domain and range Map to this "
        "method.");
    }
#endif // HAVE_TPETRA_DEBUG

    // What parameters do we pass to C's constructor?  Do we call
    // fillComplete on C after filling it?  And if so, what parameters
    // do we pass to C's fillComplete call?
    bool callFillComplete = true;
    RCP<ParameterList> constructorSublist;
    RCP<ParameterList> fillCompleteSublist;
    if (! params.is_null ()) {
      callFillComplete = params->get ("Call fillComplete", callFillComplete);
      constructorSublist = sublist (params, "Constructor parameters");
      fillCompleteSublist = sublist (params, "fillComplete parameters");
    }

    RCP<const map_type> A_rowMap = A.getRowMap ();
    RCP<const map_type> B_rowMap = B.getRowMap ();
    RCP<const map_type> C_rowMap = B_rowMap; // see discussion in documentation
    RCP<crs_matrix_type> C; // The result matrix.

    // If A and B's row Maps are the same, we can compute an upper
    // bound on the number of entries in each row of C, before
    // actually computing the sum.  A reasonable upper bound is the
    // sum of the two entry counts in each row.  If we choose this as
    // the actual per-row upper bound, we can use static profile.
    if (A_rowMap->isSameAs (*B_rowMap)) {
      const LO localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
      ArrayRCP<size_t> C_maxNumEntriesPerRow (localNumRows, 0);

      // Get the number of entries in each row of A.
      if (alpha != ZERO) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += A_numEntries;
        }
      }
      // Get the number of entries in each row of B.
      if (beta != ZERO) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += B_numEntries;
        }
      }
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile, constructorSublist));
      }
      // Since A and B have the same row Maps, we could add them
      // together all at once and merge values before we call
      // insertGlobalValues.  However, we don't really need to, since
      // we've already allocated enough space in each row of C for C
      // to do the merge itself.
    }
    else { // the row Maps of A and B are not the same
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile,
                                      constructorSublist));
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
      "Tpetra::RowMatrix::add: C should not be null at this point.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    //
    // Compute C = alpha*A + beta*B.
    //
    Array<GO> ind;
    Array<Scalar> val;

    if (alpha != ZERO) {
      const LO A_localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (A_numEntries);
          val.resize (A_numEntries);
        }
        ArrayView<GO> indView = ind (0, A_numEntries);
        ArrayView<Scalar> valView = val (0, A_numEntries);
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != ONE) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (beta != ZERO) {
      const LO B_localNumRows = static_cast<LO> (B_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (B_numEntries);
          val.resize (B_numEntries);
        }
        ArrayView<GO> indView = ind (0, B_numEntries);
        ArrayView<Scalar> valView = val (0, B_numEntries);
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != ONE) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (callFillComplete) {
      if (fillCompleteSublist.is_null ()) {
        C->fillComplete (theDomainMap, theRangeMap);
      } else {
        C->fillComplete (theDomainMap, theRangeMap, fillCompleteSublist);
      }
    }
    return rcp_implicit_cast<row_matrix_type> (C);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  transferAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> > & destMat,
                           const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, node_type>& rowTransfer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef node_type NT;
    typedef CrsMatrix<Scalar, LO, GO, NT> this_type;
    typedef Vector<int, LO, GO, NT> IntVectorType;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string label;
    if(!params.is_null())
      label = params->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Pack-1"))));
#endif

    // Make sure that the input argument rowTransfer is either an
    // Import or an Export.  Import and Export are the only two
    // subclasses of Transfer that we defined, but users might
    // (unwisely, for now at least) decide to implement their own
    // subclasses.  Exclude this possibility.
    const import_type* xferAsImport = dynamic_cast<const import_type*> (&rowTransfer);
    const export_type* xferAsExport = dynamic_cast<const export_type*> (&rowTransfer);
    TEUCHOS_TEST_FOR_EXCEPTION(
      xferAsImport == NULL && xferAsExport == NULL, std::invalid_argument,
      "Tpetra::CrsMatrix::transferAndFillComplete: The 'rowTransfer' input "
      "argument must be either an Import or an Export, and its template "
      "parameters must match the corresponding template parameters of the "
      "CrsMatrix.");

    // FIXME (mfh 15 May 2014) Wouldn't communication still be needed,
    // if the source Map is not distributed but the target Map is?
    const bool communication_needed = rowTransfer.getSourceMap ()->isDistributed ();

    //
    // Get the caller's parameters
    //

    bool reverseMode = false; // Are we in reverse mode?
    bool restrictComm = false; // Do we need to restrict the communicator?
    RCP<ParameterList> matrixparams; // parameters for the destination matrix
    if (! params.is_null ()) {
      reverseMode = params->get ("Reverse Mode", reverseMode);
      restrictComm = params->get ("Restrict Communicator", restrictComm);
      matrixparams = sublist (params, "CrsMatrix");
    }

    // Get the new domain and range Maps.  We need some of them for
    // error checking, now that we have the reverseMode parameter.
    RCP<const map_type> MyRowMap = reverseMode ?
      rowTransfer.getSourceMap () : rowTransfer.getTargetMap ();
    RCP<const map_type> MyColMap; // create this below
    RCP<const map_type> MyDomainMap = ! domainMap.is_null () ?
      domainMap : getDomainMap ();
    RCP<const map_type> MyRangeMap = ! rangeMap.is_null () ?
      rangeMap : getRangeMap ();
    RCP<const map_type> BaseRowMap = MyRowMap;
    RCP<const map_type> BaseDomainMap = MyDomainMap;

    // If the user gave us a nonnull destMat, then check whether it's
    // "pristine."  That means that it has no entries.
    //
    // FIXME (mfh 15 May 2014) If this is not true on all processes,
    // then this exception test may hang.  It would be better to
    // forward an error flag to the next communication phase.
    if (! destMat.is_null ()) {
      // FIXME (mfh 15 May 2014): The classic Petra idiom for checking
      // whether a graph or matrix has no entries on the calling
      // process, is that it is neither locally nor globally indexed.
      // This may change eventually with the Kokkos refactor version
      // of Tpetra, so it would be better just to check the quantity
      // of interest directly.  Note that with the Kokkos refactor
      // version of Tpetra, asking for the total number of entries in
      // a graph or matrix that is not fill complete might require
      // computation (kernel launch), since it is not thread scalable
      // to update a count every time an entry is inserted.
      const bool NewFlag = ! destMat->getGraph ()->isLocallyIndexed () &&
        ! destMat->getGraph ()->isGloballyIndexed ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! NewFlag, std::invalid_argument, "Tpetra::CrsMatrix::"
        "transferAndFillComplete: The input argument 'destMat' is only allowed "
        "to be nonnull, if its graph is empty (neither locally nor globally "
        "indexed).");
      // FIXME (mfh 15 May 2014) At some point, we want to change
      // graphs and matrices so that their DistObject Map
      // (this->getMap()) may differ from their row Map.  This will
      // make redistribution for 2-D distributions more efficient.  I
      // hesitate to change this check, because I'm not sure how much
      // the code here depends on getMap() and getRowMap() being the
      // same.
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destMat->getRowMap ()->isSameAs (*MyRowMap), std::invalid_argument,
        "Tpetra::CrsMatrix::transferAndFillComplete: The (row) Map of the "
        "input argument 'destMat' is not the same as the (row) Map specified "
        "by the input argument 'rowTransfer'.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destMat->checkSizes (*this), std::invalid_argument,
        "Tpetra::CrsMatrix::transferAndFillComplete: You provided a nonnull "
        "destination matrix, but checkSizes() indicates that it is not a legal "
        "legal target for redistribution from the source matrix (*this).  This "
        "may mean that they do not have the same dimensions.");
    }

    // If forward mode (the default), then *this's (row) Map must be
    // the same as the source Map of the Transfer.  If reverse mode,
    // then *this's (row) Map must be the same as the target Map of
    // the Transfer.
    //
    // FIXME (mfh 15 May 2014) At some point, we want to change graphs
    // and matrices so that their DistObject Map (this->getMap()) may
    // differ from their row Map.  This will make redistribution for
    // 2-D distributions more efficient.  I hesitate to change this
    // check, because I'm not sure how much the code here depends on
    // getMap() and getRowMap() being the same.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (reverseMode || getRowMap ()->isSameAs (*rowTransfer.getSourceMap ())),
      std::invalid_argument, "Tpetra::CrsMatrix::transferAndFillComplete: "
      "rowTransfer->getSourceMap() must match this->getRowMap() in forward mode.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (! reverseMode || getRowMap ()->isSameAs (*rowTransfer.getTargetMap ())),
      std::invalid_argument, "Tpetra::CrsMatrix::transferAndFillComplete: "
      "rowTransfer->getTargetMap() must match this->getRowMap() in reverse mode.");

    // The basic algorithm here is:
    //
    // 1. Call the moral equivalent of "distor.do" to handle the import.
    // 2. Copy all the Imported and Copy/Permuted data into the raw
    //    CrsMatrix / CrsGraphData pointers, still using GIDs.
    // 3. Call an optimized version of MakeColMap that avoids the
    //    Directory lookups (since the importer knows who owns all the
    //    GIDs) AND reindexes to LIDs.
    // 4. Call expertStaticFillComplete()

    // Get information from the Importer
    const size_t NumSameIDs = rowTransfer.getNumSameIDs();
    ArrayView<const LO> ExportLIDs = reverseMode ?
      rowTransfer.getRemoteLIDs () : rowTransfer.getExportLIDs ();
    ArrayView<const LO> RemoteLIDs = reverseMode ?
      rowTransfer.getExportLIDs () : rowTransfer.getRemoteLIDs ();
    ArrayView<const LO> PermuteToLIDs = reverseMode ?
      rowTransfer.getPermuteFromLIDs () : rowTransfer.getPermuteToLIDs ();
    ArrayView<const LO> PermuteFromLIDs = reverseMode ?
      rowTransfer.getPermuteToLIDs () : rowTransfer.getPermuteFromLIDs ();
    Distributor& Distor = rowTransfer.getDistributor ();

    // Owning PIDs
    Teuchos::Array<int> SourcePids;
    Teuchos::Array<int> TargetPids;
    int MyPID = getComm ()->getRank ();

    // Temp variables for sub-communicators
    RCP<const map_type> ReducedRowMap, ReducedColMap,
      ReducedDomainMap, ReducedRangeMap;
    RCP<const Comm<int> > ReducedComm;

    // If the user gave us a null destMat, then construct the new
    // destination matrix.  We will replace its column Map later.
    if (destMat.is_null ()) {
      destMat = rcp (new this_type (MyRowMap, 0, StaticProfile, matrixparams));
    }

    /***************************************************/
    /***** 1) First communicator restriction phase ****/
    /***************************************************/
    if (restrictComm) {
      ReducedRowMap = MyRowMap->removeEmptyProcesses ();
      ReducedComm = ReducedRowMap.is_null () ?
        Teuchos::null :
        ReducedRowMap->getComm ();
      destMat->removeEmptyProcessesInPlace (ReducedRowMap);

      ReducedDomainMap = MyRowMap.getRawPtr () == MyDomainMap.getRawPtr () ?
        ReducedRowMap :
        MyDomainMap->replaceCommWithSubset (ReducedComm);
      ReducedRangeMap = MyRowMap.getRawPtr () == MyRangeMap.getRawPtr () ?
        ReducedRowMap :
        MyRangeMap->replaceCommWithSubset (ReducedComm);

      // Reset the "my" maps
      MyRowMap    = ReducedRowMap;
      MyDomainMap = ReducedDomainMap;
      MyRangeMap  = ReducedRangeMap;

      // Update my PID, if we've restricted the communicator
      if (! ReducedComm.is_null ()) {
        MyPID = ReducedComm->getRank ();
      }
      else {
        MyPID = -2; // For debugging
      }
    }
    else {
      ReducedComm = MyRowMap->getComm ();
    }

    /***************************************************/
    /***** 2) From Tpera::DistObject::doTransfer() ****/
    /***************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC ImportSetup"))));
#endif
    // Get the owning PIDs
    RCP<const import_type> MyImporter = getGraph ()->getImporter ();

    if (! restrictComm && ! MyImporter.is_null () &&
        BaseDomainMap->isSameAs (*getDomainMap ())) {
      // Same domain map as source matrix
      //
      // NOTE: This won't work for restrictComm (because the Import
      // doesn't know the restricted PIDs), though writing an
      // optimized version for that case would be easy (Import an
      // IntVector of the new PIDs).  Might want to add this later.
      Import_Util::getPids (*MyImporter, SourcePids, false);
    }
    else if (MyImporter.is_null () && BaseDomainMap->isSameAs (*getDomainMap ())) {
      // Matrix has no off-process entries
      SourcePids.resize (getColMap ()->getNodeNumElements ());
      SourcePids.assign (getColMap ()->getNodeNumElements (), MyPID);
    }
    else if (BaseDomainMap->isSameAs (*BaseRowMap) &&
             getDomainMap ()->isSameAs (*getRowMap ())) {
      // We can use the rowTransfer + SourceMatrix's Import to find out who owns what.
      IntVectorType TargetRow_pids (domainMap);
      IntVectorType SourceRow_pids (getRowMap ());
      IntVectorType SourceCol_pids (getColMap ());

      TargetRow_pids.putScalar (MyPID);
      if (! reverseMode && xferAsImport != NULL) {
        SourceRow_pids.doExport (TargetRow_pids, *xferAsImport, INSERT);
      }
      else if (reverseMode && xferAsExport != NULL) {
        SourceRow_pids.doExport (TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (! reverseMode && xferAsExport != NULL) {
        SourceRow_pids.doImport (TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (reverseMode && xferAsImport != NULL) {
        SourceRow_pids.doImport (TargetRow_pids, *xferAsImport, INSERT);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::CrsMatrix::"
          "transferAndFillComplete: Should never get here!  "
          "Please report this bug to a Tpetra developer.");
      }
      SourceCol_pids.doImport (SourceRow_pids, *MyImporter, INSERT);
      SourcePids.resize (getColMap ()->getNodeNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::CrsMatrix::"
        "transferAndFillComplete: This method only allows either domainMap == "
        "getDomainMap (), or (domainMap == rowTransfer.getTargetMap () and "
        "getDomainMap () == getRowMap ()).");
    }
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Pack-2"))));
#endif

    // Tpetra-specific stuff
    //
    // FIXME (mfh 15 May 2014) This should work fine if CrsMatrix
    // inherits from DistObject (in which case all arrays that get
    // resized here are Teuchos::Array), but it won't work if
    // CrsMatrix inherits from DistObjectKA (in which case all arrays
    // that get resized here are Kokkos::View).  In the latter case,
    // imports_ and numExportPacketsPerLID_ each have only a device
    // view, but numImportPacketsPerLID_ has a device view and a host
    // view (host_numImportPacketsPerLID_).
    //
    // Currently, CrsMatrix inherits from DistObject, not
    // DistObjectKA, so the code below should be fine for the Kokkos
    // refactor version of CrsMatrix.
    //
    // For this and for all other cases in this function that want to
    // resize the DistObject's communication arrays, it would make
    // sense to give DistObject (and DistObjectKA) methods for
    // resizing that don't expose the details of whether these are
    // Teuchos::Array or Kokkos::View.
    size_t constantNumPackets = destMat->constantNumberOfPackets ();
    if (constantNumPackets == 0) {
      destMat->numExportPacketsPerLID_old_.resize (ExportLIDs.size ());
      destMat->numImportPacketsPerLID_old_.resize (RemoteLIDs.size ());
    }
    else {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = RemoteLIDs.size() * constantNumPackets;
      if (static_cast<size_t> (destMat->imports_old_.size ()) != rbufLen) {
        destMat->imports_old_.resize (rbufLen);
      }
    }

    // Pack & Prepare w/ owning PIDs
    //
    // FIXME (mfh 15 May 2014) This should work fine if CrsMatrix
    // inherits from DistObject (in which case all arrays that get
    // passed in here are Teuchos::Array), but it won't work if
    // CrsMatrix inherits from DistObjectKA (in which case all arrays
    // that get passed in here are Kokkos::View).  In the latter case,
    // exports_ and numExportPacketsPerLID_ each have only a device
    // view.
    //
    // Currently, CrsMatrix inherits from DistObject, not
    // DistObjectKA, so the code below should be fine for the Kokkos
    // refactor version of CrsMatrix.
#ifdef HAVE_TPETRA_DEBUG
    {
      using Teuchos::outArg;
      using Teuchos::REDUCE_MAX;
      using Teuchos::reduceAll;
      using std::cerr;
      using std::endl;
      RCP<const Teuchos::Comm<int> > comm = this->getComm ();
      const int myRank = comm->getRank ();
      const int numProcs = comm->getSize ();

      std::ostringstream os;
      int lclErr = 0;
      try {
        Import_Util::packAndPrepareWithOwningPIDs (*this, ExportLIDs,
                                                   destMat->exports_old_,
                                                   destMat->numExportPacketsPerLID_old_ (),
                                                   constantNumPackets, Distor,
                                                   SourcePids);
      }
      catch (std::exception& e) {
        os << "Proc " << myRank << ": " << e.what ();
        lclErr = 1;
      }
      int gblErr = 0;
      if (! comm.is_null ()) {
        reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
      }
      if (gblErr != 0) {
        if (myRank == 0) {
          cerr << "packAndPrepareWithOwningPIDs threw an exception: " << endl;
        }
        std::ostringstream err;
        for (int r = 0; r < numProcs; ++r) {
          if (r == myRank && lclErr != 0) {
            cerr << os.str () << endl;
          }
          comm->barrier ();
          comm->barrier ();
          comm->barrier ();
        }

        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "packAndPrepareWithOwningPIDs threw an "
          "exception.");
      }
    }

#else
    Import_Util::packAndPrepareWithOwningPIDs (*this, ExportLIDs,
                                               destMat->exports_old_,
                                               destMat->numExportPacketsPerLID_old_ (),
                                               constantNumPackets, Distor,
                                               SourcePids);
#endif // HAVE_TPETRA_DEBUG

    // Do the exchange of remote data.
    //
    // FIXME (mfh 15 May 2014) This should work fine if CrsMatrix
    // inherits from DistObject (in which case all arrays that get
    // passed in here are Teuchos::Array), but it won't work if
    // CrsMatrix inherits from DistObjectKA (in which case all arrays
    // that get passed in here are Kokkos::View).
    //
    // In the latter case, imports_, exports_, and
    // numExportPacketsPerLID_ each have only a device view.
    // numImportPacketsPerLIDs_ is a device view, and also has a host
    // view (host_numImportPacketsPerLID_).
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Transfer"))));
#endif

    if (communication_needed) {
      if (reverseMode) {
        if (constantNumPackets == 0) { // variable number of packets per LID
          Distor.doReversePostsAndWaits (destMat->numExportPacketsPerLID_old_ ().getConst (), 1,
                                         destMat->numImportPacketsPerLID_old_ ());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < destMat->numImportPacketsPerLID_old_.size (); ++i) {
            totalImportPackets += destMat->numImportPacketsPerLID_old_[i];
          }
          destMat->imports_old_.resize (totalImportPackets);
          Distor.doReversePostsAndWaits (destMat->exports_old_ ().getConst (),
                                         destMat->numExportPacketsPerLID_old_ (),
                                         destMat->imports_old_ (),
                                         destMat->numImportPacketsPerLID_old_ ());
        }
        else { // constant number of packets per LID
          Distor.doReversePostsAndWaits (destMat->exports_old_ ().getConst (),
                                         constantNumPackets,
                                         destMat->imports_old_ ());
        }
      }
      else { // forward mode (the default)
        if (constantNumPackets == 0) { // variable number of packets per LID
          Distor.doPostsAndWaits (destMat->numExportPacketsPerLID_old_ ().getConst (), 1,
                                  destMat->numImportPacketsPerLID_old_ ());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < destMat->numImportPacketsPerLID_old_.size (); ++i) {
            totalImportPackets += destMat->numImportPacketsPerLID_old_[i];
          }
          destMat->imports_old_.resize (totalImportPackets);
          Distor.doPostsAndWaits (destMat->exports_old_ ().getConst (),
                                  destMat->numExportPacketsPerLID_old_ (),
                                  destMat->imports_old_ (),
                                  destMat->numImportPacketsPerLID_old_ ());
        }
        else { // constant number of packets per LID
          Distor.doPostsAndWaits (destMat->exports_old_ ().getConst (),
                                  constantNumPackets,
                                  destMat->imports_old_ ());
        }
      }
    }

    /*********************************************************************/
    /**** 3) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
    /*********************************************************************/

    // FIXME (mfh 15 May 2014) This should work fine if CrsMatrix
    // inherits from DistObject (in which case all arrays that get
    // passed in here are Teuchos::Array), but it won't work if
    // CrsMatrix inherits from DistObjectKA (in which case all arrays
    // that get passed in here are Kokkos::View).
    //
    // In the latter case, imports_ only has a device view.
    // numImportPacketsPerLIDs_ is a device view, and also has a host
    // view (host_numImportPacketsPerLID_).
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Unpack-1"))));
#endif
    size_t mynnz =
      Import_Util::unpackAndCombineWithOwningPIDsCount (*this, RemoteLIDs,
                                                        destMat->imports_old_ (),
                                                        destMat->numImportPacketsPerLID_old_ (),
                                                        constantNumPackets, Distor, INSERT,
                                                        NumSameIDs, PermuteToLIDs,
                                                        PermuteFromLIDs);
    size_t N = BaseRowMap->getNodeNumElements ();

    // Allocations
    ArrayRCP<size_t> CSR_rowptr(N+1);
    ArrayRCP<GO> CSR_colind_GID;
    ArrayRCP<LO> CSR_colind_LID;
    ArrayRCP<Scalar> CSR_vals;
    CSR_colind_GID.resize (mynnz);
    CSR_vals.resize (mynnz);

    // If LO and GO are the same, we can reuse memory when
    // converting the column indices from global to local indices.
    if (typeid (LO) == typeid (GO)) {
      CSR_colind_LID = Teuchos::arcp_reinterpret_cast<LO> (CSR_colind_GID);
    }
    else {
      CSR_colind_LID.resize (mynnz);
    }

    // FIXME (mfh 15 May 2014) This should work fine if CrsMatrix
    // inherits from DistObject (in which case all arrays that get
    // passed in here are Teuchos::Array), but it won't work if
    // CrsMatrix inherits from DistObjectKA (in which case all arrays
    // that get passed in here are Kokkos::View).
    //
    // In the latter case, imports_ only has a device view.
    // numImportPacketsPerLIDs_ is a device view, and also has a host
    // view (host_numImportPacketsPerLID_).
    //
    // FIXME (mfh 15 May 2014) Why can't we abstract this out as an
    // unpackAndCombine method on a "CrsArrays" object?  This passing
    // in a huge list of arrays is icky.  Can't we have a bit of an
    // abstraction?  Implementing a concrete DistObject subclass only
    // takes five methods.
    Import_Util::unpackAndCombineIntoCrsArrays (*this, RemoteLIDs, destMat->imports_old_ (),
                                                destMat->numImportPacketsPerLID_old_ (),
                                                constantNumPackets, Distor, INSERT, NumSameIDs,
                                                PermuteToLIDs, PermuteFromLIDs, N, mynnz, MyPID,
                                                CSR_rowptr (), CSR_colind_GID (), CSR_vals (),
                                                SourcePids (), TargetPids);

    /**************************************************************/
    /**** 4) Call Optimized MakeColMap w/ no Directory Lookups ****/
    /**************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Unpack-2"))));
#endif
    // Call an optimized version of makeColMap that avoids the
    // Directory lookups (since the Import object knows who owns all
    // the GIDs).
    Teuchos::Array<int> RemotePids;
    Import_Util::lowCommunicationMakeColMapAndReindex (CSR_rowptr (),
                                                       CSR_colind_LID (),
                                                       CSR_colind_GID (),
                                                       BaseDomainMap,
                                                       TargetPids, RemotePids,
                                                       MyColMap);

    /*******************************************************/
    /**** 4) Second communicator restriction phase      ****/
    /*******************************************************/
    if (restrictComm) {
      ReducedColMap = (MyRowMap.getRawPtr () == MyColMap.getRawPtr ()) ?
        ReducedRowMap :
        MyColMap->replaceCommWithSubset (ReducedComm);
      MyColMap = ReducedColMap; // Reset the "my" maps
    }

    // Replace the col map
    destMat->replaceColMap (MyColMap);

    // Short circuit if the processor is no longer in the communicator
    //
    // NOTE: Epetra replaces modifies all "removed" processes so they
    // have a dummy (serial) Map that doesn't touch the original
    // communicator.  Duplicating that here might be a good idea.
    if (ReducedComm.is_null ()) {
      return;
    }

    /***************************************************/
    /**** 5) Sort                                   ****/
    /***************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC Unpack-3"))));
#endif
    Import_Util::sortCrsEntries (CSR_rowptr (),
                                 CSR_colind_LID (),
                                 CSR_vals ());
    if ((! reverseMode && xferAsImport != NULL) ||
        (reverseMode && xferAsExport != NULL)) {
      Import_Util::sortCrsEntries (CSR_rowptr (),
                                   CSR_colind_LID (),
                                   CSR_vals ());
    }
    else if ((! reverseMode && xferAsExport != NULL) ||
             (reverseMode && xferAsImport != NULL)) {
      Import_Util::sortAndMergeCrsEntries (CSR_rowptr (),
                                           CSR_colind_LID (),
                                           CSR_vals ());
      if (CSR_rowptr[N] != mynnz) {
        CSR_colind_LID.resize (CSR_rowptr[N]);
        CSR_vals.resize (CSR_rowptr[N]);
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Tpetra::CrsMatrix::"
        "transferAndFillComplete: Should never get here!  "
        "Please report this bug to a Tpetra developer.");
    }
    /***************************************************/
    /**** 6) Reset the colmap and the arrays        ****/
    /***************************************************/

    // Call constructor for the new matrix (restricted as needed)
    //
    // NOTE (mfh 15 May 2014) This should work fine for the Kokkos
    // refactor version of CrsMatrix, though it reserves the right to
    // make a deep copy of the arrays.
    destMat->setAllValues (CSR_rowptr, CSR_colind_LID, CSR_vals);

    /***************************************************/
    /**** 7) Build Importer & Call ESFC             ****/
    /***************************************************/
    // Pre-build the importer using the existing PIDs
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC ESFC"))));
#endif
    RCP<import_type> MyImport = rcp (new import_type (MyDomainMap, MyColMap, RemotePids));
    destMat->expertStaticFillComplete (MyDomainMap, MyRangeMap, MyImport);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >& destMatrix,
                         const import_type& importer,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, importer, domainMap, rangeMap, params);
  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  CrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >& destMatrix,
                         const export_type& exporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, exporter, domainMap, rangeMap, params);
  }

} // namespace Tpetra

#endif
