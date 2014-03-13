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

/// \file Tpetra_Map_def.hpp
///
/// Implementation of the methods of Tpetra::Map, and of related
/// nonmember constructors for Tpetra::Map.

#ifndef TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP
#define TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP

#include <Tpetra_Directory.hpp> // must include for implicit instantiation to work
#include <Tpetra_Util.hpp>
#include <Teuchos_as.hpp>
#include <stdexcept>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Map_decl.hpp"
#endif

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const LocalGlobal lOrG,
       const Teuchos::RCP<node_type>& node) :
    comm_ (comm),
    node_ (node)
  {
    // Start with a host Map.  We could create the device Map on
    // demand, but it's easier to create it here.
    //
    // FIXME (mfh 06 Feb 2014) If we're using UVM, we don't really
    // need the host and device Maps to be separate.
    mapHost_ = host_impl_type (Teuchos::as<GlobalOrdinal> (globalNumIndices),
                               indexBase, *comm, lOrG);
    mapDevice_.create_copy_view (mapHost_);

    // Create the Directory on demand in getRemoteIndexList().
    //setupDirectory ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const size_t myNumIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const Teuchos::RCP<node_type>& node) :
    comm_ (comm),
    node_ (node)
  {
    typedef GlobalOrdinal GO;
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();

    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);
    const GO myNumInds = (myNumIndices == GSTI) ?
      getInvalidLocalIndex () : Teuchos::as<GO> (myNumIndices);

    // Start with a host Map.  We could create the device Map on
    // demand, but it's easier to create it here.
    //
    // FIXME (mfh 06 Feb 2014) If we're using UVM, we don't really
    // need the host and device Maps to be separate.
    mapHost_ = host_impl_type (globalNumInds, myNumInds, indexBase, *comm);
    mapDevice_.create_copy_view (mapHost_);

    // Create the Directory on demand in getRemoteIndexList().
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const Kokkos::View<const GlobalOrdinal*, device_type>& myGlobalIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const Teuchos::RCP<node_type>& node) :
    comm_ (comm),
    node_ (node)
  {
    typedef GlobalOrdinal GO;
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);

    // Since we already have device data, start here with a device
    // Map, and then create the host Map.
    //
    // FIXME (mfh 06 Feb 2014) If we're using UVM, we don't really
    // need the host and device Maps to be separate.
    mapDevice_ = device_impl_type (globalNumInds, myGlobalIndices, indexBase, *comm);
    mapHost_.create_copy_view (mapDevice_);

    // Create the Directory on demand in getRemoteIndexList().
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const Teuchos::RCP<node_type> &node) :
    comm_ (comm),
    node_ (node)
  {
    typedef GlobalOrdinal GO;
    typedef Kokkos::View<const GlobalOrdinal*, device_type,
                         Kokkos::MemoryUnmanaged> host_view_type;
    typedef Kokkos::View<GlobalOrdinal*, device_type> device_view_type;

    // Copy the input GID list from host (we assume that
    // Teuchos::ArrayView should only view host memory) to device.
    //
    // FIXME (mfh 06 Feb 2014) We could use the CUDA API function here
    // that can look at a pointer and tell whether it lives on host or
    // device, to tell whether the Teuchos::ArrayView is viewing host
    // or device memory.
    host_view_type gidsHost (myGlobalIndices.getRawPtr (), myGlobalIndices.size ());
    device_view_type gidsDevice ("GIDs", myGlobalIndices.size ());
    Kokkos::deep_copy (gidsDevice, gidsHost);

    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);

    mapHost_ = host_impl_type (globalNumInds, gidsDevice, indexBase, *comm);
    // FIXME (mfh 06 Feb 2014) If we're using UVM, we don't really
    // need the host and device Maps to be separate.
    mapDevice_.template create_copy_view<host_mirror_device_type> (mapHost_);

    // Create the Directory on demand in getRemoteIndexList().
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalNumElements () const {
    return mapHost_.getGlobalNumIndices ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeNumElements () const {
    return mapHost_.getMyNumIndices ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getIndexBase () const {
    return mapHost_.getIndexBase ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getInvalidGlobalIndex () const {
    return mapHost_.getInvalidGlobalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LocalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getInvalidLocalIndex () const {
    return mapHost_.getInvalidLocalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LocalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMinLocalIndex () const {
    return static_cast<LocalOrdinal> (0);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LocalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMaxLocalIndex () const {
    return mapHost_.getMaxLocalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMinGlobalIndex () const {
    return mapHost_.getMinGlobalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMaxGlobalIndex () const {
    return mapHost_.getMaxGlobalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMinAllGlobalIndex () const {
    return mapHost_.getMinAllGlobalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getMaxAllGlobalIndex () const {
    return mapHost_.getMaxAllGlobalIndex ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LocalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalElement (const GlobalOrdinal globalIndex) const
  {
    return mapHost_.getLocalIndex (globalIndex);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalElement (const LocalOrdinal localIndex) const
  {
    return mapHost_.getGlobalIndex (localIndex);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isNodeLocalElement (const LocalOrdinal localIndex) const {
    return mapHost_.isOwnedLocalIndex (localIndex);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isNodeGlobalElement (const GlobalOrdinal globalIndex) const {
    return mapHost_.isOwnedGlobalIndex (globalIndex);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isUniform () const {
    return mapHost_.isUniform ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isContiguous () const {
    return mapHost_.isContiguous ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, bail out with an exception if the two
    // communicators don't have the same numbers of processes.
    // This is explicitly forbidden by the public documentation.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getComm ()->getSize () != map.getComm ()->getSize (),
      std::invalid_argument, "Tpetra::Map::isCompatibile: The two Maps' "
      "communicators must have the same numbers of processes in order to call "
      "this method.");
#endif // HAVE_TPETRA_DEBUG

    // Pointer equality on one process always implies pointer equality
    // on all processes, since Map is immutable.
    if (this == &map) {
      return true;
    }

    // Do both Maps have the same number of elements, both globally
    // and on the calling process?
    int locallyCompat = 0;
    if (getGlobalNumElements() != map.getGlobalNumElements() ||
        getNodeNumElements() != map.getNodeNumElements()) {
      locallyCompat = 0; // NOT compatible on this process
    }
    else {
      locallyCompat = 1; // compatible on this process
    }

    int globallyCompat = 0;
    reduceAll<int, int> (*comm_, REDUCE_MIN, locallyCompat, outArg (globallyCompat));
    return (globallyCompat == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, bail out with an exception if the two
    // communicators don't have the same numbers of processes.
    // This is explicitly forbidden by the public documentation.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getComm ()->getSize () != map.getComm ()->getSize (),
      std::invalid_argument, "Tpetra::Map::isSameAs: The two Maps' communicators"
      "must have the same numbers of processes in order to call this method.");
#endif // HAVE_TPETRA_DEBUG

    if (this == &map) {
      // If the input Map is the same object (has the same address) as
      // *this, then the Maps are the same.  We assume that if this
      // holds on one process, then it holds on all processes.
      return true;
    }

    // Check all other known properties that are the same on all
    // processes.  If the two Maps do not share any of these
    // properties, then they cannot be the same.
    if (getMinAllGlobalIndex () != map.getMinAllGlobalIndex () ||
        getMaxAllGlobalIndex () != map.getMaxAllGlobalIndex () ||
        getGlobalNumElements () != map.getGlobalNumElements () ||
        isDistributed () != map.isDistributed () ||
        getIndexBase () != map.getIndexBase ()) {
      return false;
    }

    // If we get this far, we need to check local properties and the
    // communicate same-ness across all processes
    // we prefer local work over communication, ergo, we will perform all
    // comparisons and conduct a single communication
    int isSame_lcl = 1;

    // The two communicators must have the same number of processes,
    // with process ranks occurring in the same order.
    if (! Details::congruent (*comm_, * (map.getComm ()))) {
      isSame_lcl = 0;
    }

    // Check the number of entries owned by this process.
    if (getNodeNumElements () != map.getNodeNumElements ()) {
      isSame_lcl = 0;
    }

    if (getNodeNumElements () == 0 && map.getNodeNumElements () == 0) {
      isSame_lcl = 1; // both Maps have no GIDs on the calling process
    }
    else if (isSame_lcl != 0) {
      // Check that the entries owned by this process in both Maps are
      // the same.  Only do this if we haven't already determined that
      // the Maps aren't the same.  If the Maps are contiguous, we can
      // check the ranges easily by looking at the min and max GID on
      // this process.  Otherwise, we'll compare their GID lists.
      if (isContiguous () && map.isContiguous ()) {
        if (getMinGlobalIndex () != map.getMinGlobalIndex () ||
            getMaxGlobalIndex() != map.getMaxGlobalIndex()) {
          isSame_lcl = 0;
        }
      }
      else {
        // We could be more clever here to avoid calling
        // getNodeElementList() on either of the two Maps has
        // contiguous GIDs.  For now, we call it regardless of whether
        // the Map is contiguous, as long as one of the Maps is not
        // contiguous.
        //
        // std::equal requires that the latter range is as large as
        // the former.  We know the ranges have equal length, because
        // they have the same number of local entries.
        Teuchos::ArrayView<const GlobalOrdinal> ge1 =     getNodeElementList ();
        Teuchos::ArrayView<const GlobalOrdinal> ge2 = map.getNodeElementList ();

        TEUCHOS_TEST_FOR_EXCEPTION(
          ge1.getRawPtr () == NULL && this->getNodeNumElements () != 0,
          std::logic_error,
          "Tpetra::Map::isSameAs: this->getNodeElementList() returns null, "
          "but *this Map on the calling process has "
          << this->getNodeNumElements () << " > 0 indices.");

        TEUCHOS_TEST_FOR_EXCEPTION(
          ge1.getRawPtr () == NULL && this->getNodeNumElements () != 0,
          std::logic_error,
          "Tpetra::Map::isSameAs: map.getNodeElementList() returns null, "
          "but the input Map on the calling process has "
          << this->getNodeNumElements () << " > 0 indices.");

        if (! std::equal (ge1.begin (), ge1.end (), ge2.begin ())) {
          isSame_lcl = 0;
        }
      }
    }

    // Return true if and only if all processes report "same-ness."
    int isSame_gbl = 0;
    reduceAll<int, int> (*comm_, REDUCE_MIN, isSame_lcl, outArg (isSame_gbl));
    return isSame_gbl == 1;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayView<const GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNodeElementList () const
  {
    typedef GlobalOrdinal GO;
    Kokkos::View<const GO*, host_mirror_device_type> myGlobalInds =
      mapHost_.getMyGlobalIndices (); // creates it if it doesn't exist

    // FIXME (mfh 06 Feb 2014) It's called something else other than getRawPtr.
    return Teuchos::ArrayView<const GO> (myGlobalInds.ptr_on_device (),
                                         myGlobalInds.dimension_0 ());
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isDistributed () const {
    return mapHost_.isDistributed ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  description () const {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;

    os << "\"Tpetra::Map\": {"
       << "LocalOrdinalType: " << TypeNameTraits<local_ordinal_type>::name ()
       << ", GlobalOrdinalType: " << TypeNameTraits<global_ordinal_type>::name ()
       << ", DeviceType: " << TypeNameTraits<device_type>::name ();
    if (this->getObjectLabel () != "") {
      os << ", Label: \"" << this->getObjectLabel () << "\"";
    }
    os << ", Global number of entries: " << getGlobalNumElements ()
       << ", Number of processes: " << getComm ()->getSize ()
       << ", Uniform: " << (isUniform () ? "true" : "false")
       << ", Contiguous: " << (isContiguous () ? "true" : "false")
       << ", Distributed: " << (isDistributed () ? "true" : "false")
       << "}";
    return os.str ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::OSTab;
    using Teuchos::toString;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    typedef typename ArrayView<const GlobalOrdinal>::size_type size_type;

    const size_t nME = getNodeNumElements ();
    ArrayView<const GlobalOrdinal> myEntries = getNodeElementList ();
    const int myRank = comm_->getRank ();
    const int numProcs = comm_->getSize ();

    const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    // By convention, describe() always begins with a tab before printing.
    OSTab tab0 (out);

    if (vl == VERB_NONE) {
      // do nothing
    }
    else if (vl == VERB_LOW) {
      if (myRank == 0) {
        out << "\"Tpetra::Map\":" << endl;
        OSTab tab1 (out);
        out << "LocalOrdinalType: " << TypeNameTraits<LocalOrdinal>::name () << endl
            << "GlobalOrdinalType: " << TypeNameTraits<GlobalOrdinal>::name () << endl
            << "DeviceType: " << TypeNameTraits<device_type>::name () << endl;
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
        out << "Global number of entries: " << getGlobalNumElements () << endl
            << "Minimum global index: " << getMinAllGlobalIndex () << endl
            << "Maximum global index: " << getMaxAllGlobalIndex () << endl
            << "Index base: " << getIndexBase () << endl
            << "Number of processes: " << getComm ()->getSize () << endl
            << "Uniform: " << (isUniform () ? "true" : "false") << endl
            << "Contiguous: " << (isContiguous () ? "true" : "false") << endl
            << "Distributed: " << (isDistributed () ? "true" : "false") << endl;
      }
    }

    if (vl >= VERB_HIGH) { // HIGH or EXTREME
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          out << "Process " << myRank << ":" << endl;
          OSTab tab1 (out);
          out << "My number of entries: " << nME << endl
              << "My minimum global index: " << getMinGlobalIndex () << endl
              << "My maximum global index: " << getMaxGlobalIndex () << endl;
          if (vl == VERB_EXTREME) {
            out << "My global indices: [";
            for (size_type k = 0; k < myEntries.size (); ++k) {
              out << myEntries[k];
              if (k + 1 < myEntries.size ()) {
                out << ", ";
              }
            }
            out << "]" << endl;
          }
          std::flush (out);
        }
        // Do a few global ops to give I/O a chance to complete
        comm_->barrier ();
        comm_->barrier ();
        comm_->barrier ();
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  replaceCommWithSubset (const Teuchos::RCP<const Teuchos::Comm<int> >& newComm) const
  {
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef global_size_t GST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Map<LO, GO, node_type> map_type;

    // mfh 26 Mar 2013: The lazy way to do this is simply to recreate
    // the Map by calling its ordinary public constructor, using the
    // original Map's data.  This only involves O(1) all-reduces over
    // the new communicator, which in the common case only includes a
    // small number of processes.

    // Create the Map to return.
    if (newComm.is_null ()) {
      return Teuchos::null; // my process does not participate in the new Map
    } else {
      // Map requires that the index base equal the global min GID.
      // Figuring out the global min GID requires a reduction over all
      // processes in the new communicator.  It could be that some (or
      // even all) of these processes contain zero entries.  (Recall
      // that this method, unlike removeEmptyProcesses(), may remove
      // an arbitrary subset of processes.)  We deal with this by
      // doing a min over the min GID on each process if the process
      // has more than zero entries, or the global max GID, if that
      // process has zero entries.  If no processes have any entries,
      // then the index base doesn't matter anyway.
      const GO myMinGid = (this->getNodeNumElements () == 0) ?
        this->getMaxAllGlobalIndex () : this->getMinGlobalIndex ();
      GO newIndexBase = this->getInvalidGlobalIndex ();
      reduceAll<int, GO> (*newComm, REDUCE_MIN, myMinGid, outArg (newIndexBase));

      // Make Map's constructor compute the global number of indices.
      const GST globalNumInds = Teuchos::OrdinalTraits<GST>::invalid ();
      Kokkos::View<const GlobalOrdinal*, DeviceType> myGIDs =
        mapDevice_.getMyGlobalIndices ();
      return rcp (new map_type (globalNumInds, myGIDs, newIndexBase,
                                newComm, this->getNode ()));
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  removeEmptyProcesses () const
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    // Create the new communicator.  split() returns a valid
    // communicator on all processes.  On processes where color == 0,
    // ignore the result.  Passing key == 0 tells MPI to order the
    // processes in the new communicator by their rank in the old
    // communicator.
    const int color = (getNodeNumElements () == 0) ? 0 : 1;
    // MPI_Comm_split must be called collectively over the original
    // communicator.  We can't just call it on processes with color
    // one, even though we will ignore its result on processes with
    // color zero.
    RCP<const Comm<int> > newComm = comm_->split (color, 0);
    if (color == 0) {
      newComm = null;
    }

    // Create the Map to return.
    if (newComm.is_null ()) {
      return null; // my process does not participate in the new Map
    } else {
      // The default constructor that's useful for clone() above is
      // also useful here.
      RCP<Map> map    = rcp (new Map ());
      map->comm_      = newComm;
      map->node_      = node_;
      map->mapHost_   = mapHost_;
      map->mapDevice_ = mapDevice_;

      // Uniformity and contiguity have not changed.  The directory
      // has changed, but we've taken care of that above.  However,
      // distributed-ness may have changed, since the communicator has
      // changed.
      //
      // If the original Map was NOT distributed, then the new Map
      // cannot be distributed.  If the number of processes in the new
      // communicator is 1, then the new Map is not distributed.
      // Otherwise, we have to check the new Map using an all-reduce
      // (over the new communicator).  For example, the original Map
      // may have had some processes with zero elements, and all other
      // processes with the same number of elements as in the whole
      // Map.  That Map is technically distributed, because of the
      // processes with zero elements.  Removing those processes would
      // make the new Map locally replicated.
      if (! isDistributed () || newComm->getSize () == 1) {
        map->mapHost_.setDistributed (false);
        map->mapDevice_.setDistributed (false);
      } else {
        const int iOwnAllGids =
          (getNodeNumElements () == getGlobalNumElements ()) ? 1 : 0;
        int allProcsOwnAllGids = 0;
        reduceAll<int, int> (*newComm, REDUCE_MIN, iOwnAllGids,
                             outArg (allProcsOwnAllGids));
        map->mapHost_.setDistributed (allProcsOwnAllGids != 1);
        map->mapDevice_.setDistributed (allProcsOwnAllGids != 1);
      }

      // The Directory will be created on demand in getRemoteIndexList().
      //
      // FIXME (mfh 26 Mar 2013) It should be possible to "filter" the
      // directory more efficiently than just recreating it.  If
      // directory recreation proves a bottleneck, we can always
      // revisit this.  On the other hand, Directory creation is only
      // collective over the new, presumably much smaller
      // communicator, so it may not be worth the effort to optimize.
      map->directory_ = null;
      return map;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  setupDirectory () const
  {
    using Teuchos::rcp;
    typedef Directory<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > directory_type;
    // Only create the Directory if it hasn't been created yet.
    // This is a collective operation.
    if (directory_.is_null ()) {
      directory_ = rcp (new directory_type (*this));
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LookupStatus
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal> & GIDs,
                      const Teuchos::ArrayView<int> & PIDs,
                      const Teuchos::ArrayView<LocalOrdinal> & LIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getComm ().is_null (), std::logic_error, "Tpetra::Map (Kokkos "
      "refactor)::getRemoteIndexList (3 args): getComm() returns null.  "
      "Please report this bug to the Tpetra developers.");

    // mfh 03 Mar 2014: The exception test below should be commented
    // out; it's valid to give getRemoteIndexList GIDs that the Map
    // doesn't own, and a Map with zero GIDs doesn't own any GIDs.
    //
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   GIDs.size () > 0 && getGlobalNumElements () == 0, std::runtime_error,
    //   "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): The Map has "
    //   "zero entries (globally), so you may not call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != PIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): GIDs.size ()"
      " = " << GIDs.size () << " != PIDs.size () = " << PIDs.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != LIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): GIDs.size ()"
      " = " << GIDs.size () << " != LIDs.size () = " << LIDs.size () << ".");

    // getRemoteIndexList must be called collectively, and Directory
    // creation is collective too, so it's OK to create the Directory
    // on demand.
    setupDirectory ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      directory_.is_null (), std::logic_error,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): "
      "setupDirectory() failed to construct the directory.  "
      "Please report this bug to the Tpetra developers.");
    return directory_->getDirectoryEntries (*this, GIDs, PIDs, LIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LookupStatus
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDs,
                      const Teuchos::ArrayView<int> & PIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getComm ().is_null (), std::logic_error, "Tpetra::Map (Kokkos "
      "refactor)::getRemoteIndexList (2 args): getComm() returns null.  "
      "Please report this bug to the Tpetra developers.");

    // mfh 03 Mar 2014: SubmapImport test 2 actually triggers the
    // commented-out exception test below.  If I leave it commented
    // out, the test passes.  (It should be commented out; it's valid
    // to give getRemoteIndexList GIDs that the Map doesn't own.)
    //
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   GIDs.size () > 0 && getGlobalNumElements () == 0, std::runtime_error,
    //   "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (2 args): The Map has "
    //   "zero entries (globally), so you may not call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != PIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (2 args): GIDs.size ()"
      " = " << GIDs.size () << " != PIDs.size () = " << PIDs.size () << ".");

    // getRemoteIndexList must be called collectively, and Directory
    // creation is collective too, so it's OK to create the Directory
    // on demand.
    setupDirectory ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      directory_.is_null (), std::logic_error,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (2 args): "
      "setupDirectory() failed to construct the directory.  "
      "Please report this bug to the Tpetra developers.");
    return directory_->getDirectoryEntries (*this, GIDs, PIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Teuchos::Comm<int> >
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getComm () const {
    return comm_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNode () const {
    return node_;
  }

} // namespace Tpetra

#endif // TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP
