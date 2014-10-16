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

#ifndef TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP
#define TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP

#include <Tpetra_Directory.hpp> // must include for implicit instantiation to work
#include <Tpetra_Util.hpp>
#include <Teuchos_as.hpp>
#include <stdexcept>

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_Map_decl.hpp"
#endif

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map () :
    comm_ (new Teuchos::SerialComm<int> ()),
    node_ (KokkosClassic::Details::getNode<Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > ()),
    directory_ (new Directory<LocalOrdinal, GlobalOrdinal, node_type> ())
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const LocalGlobal lOrG,
       const Teuchos::RCP<node_type>& node) :
    comm_ (comm),
    node_ (node),
    directory_ (new Directory<LocalOrdinal, GlobalOrdinal, node_type> ())
  {
    // Start with a host Map implementation, since this will make this
    // class' public (host) methods work.  If users want device
    // methods, they will call getDeviceView(), which will initialize
    // the device Map implementation.
    //
    // NOTE (mfh 06 Feb 2014) If we're using UVM, we don't really need
    // the host and device Maps to be separate.
    mapHost_ = host_impl_type (Teuchos::as<GlobalOrdinal> (globalNumIndices),
                               indexBase, *comm, lOrG);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Map (const global_size_t globalNumIndices,
       const size_t myNumIndices,
       const GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
       const Teuchos::RCP<node_type>& node) :
    comm_ (comm),
    node_ (node),
    directory_ (new Directory<LocalOrdinal, GlobalOrdinal, node_type> ())
  {
    typedef GlobalOrdinal GO;
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();

    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);
    const GO myNumInds = (myNumIndices == GSTI) ?
      getInvalidLocalIndex () : Teuchos::as<GO> (myNumIndices);

    // Start with a host Map implementation, since this will make this
    // class' public (host) methods work.  If users want device
    // methods, they will call getDeviceView(), which will initialize
    // the device Map implementation.
    //
    // NOTE (mfh 06 Feb 2014) If we're using UVM, we don't really need
    // the host and device Maps to be separate.
    mapHost_ = host_impl_type (globalNumInds, myNumInds, indexBase, *comm);

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
    node_ (node),
    directory_ (new Directory<LocalOrdinal, GlobalOrdinal, node_type> ())
  {
    typedef GlobalOrdinal GO;
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);

    // Since we already have device data, start here with a device Map
    // implementation.  In order to make this class' host methods work
    // and still be const (that is, legal to call in a host parallel
    // operation), we can't create the host mirror on demand; we must
    // create it here, in advance.
    //
    // NOTE (mfh 06 Feb 2014) If we're using UVM, we don't really need
    // the host and device Maps to be separate.
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
    node_ (node),
    directory_ (new Directory<LocalOrdinal, GlobalOrdinal, node_type> ())
  {
    typedef GlobalOrdinal GO;
    typedef Kokkos::View<GlobalOrdinal*, device_type> device_view_type;
    typedef Kokkos::View<GlobalOrdinal*, typename device_view_type::host_mirror_space> host_view_type;
    typedef Kokkos::View<const GlobalOrdinal*, typename device_view_type::host_mirror_space> const_host_view_type;

    // Copy the input GID list from host (we assume that
    // Teuchos::ArrayView should only view host memory) to device.
    //
    // FIXME (mfh 06 Feb, 24 Mar 2014) We could use the CUDA API
    // function here that can look at a pointer and tell whether it
    // lives on host or device, to tell whether the Teuchos::ArrayView
    // is viewing host or device memory.  Regardless, we don't own the
    // data and we will need a deep copy anyway, so we might as well
    // copy it.
    const_host_view_type gidsHost_view (Kokkos::ViewWithoutManaging(),myGlobalIndices.getRawPtr (), myGlobalIndices.size ());
    host_view_type gidsHost ("GIDS_Host", myGlobalIndices.size ());
    Kokkos::deep_copy (gidsHost, gidsHost_view);
    device_view_type gidsDevice ("GIDs", myGlobalIndices.size ());
    Kokkos::deep_copy (gidsDevice, gidsHost);

    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const GO globalNumInds = (globalNumIndices == GSTI) ?
      getInvalidGlobalIndex () : Teuchos::as<GO> (globalNumIndices);

    // Start with a host Map implementation, since this will make this
    // class' public (host) methods work.  If users want device
    // methods, they will call getDeviceView(), which will initialize
    // the device Map implementation.
    //
    // NOTE (mfh 06 Feb 2014) If we're using UVM, we don't really need
    // the host and device Maps to be separate, but we need to create a compatible
    // host view of the GIDs using the device ptr.
    mapDevice_ = device_impl_type (globalNumInds, gidsDevice , indexBase, *comm);
    #ifdef KOKKOS_USE_CUDA_UVM
      const GlobalOrdinal* ptr = gidsDevice.ptr_on_device();
      mapHost_ = host_impl_type (globalNumInds, const_host_view_type(Kokkos::ViewWithoutManaging(),ptr,myGlobalIndices.size ()), indexBase, *comm);
    #else
      mapHost_ = host_impl_type (globalNumInds, gidsHost, indexBase, *comm);
    #endif
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
  isCompatible (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    //
    // Tests that avoid the Boolean all-reduce below by using
    // globally consistent quantities.
    //
    if (this == &map) {
      // Pointer equality on one process always implies pointer
      // equality on all processes, since Map is immutable.
      return true;
    }
    else if (getComm ()->getSize () != map.getComm ()->getSize ()) {
      // The two communicators have different numbers of processes.
      // It's not correct to call isCompatible() in that case.  This
      // may result in the all-reduce hanging below.
      return false;
    }
    else if (getGlobalNumElements () != map.getGlobalNumElements ()) {
      // Two Maps are definitely NOT compatible if they have different
      // global numbers of indices.
      return false;
    }
    else if (isContiguous () && isUniform () &&
             map.isContiguous () && map.isUniform ()) {
      // Contiguous uniform Maps with the same number of processes in
      // their communicators, and with the same global numbers of
      // indices, are always compatible.
      return true;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      getGlobalNumElements () != map.getGlobalNumElements (), std::logic_error,
      "Tpetra::Map::isCompatible: There's a bug in this method.  We've already "
      "checked that this condition is true above, but it's false here.  "
      "Please report this bug to the Tpetra developers.");

    // Do both Maps have the same number of indices on each process?
    const int locallyCompat =
      (getNodeNumElements () == map.getNodeNumElements ()) ? 1 : 0;

    int globallyCompat = 0;
    reduceAll<int, int> (*comm_, REDUCE_MIN, locallyCompat, outArg (globallyCompat));
    return (globallyCompat == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  locallySameAs (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const
  {
    using Teuchos::ArrayView;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    // If both Maps are contiguous, we can compare their GID ranges
    // easily by looking at the min and max GID on this process.
    // Otherwise, we'll compare their GID lists.  If only one Map is
    // contiguous, then we only have to call getNodeElementList() on
    // the noncontiguous Map.  (It's best to avoid calling it on a
    // contiguous Map, since it results in unnecessary storage that
    // persists for the lifetime of the Map.)

    if (getNodeNumElements () != map.getNodeNumElements ()) {
      return false;
    }
    else if (getMinGlobalIndex () != map.getMinGlobalIndex () ||
             getMaxGlobalIndex () != map.getMaxGlobalIndex ()) {
      return false;
    }
    else {
      if (isContiguous ()) {
        if (map.isContiguous ()) {
          return true; // min and max match, so the ranges match.
        }
        else { // *this is contiguous, but map is not contiguous
          TEUCHOS_TEST_FOR_EXCEPTION(
            ! this->isContiguous () || map.isContiguous (), std::logic_error,
            "Tpetra::Map::locallySameAs: BUG");
          ArrayView<const GO> rhsElts = map.getNodeElementList ();
          const GO minLhsGid = this->getMinGlobalIndex ();
          const size_type numRhsElts = rhsElts.size ();
          for (size_type k = 0; k < numRhsElts; ++k) {
            const GO curLhsGid = minLhsGid + static_cast<GO> (k);
            if (curLhsGid != rhsElts[k]) {
              return false; // stop on first mismatch
            }
          }
          return true;
        }
      }
      else if (map.isContiguous ()) { // *this is not contiguous, but map is
        TEUCHOS_TEST_FOR_EXCEPTION(
          this->isContiguous () || ! map.isContiguous (), std::logic_error,
          "Tpetra::Map::locallySameAs: BUG");
        ArrayView<const GO> lhsElts = this->getNodeElementList ();
        const GO minRhsGid = map.getMinGlobalIndex ();
        const size_type numLhsElts = lhsElts.size ();
        for (size_type k = 0; k < numLhsElts; ++k) {
          const GO curRhsGid = minRhsGid + static_cast<GO> (k);
          if (curRhsGid != lhsElts[k]) {
            return false; // stop on first mismatch
          }
        }
        return true;
      }
      else { // neither *this nor map are contiguous
        // std::equal requires that the latter range is as large as
        // the former.  We know the ranges have equal length, because
        // they have the same number of local entries.
        ArrayView<const GO> lhsElts =     getNodeElementList ();
        ArrayView<const GO> rhsElts = map.getNodeElementList ();
        return std::equal (lhsElts.begin (), lhsElts.end (), rhsElts.begin ());
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    //
    // Tests that avoid the Boolean all-reduce below by using
    // globally consistent quantities.
    //
    if (this == &map) {
      // Pointer equality on one process always implies pointer
      // equality on all processes, since Map is immutable.
      return true;
    }
    else if (getComm ()->getSize () != map.getComm ()->getSize ()) {
      // The two communicators have different numbers of processes.
      // It's not correct to call isSameAs() in that case.  This
      // may result in the all-reduce hanging below.
      return false;
    }
    else if (getGlobalNumElements () != map.getGlobalNumElements ()) {
      // Two Maps are definitely NOT the same if they have different
      // global numbers of indices.
      return false;
    }
    else if (getMinAllGlobalIndex () != map.getMinAllGlobalIndex () ||
             getMaxAllGlobalIndex () != map.getMaxAllGlobalIndex () ||
             getIndexBase () != map.getIndexBase ()) {
      // If the global min or max global index doesn't match, or if
      // the index base doesn't match, then the Maps aren't the same.
      return false;
    }
    else if (isDistributed () != map.isDistributed ()) {
      // One Map is distributed and the other is not, which means that
      // the Maps aren't the same.
      return false;
    }
    else if (isContiguous () && isUniform () &&
             map.isContiguous () && map.isUniform ()) {
      // Contiguous uniform Maps with the same number of processes in
      // their communicators, with the same global numbers of indices,
      // and with matching index bases and ranges, must be the same.
      return true;
    }

    // The two communicators must have the same number of processes,
    // with process ranks occurring in the same order.  This uses
    // MPI_COMM_COMPARE.  The MPI 3.1 standard (Section 6.4) says:
    // "Operations that access communicators are local and their
    // execution does not require interprocess communication."
    // However, just to be sure, I'll put this call after the above
    // tests that don't communicate.
    if (! Details::congruent (*comm_, * (map.getComm ()))) {
      return false;
    }

    // If we get this far, we need to check local properties and then
    // communicate local sameness across all processes.
    const int isSame_lcl = locallySameAs (map) ? 1 : 0;

    // Return true if and only if all processes report local sameness.
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
    Kokkos::View<const GO*, host_mirror_space> myGlobalInds =
      mapHost_.getMyGlobalIndices (); // creates it if it doesn't exist

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
    using Teuchos::ArrayView;
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

      if (mapDevice_.initialized ()) {
        Kokkos::View<const GO*, DeviceType> myGIDs =
          mapDevice_.getMyGlobalIndices ();
        return rcp (new map_type (globalNumInds, myGIDs, newIndexBase,
                                  newComm, this->getNode ()));
      }
      else {
        Kokkos::View<const GO*, host_mirror_space> myGidsHostView =
          mapHost_.getMyGlobalIndices ();
        ArrayView<const GO> myGidsArrayView (myGidsHostView.ptr_on_device (),
                                             myGidsHostView.dimension_0 ());
        return rcp (new map_type (globalNumInds, myGidsArrayView, newIndexBase,
                                  newComm, this->getNode ()));
      }
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

      // Map's default constructor creates an uninitialized Directory.
      // The Directory will be initialized on demand in
      // getRemoteIndexList().
      //
      // FIXME (mfh 26 Mar 2013) It should be possible to "filter" the
      // directory more efficiently than just recreating it.  If
      // directory recreation proves a bottleneck, we can always
      // revisit this.  On the other hand, Directory creation is only
      // collective over the new, presumably much smaller
      // communicator, so it may not be worth the effort to optimize.

      return map;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  setupDirectory () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      directory_.is_null (), std::logic_error, "Tpetra::Map::setupDirectory: "
      "The Directory is null.  "
      "Please report this bug to the Tpetra developers.");

    // Only create the Directory if it hasn't been created yet.
    // This is a collective operation.
    if (! directory_->initialized ()) {
      directory_->initialize (*this);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LookupStatus
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDs,
                      const Teuchos::ArrayView<int>& PIDs,
                      const Teuchos::ArrayView<LocalOrdinal>& LIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != PIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): GIDs.size ()"
      " = " << GIDs.size () << " != PIDs.size () = " << PIDs.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != LIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (3 args): GIDs.size ()"
      " = " << GIDs.size () << " != LIDs.size () = " << LIDs.size () << ".");

    // Empty Maps (i.e., containing no indices on any processes in the
    // Map's communicator) are perfectly valid.  In that case, if the
    // input GID list is nonempty, we fill the output arrays with
    // invalid values, and return IDNotPresent to notify the caller.
    // It's perfectly valid to give getRemoteIndexList GIDs that the
    // Map doesn't own.  SubmapImport test 2 needs this functionality.
    if (getGlobalNumElements () == 0) {
      if (GIDs.size () == 0) {
        return AllIDsPresent; // trivially
      } else {
        for (Teuchos::ArrayView<int>::size_type k = 0; k < PIDs.size (); ++k) {
          PIDs[k] = Teuchos::OrdinalTraits<int>::invalid ();
        }
        for (typename Teuchos::ArrayView<LocalOrdinal>::size_type k = 0;
             k < LIDs.size (); ++k) {
          LIDs[k] = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
        }
        return IDNotPresent;
      }
    }

    // getRemoteIndexList must be called collectively, and Directory
    // initialization is collective too, so it's OK to initialize the
    // Directory on demand.
    setupDirectory ();
    return directory_->getDirectoryEntries (*this, GIDs, PIDs, LIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  LookupStatus
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDs,
                      const Teuchos::ArrayView<int>& PIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDs.size () != PIDs.size (), std::invalid_argument,
      "Tpetra::Map (Kokkos refactor)::getRemoteIndexList (2 args): GIDs.size ()"
      " = " << GIDs.size () << " != PIDs.size () = " << PIDs.size () << ".");

    // Empty Maps (i.e., containing no indices on any processes in the
    // Map's communicator) are perfectly valid.  In that case, if the
    // input GID list is nonempty, we fill the output array with
    // invalid values, and return IDNotPresent to notify the caller.
    // It's perfectly valid to give getRemoteIndexList GIDs that the
    // Map doesn't own.  SubmapImport test 2 needs this functionality.
    if (getGlobalNumElements () == 0) {
      if (GIDs.size () == 0) {
        return AllIDsPresent; // trivially
      } else {
        // The Map contains no indices, so all output PIDs are invalid.
        for (Teuchos::ArrayView<int>::size_type k = 0; k < PIDs.size (); ++k) {
          PIDs[k] = Teuchos::OrdinalTraits<int>::invalid ();
        }
        return IDNotPresent;
      }
    }

    // getRemoteIndexList must be called collectively, and Directory
    // creation is collective too, so it's OK to create the Directory
    // on demand.
    setupDirectory ();
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

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isOneToOne () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      getComm ().is_null (), std::logic_error, "Tpetra::Map::isOneToOne: "
      "getComm() returns null.  Please report this bug to the Tpetra "
      "developers.");

    // This is a collective operation, if it hasn't been called before.
    setupDirectory ();
    return directory_->isOneToOne (*this);
  }

} // namespace Tpetra

#endif // TPETRA_KOKKOSREFACTOR_MAP_DEF_HPP
