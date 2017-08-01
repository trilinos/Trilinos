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

#ifndef TPETRA_DISTOBJECT_DEF_HPP
#define TPETRA_DISTOBJECT_DEF_HPP

/// \file Tpetra_DistObject_def.hpp
/// \brief Definition of the Tpetra::DistObject class
///
/// If you want to use Tpetra::DistObject, include
/// "Tpetra_DistObject.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::DistObject,
/// include "Tpetra_DistObject_decl.hpp".

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"

namespace Tpetra {

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  DistObject (const Teuchos::RCP<const map_type>& map) :
    map_ (map)
  {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    using Teuchos::RCP;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;

    RCP<Time> doXferTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::doTransfer");
    if (doXferTimer.is_null ()) {
      doXferTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::doTransfer");
    }
    doXferTimer_ = doXferTimer;

    RCP<Time> copyAndPermuteTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::copyAndPermute");
    if (copyAndPermuteTimer.is_null ()) {
      copyAndPermuteTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::copyAndPermute");
    }
    copyAndPermuteTimer_ = copyAndPermuteTimer;

    RCP<Time> packAndPrepareTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::packAndPrepare");
    if (packAndPrepareTimer.is_null ()) {
      packAndPrepareTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::packAndPrepare");
    }
    packAndPrepareTimer_ = packAndPrepareTimer;

    RCP<Time> doPostsAndWaitsTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::doPostsAndWaits");
    if (doPostsAndWaitsTimer.is_null ()) {
      doPostsAndWaitsTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::doPostsAndWaits");
    }
    doPostsAndWaitsTimer_ = doPostsAndWaitsTimer;

    RCP<Time> unpackAndCombineTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::unpackAndCombine");
    if (unpackAndCombineTimer.is_null ()) {
      unpackAndCombineTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::unpackAndCombine");
    }
    unpackAndCombineTimer_ = unpackAndCombineTimer;
#endif // HAVE_TPETRA_TRANSFER_TIMERS
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  DistObject (const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>& rhs)
    : map_ (rhs.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  ~DistObject ()
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  std::string
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "\"Tpetra::DistObject\": {"
       << "Packet: " << TypeNameTraits<packet_type>::name ()
       << ", LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name ()
       << ", GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name ()
       << ", Node: " << TypeNameTraits<Node>::name ();
    if (this->getObjectLabel () != "") {
      os << "Label: \"" << this->getObjectLabel () << "\"";
    }
    os << "}";
    return os.str ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::rcpFromRef;
    using Teuchos::TypeNameTraits;
    using std::endl;
    const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ?
      Teuchos::VERB_LOW : verbLevel;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap ()->getComm ();
    const int myRank = comm.is_null () ? 0 : comm->getRank ();
    const int numProcs = comm.is_null () ? 1 : comm->getSize ();

    if (vl != Teuchos::VERB_NONE) {
      Teuchos::OSTab tab0 (out);
      if (myRank == 0) {
        out << "\"Tpetra::DistObject\":" << endl;
      }
      Teuchos::OSTab tab1 (out);
      if (myRank == 0) {
        out << "Template parameters:" << endl;
        {
          Teuchos::OSTab tab2 (out);
          out << "Packet: " << TypeNameTraits<packet_type>::name () << endl
              << "LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name () << endl
              << "GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name () << endl
              << "Node: " << TypeNameTraits<node_type>::name () << endl;
        }
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
      } // if myRank == 0

      // Describe the Map.
      {
        if (myRank == 0) {
          out << "Map:" << endl;
        }
        Teuchos::OSTab tab2 (out);
        map_->describe (out, vl);
      }

      // At verbosity > VERB_LOW, each process prints something.
      if (vl > Teuchos::VERB_LOW) {
        for (int p = 0; p < numProcs; ++p) {
          if (myRank == p) {
            out << "Process " << myRank << ":" << endl;
            Teuchos::OSTab tab2 (out);
            out << "Export buffer size (in packets): "
                << exports_.dimension_0 ()
                << endl
                << "Import buffer size (in packets): "
                << imports_.dimension_0 ()
                << endl;
          }
          if (! comm.is_null ()) {
            comm->barrier (); // give output time to finish
            comm->barrier ();
            comm->barrier ();
          }
        } // for each process rank p
      } // if vl > VERB_LOW
    } // if vl != VERB_NONE
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Tpetra::DistObject::removeEmptyProcessesInPlace: Not implemented");
  }

  /* These are provided in base DistObject template
  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input,
                               const Teuchos::RCP<const Map<typename DistObjectType::local_ordinal_type,
                                                            typename DistObjectType::global_ordinal_type,
                                                            typename DistObjectType::node_type> >& newMap)
  {
    input->removeEmptyProcessesInPlace (newMap);
    if (newMap.is_null ()) { // my process is excluded
      input = Teuchos::null;
    }
  }

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input)
  {
    using Teuchos::RCP;
    typedef typename DistObjectType::local_ordinal_type LO;
    typedef typename DistObjectType::global_ordinal_type GO;
    typedef typename DistObjectType::node_type NT;
    typedef Map<LO, GO, NT> map_type;

    RCP<const map_type> newMap = input->getMap ()->removeEmptyProcesses ();
    removeEmptyProcessesInPlace<DistObjectType> (input, newMap);
  }
  */

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doImport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
            CombineMode CM)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(*getMap() != *importer.getTargetMap(),
      std::invalid_argument, "doImport: The target DistObject's Map is not "
      "identical to the Import's target Map.");
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *importer.getSourceMap(),
        std::invalid_argument, "doImport: The source is a DistObject, yet its "
        "Map is not identical to the Import's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = importer.getNumSameIDs ();

    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    const view_type exportLIDs      = importer.getExportLIDs();
    const view_type remoteLIDs      = importer.getRemoteLIDs();
    const view_type permuteToLIDs   = importer.getPermuteToLIDs();
    const view_type permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs,
                      remoteLIDs, exportLIDs, importer.getDistributor (),
                      DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doExport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            CombineMode CM)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doExport: The target DistObject's Map is not identical to the Export's "
      "target Map.");
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *exporter.getSourceMap(),
        std::invalid_argument, "doExport: The source is a DistObject, yet its "
        "Map is not identical to the Export's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getExportLIDs();
    view_type remoteLIDs      = exporter.getRemoteLIDs();
    view_type permuteToLIDs   = exporter.getPermuteToLIDs();
    view_type permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doImport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node> & exporter,
            CombineMode CM)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doImport (reverse mode): The target DistObject's Map is not identical "
      "to the Export's source Map.");
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *exporter.getTargetMap(),
        std::invalid_argument,
        "doImport (reverse mode): The source is a DistObject, yet its "
        "Map is not identical to the Export's target Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getRemoteLIDs();
    view_type remoteLIDs      = exporter.getExportLIDs();
    view_type permuteToLIDs   = exporter.getPermuteFromLIDs();
    view_type permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doExport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node> & importer,
            CombineMode CM)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *importer.getSourceMap(), std::invalid_argument,
      "doExport (reverse mode): The target object's Map "
      "is not identical to the Import's source Map.");
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *importer.getTargetMap(),
        std::invalid_argument,
        "doExport (reverse mode): The source is a DistObject, yet its "
        "Map is not identical to the Import's target Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = importer.getNumSameIDs();

    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = importer.getRemoteLIDs();
    view_type remoteLIDs      = importer.getExportLIDs();
    view_type permuteToLIDs   = importer.getPermuteFromLIDs();
    view_type permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, importer.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  isDistributed () const {
    return map_->isDistributed ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  constantNumberOfPackets () const {
    return 0; // default implementation; subclasses may override
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransfer (const SrcDistObject& src,
              CombineMode CM,
              size_t numSameIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs_,
              Distributor& distor,
              ReverseOption revOp)
  {
    using Tpetra::Details::getDualViewCopyFromArrayView;
    typedef LocalOrdinal LO;
    typedef device_type DT;

    if (this->useNewInterface ()) {
      const bool commOnHost = false;

      // Convert arguments to Kokkos::DualView.  This currently
      // involves deep copy, either to host or to device (depending on
      // commOnHost).  At some point, we need to change the interface
      // of doTransfer so it takes DualView (or just View) rather than
      // Teuchos::ArrayView, so that we won't need this deep copy.
      //
      // We don't need to sync the arguments.  commOnHost determines
      // where the most recent version lives.
      Kokkos::DualView<LO*, DT> permuteToLIDs =
        getDualViewCopyFromArrayView<LO, DT> (permuteToLIDs_,
                                              "permuteToLIDs",
                                              commOnHost);
      Kokkos::DualView<LO*, DT> permuteFromLIDs =
        getDualViewCopyFromArrayView<LO, DT> (permuteFromLIDs_,
                                              "permuteFromLIDs",
                                              commOnHost);
      // No need to sync this.  packAndPrepareNew will use it to
      // determine where to pack (in host or device memory).
      Kokkos::DualView<LO*, DT> remoteLIDs =
        getDualViewCopyFromArrayView<LO, DT> (remoteLIDs_,
                                              "remoteLIDs",
                                              commOnHost);
      Kokkos::DualView<LO*, DT> exportLIDs =
        getDualViewCopyFromArrayView<LO, DT> (exportLIDs_,
                                              "exportLIDs",
                                              commOnHost);

      doTransferNew (src, CM, numSameIDs, permuteToLIDs, permuteFromLIDs,
                     remoteLIDs, exportLIDs, distor, revOp, commOnHost);
    }
    else {
      doTransferOld (src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_,
                     remoteLIDs_, exportLIDs_, distor, revOp);
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  reallocImportsIfNeeded (const size_t newSize, const bool debug)
  {
    if (debug) {
      std::ostringstream os;
      os << "*** Reallocate (if needed) imports_ from "
         << imports_.dimension_0 () << " to " << newSize << std::endl;
      std::cerr << os.str ();
    }
    using Details::reallocDualViewIfNeeded;
    return reallocDualViewIfNeeded (this->imports_, newSize, "imports");
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  reallocArraysForNumPacketsPerLid (const size_t numExportLIDs,
                                    const size_t numImportLIDs)
  {
    using Details::reallocDualViewIfNeeded;

    // If an array is already allocated, and if is at least
    // tooBigFactor times bigger than it needs to be, free it and
    // reallocate to the size we need, in order to save space.
    // Otherwise, take subviews to reduce allocation size.
    constexpr size_t tooBigFactor = 10;

    // Reallocate numExportPacketsPerLID_ if needed.
    const bool firstReallocated =
      reallocDualViewIfNeeded (this->numExportPacketsPerLID_,
                               numExportLIDs,
                               "numExportPacketsPerLID",
                               tooBigFactor,
                               true); // need fence before, if realloc'ing

    // If we reallocated above, then we fenced after that
    // reallocation.  This means that we don't need to fence again,
    // before the next reallocation.
    const bool needFenceBeforeNextAlloc = ! firstReallocated;
    const bool secondReallocated =
      reallocDualViewIfNeeded (this->numImportPacketsPerLID_,
                               numImportLIDs,
                               "numExportPacketsPerLID",
                               tooBigFactor,
                               needFenceBeforeNextAlloc);
    return firstReallocated || secondReallocated;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransferOld (const SrcDistObject& src,
              CombineMode CM,
              size_t numSameIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
              Distributor &distor,
              ReverseOption revOp)
  {
    using Tpetra::Details::getArrayViewFromDualView;
    using Tpetra::Details::reallocDualViewIfNeeded;
    const bool debug = false;

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! checkSizes (src), std::invalid_argument,
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");
    KokkosClassic::ReadWriteOption rwo = KokkosClassic::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      const size_t numIDsToWrite = numSameIDs +
        static_cast<size_t> (permuteToLIDs.size ()) +
        static_cast<size_t> (remoteLIDs.size ());
      if (numIDsToWrite == this->getMap ()->getNodeNumElements ()) {
        // We're overwriting all of our local data in the destination
        // object, so a write-only view suffices.
        //
        // FIXME (mfh 10 Apr 2012) This doesn't make sense for a
        // CrsMatrix with a dynamic graph.  INSERT mode could mean
        // that we're adding new entries to the object, but we don't
        // want to get rid of the old ones.
        rwo = KokkosClassic::WriteOnly;
      }
    }
    // Tell the source to create a read-only view of its data.  On a
    // discrete accelerator such as a GPU, this brings EVERYTHING from
    // device memory to host memory.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs to send) and packet counts, createViews()
    // could create a "sparse view" that only brings in the necessary
    // data from device to host memory.
    const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);
    if (srcDistObj != NULL) {
      srcDistObj->createViews ();
    }

    // Tell the target to create a view of its data.  Depending on
    // rwo, this could be a write-only view or a read-and-write view.
    // On a discrete accelerator such as a GPU, a write-only view only
    // requires a transfer from host to device memory.  A
    // read-and-write view requires a two-way transfer.  This has the
    // same problem as createViews(): it transfers EVERYTHING, not
    // just the necessary data.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs into which to receive) and packet counts,
    // createViewsNonConst() could create a "sparse view" that only
    // transfers the necessary data.
    this->createViewsNonConst (rwo);

    if (numSameIDs + permuteToLIDs.size()) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
      // There is at least one GID to copy or permute.
      copyAndPermute (src, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }

    // The method may return zero even if the implementation actually
    // does have a constant number of packets per LID.  However, if it
    // returns nonzero, we may use this information to avoid
    // (re)allocating num{Ex,Im}portPacketsPerLID_.  packAndPrepare()
    // will set this to its final value.
    //
    // We only need this if CM != ZERO, but it has to be lifted out of
    // that scope because there are multiple tests for CM != ZERO.
    size_t constantNumPackets = this->constantNumberOfPackets ();

    // We only need to pack communication buffers if the combine mode
    // is not ZERO. A "ZERO combine mode" means that the results are
    // the same as if we had received all zeros, and added them to the
    // existing values. That means we don't need to communicate.
    if (CM != ZERO) {
      if (constantNumPackets == 0) {
        this->reallocArraysForNumPacketsPerLid (exportLIDs.size (),
                                                remoteLIDs.size ());
      }

      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        // Ask the source to pack data.  Also ask it whether there are a
        // constant number of packets per element (constantNumPackets is
        // an output argument).  If there are, constantNumPackets will
        // come back nonzero.  Otherwise, the source will fill the
        // numExportPacketsPerLID_ array.
        numExportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
        Teuchos::ArrayView<size_t> numExportPacketsPerLID =
          getArrayViewFromDualView (numExportPacketsPerLID_);

        // FIXME (mfh 26 Apr 2016) For backwards compatibility, use
        // the old packAndPrepare interface that takes and resizes the
        // exports buffer as a Teuchos::Array<packet_type>.  Then,
        // copy out that buffer into the host version of exports_.

        Teuchos::Array<packet_type> exportsOld;
        packAndPrepare (src, exportLIDs, exportsOld, numExportPacketsPerLID,
                        constantNumPackets, distor);
        const size_t exportsLen = static_cast<size_t> (exportsOld.size ());
        reallocDualViewIfNeeded (this->exports_, exportsLen, "exports");
        Kokkos::View<const packet_type*, Kokkos::HostSpace,
          Kokkos::MemoryUnmanaged> exportsOldK (exportsOld.getRawPtr (),
                                                exportsLen);
        exports_.template modify<Kokkos::HostSpace> ();
        Kokkos::deep_copy (exports_.template view<Kokkos::HostSpace> (),
                           exportsOldK);
      }
    }

    // We don't need the source's data anymore, so it can let go of
    // its views.  On an accelerator device with a separate memory
    // space (like a GPU), this frees host memory, since device memory
    // has the "master" version of the data.
    if (srcDistObj != NULL) {
      srcDistObj->releaseViews ();
    }

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.size() * constantNumPackets;
        if (debug) {
          std::ostringstream os;
          os << "*** doTransferOld: Const # packets: imports_.dimension_0() = "
             << imports_.dimension_0 () << ", rbufLen = " << rbufLen
             << std::endl;
          std::cerr << os.str ();
        }
        reallocImportsIfNeeded (rbufLen, debug);
      }

      // Do we need to do communication (via doPostsAndWaits)?
      bool needCommunication = true;
      if (revOp == DoReverse && ! isDistributed ()) {
        needCommunication = false;
      }
      // FIXME (mfh 30 Jun 2013): Checking whether the source object
      // is distributed requires a cast to DistObject.  If it's not a
      // DistObject, then I'm not quite sure what to do.  Perhaps it
      // would be more appropriate for SrcDistObject to have an
      // isDistributed() method.  For now, I'll just assume that we
      // need to do communication unless the cast succeeds and the
      // source is not distributed.
      else if (revOp == DoForward && srcDistObj != NULL &&
               ! srcDistObj->isDistributed ()) {
        needCommunication = false;
      }

      if (needCommunication) {
        if (revOp == DoReverse) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            // First communicate the number of packets per LID to receive.

            // Make sure that host has the latest version, since we're
            // using the version on host.  If host has the latest
            // version already, syncing to host does nothing.
            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
              getArrayViewFromDualView (numExportPacketsPerLID_);

            // numImportPacketsPerLID_ is the output array here, so
            // mark it as modified.  It's strictly output, so we don't
            // have to sync from device.
            //numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<size_t> numImportPacketsPerLID =
              getArrayViewFromDualView (numImportPacketsPerLID_);
            distor.doReversePostsAndWaits (numExportPacketsPerLID, 1,
                                           numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            {
              typedef typename Kokkos::DualView<size_t*,
                device_type>::t_host::execution_space host_exec_space;
              typedef Kokkos::RangePolicy<host_exec_space, Array_size_type> range_type;
              const size_t* const arrayToSum = numImportPacketsPerLID.getRawPtr ();
              Kokkos::parallel_reduce ("Count import packets",
                                       range_type (0, numImportPacketsPerLID.size ()),
                                       [=] (const Array_size_type& i, size_t& lclSum) {
                                         lclSum += arrayToSum[i];
                                       }, totalImportPackets);
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doReversePostsAndWaits (hostExports,
                                           numExportPacketsPerLID,
                                           hostImports,
                                           numImportPacketsPerLID);
          }
          else {
            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doReversePostsAndWaits (hostExports,
                                           constantNumPackets,
                                           hostImports);
          }
        }
        else { // revOp == DoForward
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            // First communicate the number of packets per LID to receive.

            // Make sure that host has the latest version, since we're
            // using the version on host.  If host has the latest
            // version already, syncing to host does nothing.
            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
              getArrayViewFromDualView (numExportPacketsPerLID_);

            // numImportPacketsPerLID_ is the output array here, so
            // mark it as modified.  It's strictly output, so we don't
            // have to sync from device.
            //numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<size_t> numImportPacketsPerLID =
              getArrayViewFromDualView (numImportPacketsPerLID_);
            distor.doPostsAndWaits (numExportPacketsPerLID, 1,
                                    numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            {
              typedef typename Kokkos::DualView<size_t*,
                device_type>::t_host::execution_space host_exec_space;
              typedef Kokkos::RangePolicy<host_exec_space, Array_size_type> range_type;
              const size_t* const arrayToSum = numImportPacketsPerLID.getRawPtr ();
              Kokkos::parallel_reduce ("Count import packets",
                                       range_type (0, numImportPacketsPerLID.size ()),
                                       [=] (const Array_size_type& i, size_t& lclSum) {
                                         lclSum += arrayToSum[i];
                                       }, totalImportPackets);
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doPostsAndWaits (hostExports,
                                    numExportPacketsPerLID,
                                    hostImports,
                                    numImportPacketsPerLID);
          }
          else {
            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doPostsAndWaits (hostExports,
                                    constantNumPackets,
                                    hostImports);
          }
        }
        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

          // We don't need to sync imports_, because it is only for
          // output here.  This legacy version of doTransfer only uses
          // host arrays.
          imports_.template modify<Kokkos::HostSpace> ();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView (imports_);
          // NOTE (mfh 25 Apr 2016) unpackAndCombine doesn't actually
          // change its numImportPacketsPerLID argument, so we don't
          // have to mark it modified here.
          numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
          // FIXME (mfh 25 Apr 2016) unpackAndCombine doesn't actually
          // change its numImportPacketsPerLID argument, so we should
          // be able to use a const Teuchos::ArrayView here.
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView (numImportPacketsPerLID_);
          unpackAndCombine (remoteLIDs, hostImports, numImportPacketsPerLID,
                            constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    this->releaseViews ();
  }

  namespace { // (anonymous)
    template<class DeviceType, class IndexType = size_t>
    struct SumFunctor {
      SumFunctor (const Kokkos::View<const size_t*, DeviceType>& viewToSum) :
        viewToSum_ (viewToSum) {}
      KOKKOS_FUNCTION void operator() (const IndexType& i, size_t& lclSum) const {
        lclSum += viewToSum_(i);
      }
      Kokkos::View<const size_t*, DeviceType> viewToSum_;
    };

    template<class DeviceType, class IndexType = size_t>
    size_t
    countTotalImportPackets (const Kokkos::View<const size_t*, DeviceType>& numImportPacketsPerLID)
    {
      using Kokkos::parallel_reduce;
      typedef DeviceType DT;
      typedef typename DT::execution_space DES;
      typedef Kokkos::RangePolicy<DES, IndexType> range_type;

      const IndexType numOut = numImportPacketsPerLID.dimension_0 ();
      size_t totalImportPackets = 0;
      parallel_reduce ("Count import packets",
                       range_type (0, numOut),
                       SumFunctor<DeviceType, IndexType> (numImportPacketsPerLID),
                       totalImportPackets);
      return totalImportPackets;
    }
  } // namespace (anonymous)

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransferNew (const SrcDistObject& src,
                 const CombineMode CM,
                 const size_t numSameIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& permuteToLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& permuteFromLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& remoteLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& exportLIDs,
                 Distributor& distor,
                 const ReverseOption revOp,
                 const bool commOnHost)
  {
    using Tpetra::Details::getArrayViewFromDualView;
    using Kokkos::Compat::getArrayView;
    using Kokkos::Compat::getConstArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::Compat::create_const_view;
    typedef LocalOrdinal LO;
    typedef device_type DT;
    typedef typename Kokkos::DualView<LO*, DT>::t_dev::execution_space DES;
    typedef typename Kokkos::DualView<LO*, DT>::t_dev::memory_space DMS;
    //typedef typename Kokkos::DualView<LO*, DT>::t_dev::memory_space HMS;
    typedef Kokkos::HostSpace HMS; // prevent DualView with CudaUVMSpace issues
    const bool debug = false;

    if (debug) {
      std::ostringstream os;
      os << ">>> DistObject::doTransferNew: remoteLIDs.size() = "
         << remoteLIDs.dimension_0 () << std::endl;
      std::cerr << os.str ();
    }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    if (debug) {
      std::cerr << ">>> 1. checkSizes" << std::endl;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! checkSizes (src), std::invalid_argument,
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");

    // NOTE (mfh 26 Apr 2016) Chris Baker's implementation understood
    // that if CM == INSERT || CM == REPLACE, the target object could
    // be write only.  We don't optimize for that here.

    if (debug) {
      std::cerr << ">>> 2. copyAndPermuteNew" << std::endl;
    }

    if (numSameIDs + permuteToLIDs.dimension_0 () != 0) {
      // There is at least one GID to copy or permute.
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
      copyAndPermuteNew (src, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }

    // The method may return zero even if the implementation actually
    // does have a constant number of packets per LID.  However, if it
    // returns nonzero, we may use this information to avoid
    // (re)allocating num{Ex,Im}portPacketsPerLID_.  packAndPrepare()
    // will set this to its final value.
    //
    // We only need this if CM != ZERO, but it has to be lifted out of
    // that scope because there are multiple tests for CM != ZERO.
    size_t constantNumPackets = this->constantNumberOfPackets ();

    // We only need to pack communication buffers if the combine mode
    // is not ZERO. A "ZERO combine mode" means that the results are
    // the same as if we had received all zeros, and added them to the
    // existing values. That means we don't need to communicate.
    if (CM != ZERO) {
      if (constantNumPackets == 0) {
        if (debug) {
          std::cerr << ">>> 3. Allocate num{Ex,Im}portPacketsPerLID" << std::endl;
        }
        // This only reallocates if necessary, that is, if the sizes
        // don't match.
        this->reallocArraysForNumPacketsPerLid (exportLIDs.dimension_0 (),
                                                remoteLIDs.dimension_0 ());
      }

      if (debug) {
        std::cerr << ">>> 4. packAndPrepareNew" << std::endl;
      }

      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

        if (debug) {
          std::ostringstream os;
          const int myRank = this->getMap ()->getComm ()->getRank ();
          os << ">>> (Proc " << myRank << ") 5.0. Before packAndPrepareNew, "
            "exports_.dimension_0()=" << exports_.dimension_0 () << std::endl;
          std::cerr << os.str ();
        }
        // Ask the source to pack data.  Also ask it whether there are
        // a constant number of packets per element
        // (constantNumPackets is an output argument).  If there are,
        // constantNumPackets will come back nonzero.  Otherwise, the
        // source will fill the numExportPacketsPerLID_ array.
        packAndPrepareNew (src, exportLIDs, exports_, numExportPacketsPerLID_,
                           constantNumPackets, distor);
        if (debug) {
          std::ostringstream os;
          const int myRank = this->getMap ()->getComm ()->getRank ();
          os << ">>> (Proc " << myRank << ") 5.0. After packAndPrepareNew, "
            "exports_.dimension_0()=" << exports_.dimension_0 () << std::endl;
          std::cerr << os.str ();
        }
      }
    }

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        if (debug) {
          std::cerr << ">>> 6. Realloc imports_" << std::endl;
        }
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.dimension_0 () * constantNumPackets;
        reallocImportsIfNeeded (rbufLen, debug);
      }

      // Do we need to do communication (via doPostsAndWaits)?
      bool needCommunication = true;

      // This may be NULL.  It will be used below.
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);

      if (revOp == DoReverse && ! isDistributed ()) {
        needCommunication = false;
      }
      // FIXME (mfh 30 Jun 2013): Checking whether the source object
      // is distributed requires a cast to DistObject.  If it's not a
      // DistObject, then I'm not quite sure what to do.  Perhaps it
      // would be more appropriate for SrcDistObject to have an
      // isDistributed() method.  For now, I'll just assume that we
      // need to do communication unless the cast succeeds and the
      // source is not distributed.
      else if (revOp == DoForward && srcDistObj != NULL &&
               ! srcDistObj->isDistributed ()) {
        needCommunication = false;
      }

      // FIXME (mfh 17 Feb 2014) Distributor doesn't actually inspect
      // the contents of the "exports" or "imports" arrays, other than
      // to do a deep copy in the (should be technically unnecessary,
      // but isn't for some odd reason) "self-message" case.
      // Distributor thus doesn't need host views; it could do just
      // fine with device views, assuming that MPI knows how to read
      // device memory (which doesn't even require UVM).

      if (needCommunication) {
        if (revOp == DoReverse) {
          if (debug) {
            std::cerr << ">>> 7.0. Reverse mode" << std::endl;
          }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (debug) {
              std::cerr << ">>> 7.1. Variable # packets / LID: first comm"
                        << std::endl;
            }

            size_t totalImportPackets = 0;
            if (commOnHost) {
              this->numExportPacketsPerLID_.template sync<HMS> ();
              this->numImportPacketsPerLID_.template sync<HMS> ();
              this->numImportPacketsPerLID_.template modify<HMS> (); // output argument
              auto numExp_h = create_const_view (this->numExportPacketsPerLID_.template view<HMS> ());
              auto numImp_h = this->numImportPacketsPerLID_.template view<HMS> ();

              // MPI communication happens here.
              distor.doReversePostsAndWaits (numExp_h, 1, numImp_h);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_h)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_h);
            }
            else {
              this->numExportPacketsPerLID_.template sync<DMS> ();
              this->numImportPacketsPerLID_.template sync<DMS> ();
              this->numImportPacketsPerLID_.template modify<DMS> (); // output argument
              auto numExp_d = create_const_view (this->numExportPacketsPerLID_.template view<DMS> ());
              auto numImp_d = this->numImportPacketsPerLID_.template view<DMS> ();

              // MPI communication happens here.
              distor.doReversePostsAndWaits (numExp_d, 1, numImp_d);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_d)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_d);
            }

            this->reallocImportsIfNeeded (totalImportPackets, debug);

            if (debug) {
              std::cerr << ">>> 7.3. Second comm" << std::endl;
            }

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) Since we need to
            // launch MPI communication on host, we will need
            // numExportPacketsPerLID and numImportPacketsPerLID on
            // host.
            this->numExportPacketsPerLID_.template sync<HMS> ();
            this->numImportPacketsPerLID_.template sync<HMS> ();

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) doPostsAndWaits and
            // doReversePostsAndWaits currently want
            // numExportPacketsPerLID and numImportPacketsPerLID as
            // Teuchos::ArrayView, rather than as Kokkos::View.
            auto numExportPacketsPerLID_av =
              getArrayViewFromDualView (this->numExportPacketsPerLID_);
            auto numImportPacketsPerLID_av =
              getArrayViewFromDualView (this->numImportPacketsPerLID_);

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<HMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<HMS> ()),
                                             numExportPacketsPerLID_av,
                                             this->imports_.template view<HMS> (),
                                             numImportPacketsPerLID_av);
            }
            else {
              this->imports_.template modify<DMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<DMS> ()),
                                             numExportPacketsPerLID_av,
                                             this->imports_.template view<DMS> (),
                                             numImportPacketsPerLID_av);
            }
          }
          else {
            if (debug) {
              const int myRank = this->getMap ()->getComm ()->getRank ();
              std::ostringstream os;
              os << ">>> (Proc " << myRank << "): 7.1. Const # packets per "
                "LID: exports_.dimension_0() = " << exports_.dimension_0 ()
                 << ", imports_.dimension_0() = " << imports_.dimension_0 ()
                 << std::endl;
              std::cerr << os.str ();
            }

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<HMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<HMS> ()),
                                             constantNumPackets,
                                             this->imports_.template view<HMS> ());
            }
            else { // pack on device
              this->imports_.template modify<DMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<DMS> ()),
                                             constantNumPackets,
                                             this->imports_.template view<DMS> ());
            }
          }
        }
        else { // revOp == DoForward
          if (debug) {
            std::cerr << ">>> 7.0. Forward mode" << std::endl;
          }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (debug) {
              std::cerr << ">>> 7.1. Variable # packets / LID: first comm" << std::endl;
            }

            size_t totalImportPackets = 0;
            if (commOnHost) {
              this->numExportPacketsPerLID_.template sync<HMS> ();
              this->numImportPacketsPerLID_.template sync<HMS> ();
              this->numImportPacketsPerLID_.template modify<HMS> (); // output argument
              auto numExp_h = create_const_view (this->numExportPacketsPerLID_.template view<HMS> ());
              auto numImp_h = this->numImportPacketsPerLID_.template view<HMS> ();

              // MPI communication happens here.
              distor.doPostsAndWaits (numExp_h, 1, numImp_h);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_h)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_h);
            }
            else {
              this->numExportPacketsPerLID_.template sync<DMS> ();
              this->numImportPacketsPerLID_.template sync<DMS> ();
              this->numImportPacketsPerLID_.template modify<DMS> (); // output argument
              auto numExp_d = create_const_view (this->numExportPacketsPerLID_.template view<DMS> ());
              auto numImp_d = this->numImportPacketsPerLID_.template view<DMS> ();

              // MPI communication happens here.
              distor.doPostsAndWaits (numExp_d, 1, numImp_d);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_d)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_d);
            }

            this->reallocImportsIfNeeded (totalImportPackets, debug);

            if (debug) {
              std::cerr << ">>> 7.3. Second comm" << std::endl;
            }

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) Since we need to
            // launch MPI communication on host, we will need
            // numExportPacketsPerLID and numImportPacketsPerLID on
            // host.
            this->numExportPacketsPerLID_.template sync<HMS> ();
            this->numImportPacketsPerLID_.template sync<HMS> ();

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) doPostsAndWaits and
            // doReversePostsAndWaits currently want
            // numExportPacketsPerLID and numImportPacketsPerLID as
            // Teuchos::ArrayView, rather than as Kokkos::View.
            auto numExportPacketsPerLID_av =
              getArrayViewFromDualView (this->numExportPacketsPerLID_);
            auto numImportPacketsPerLID_av =
              getArrayViewFromDualView (this->numImportPacketsPerLID_);

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<HMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<HMS> ()),
                                      numExportPacketsPerLID_av,
                                      this->imports_.template view<HMS> (),
                                      numImportPacketsPerLID_av);
            }
            else { // pack on device
              this->imports_.template modify<DMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<DMS> ()),
                                      numExportPacketsPerLID_av,
                                      this->imports_.template view<DMS> (),
                                      numImportPacketsPerLID_av);
            }
          }
          else {
            if (debug) {
              const int myRank = this->getMap ()->getComm ()->getRank ();
              std::ostringstream os;
              os << ">>> (Proc " << myRank << "): 7.1. Const # packets per "
                "LID: exports_.dimension_0()=" << exports_.dimension_0 ()
                 << ", imports_.dimension_0() = " << imports_.dimension_0 ()
                 << std::endl;
              std::cerr << os.str ();
            }

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<HMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<HMS> ()),
                                      constantNumPackets,
                                      this->imports_.template view<HMS> ());
            }
            else { // pack on device
              this->imports_.template modify<DMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<DMS> ()),
                                      constantNumPackets,
                                      this->imports_.template view<DMS> ());
            }
          }
        }

        if (debug) {
          std::cerr << ">>> 8. unpackAndCombineNew" << std::endl;
        }

        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

          // NOTE (mfh 26 Apr 2016) We don't actually need to sync the
          // input DualViews, but they DO need to be most recently
          // updated in the same memory space.
          //
          // FIXME (mfh 26 Apr 2016) Check that all input DualViews
          // were most recently updated in the same memory space, and
          // sync them to the same place (based on commOnHost) if not.
          unpackAndCombineNew (remoteLIDs, imports_, numImportPacketsPerLID_,
                               constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    if (debug) {
      std::cerr << ">>> 9. Done with doTransferNew" << std::endl;
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  print (std::ostream &os) const
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using std::endl;

    RCP<FancyOStream> out = getFancyOStream (rcpFromRef (os));
    this->describe (*out, Teuchos::VERB_DEFAULT);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  createViews () const
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  createViewsNonConst (KokkosClassic::ReadWriteOption /*rwo*/)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  releaseViews () const
  {}

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input,
                               const Teuchos::RCP<const Map<typename DistObjectType::local_ordinal_type,
                                                            typename DistObjectType::global_ordinal_type,
                                                            typename DistObjectType::node_type> >& newMap)
  {
    input->removeEmptyProcessesInPlace (newMap);
    if (newMap.is_null ()) { // my process is excluded
      input = Teuchos::null;
    }
  }

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input)
  {
    using Teuchos::RCP;
    typedef typename DistObjectType::local_ordinal_type LO;
    typedef typename DistObjectType::global_ordinal_type GO;
    typedef typename DistObjectType::node_type NT;
    typedef Map<LO, GO, NT> map_type;

    RCP<const map_type> newMap = input->getMap ()->removeEmptyProcesses ();
    removeEmptyProcessesInPlace<DistObjectType> (input, newMap);
  }

// Explicit instantiation macro for general DistObject.
#define TPETRA_DISTOBJECT_INSTANT(SCALAR, LO, GO, NODE) \
  template class DistObject< SCALAR , LO , GO , NODE >;

// Explicit instantiation macro for DistObject<char, ...>.
// The "SLGN" stuff above doesn't work for Packet=char.
#define TPETRA_DISTOBJECT_INSTANT_CHAR(LO, GO, NODE) \
  template class DistObject< char , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_DISTOBJECT_DEF_HPP
