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
            out << "Export buffer size (in packets): " << exports_.size ()
                << endl
                << "Import buffer size (in packets): " << imports_.size ()
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
    TEUCHOS_TEST_FOR_EXCEPTION(*getMap() != *importer.getTargetMap(),
      std::invalid_argument, "doImport: The target DistObject's Map is not "
      "identical to the Import's target Map.");
#ifdef HAVE_TPETRA_DEBUG
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doExport: The target DistObject's Map is not identical to the Export's "
      "target Map.");
#ifdef HAVE_TPETRA_DEBUG
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *exporter.getSourceMap(),
        std::invalid_argument, "doExport: The source is a DistObject, yet its "
        "Map is not identical to the Export's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doImport (reverse mode): The target DistObject's Map is not identical "
      "to the Export's source Map.");
#ifdef HAVE_TPETRA_DEBUG
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

    typedef ArrayView<const LocalOrdinal> view_type;
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *importer.getSourceMap(), std::invalid_argument,
      "doExport (reverse mode): The target object's Map "
      "is not identical to the Import's source Map.");
#ifdef HAVE_TPETRA_DEBUG
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

    typedef ArrayView<const LocalOrdinal> view_type;
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
              Distributor &distor,
              ReverseOption revOp)
  {
    if (this->useNewInterface()) {
      doTransferNew(src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, exportLIDs_, distor, revOp);
    }
    else {
      doTransferOld(src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, exportLIDs_, distor, revOp);
    }
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
        numExportPacketsPerLID_ =
          decltype (numExportPacketsPerLID_) ("numExportPacketsPerLID",
                                              exportLIDs.size ());
        host_numExportPacketsPerLID_ =
          Kokkos::create_mirror_view (numExportPacketsPerLID_);
        numImportPacketsPerLID_ =
          decltype (numImportPacketsPerLID_) ("numImportPacketsPerLID",
                                              remoteLIDs.size ());
        host_numImportPacketsPerLID_ =
          Kokkos::create_mirror_view (numImportPacketsPerLID_);
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
        Teuchos::ArrayView<size_t> numExportPacketsPerLID (host_numExportPacketsPerLID_.ptr_on_device (), host_numExportPacketsPerLID_.dimension_0 ());
        packAndPrepare (src, exportLIDs, exports_old_, numExportPacketsPerLID,
                        constantNumPackets, distor);
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

        if (static_cast<size_t> (imports_.dimension_0 ()) != rbufLen ||
            host_imports_.dimension_0 () != imports_.dimension_0 ()) {
          Kokkos::realloc (imports_, rbufLen);
          // This is doTransferOld, so we need the host version of
          // imports_.  We'll wrap its data in a Teuchos::ArrayView
          // and use the backwards compatibility interface.
          host_imports_ = Kokkos::create_mirror_view (imports_);
        }
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
            Teuchos::ArrayView<size_t> numExportPacketsPerLID (host_numExportPacketsPerLID_.ptr_on_device (), host_numExportPacketsPerLID_.dimension_0 ());
            Teuchos::ArrayView<size_t> numImportPacketsPerLID (host_numImportPacketsPerLID_.ptr_on_device (), host_numImportPacketsPerLID_.dimension_0 ());
            distor.doReversePostsAndWaits (numExportPacketsPerLID.getConst (),
                                           1,
                                           numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
              totalImportPackets += numImportPacketsPerLID[i];
            }

            if (static_cast<size_t> (imports_.dimension_0 ()) != totalImportPackets ||
                host_imports_.dimension_0 () != imports_.dimension_0 ()) {
              Kokkos::realloc (imports_, totalImportPackets);
              // This is doTransferOld, so we need the host version of
              // imports_.  We'll wrap its data in a Teuchos::ArrayView
              // and use the backwards compatibility interface.
              host_imports_ = Kokkos::create_mirror_view (imports_);
            }
            Teuchos::ArrayView<packet_type> hostImports (host_imports_.ptr_on_device (),
                                                         host_imports_.dimension_0 ());
            distor.doReversePostsAndWaits (exports_old_().getConst(),
                                           numExportPacketsPerLID,
                                           hostImports,
                                           numImportPacketsPerLID);
          }
          else {
            Teuchos::ArrayView<packet_type> hostImports (host_imports_.ptr_on_device (),
                                                         host_imports_.dimension_0 ());
            distor.doReversePostsAndWaits (exports_old_ ().getConst (),
                                           constantNumPackets,
                                           hostImports);
          }
        }
        else { // revOp == DoForward
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            Teuchos::ArrayView<size_t> numExportPacketsPerLID (host_numExportPacketsPerLID_.ptr_on_device (), host_numExportPacketsPerLID_.dimension_0 ());
            Teuchos::ArrayView<size_t> numImportPacketsPerLID (host_numImportPacketsPerLID_.ptr_on_device (), host_numImportPacketsPerLID_.dimension_0 ());
            distor.doPostsAndWaits (numExportPacketsPerLID.getConst (), 1,
                                    numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
              totalImportPackets += numImportPacketsPerLID[i];
            }

            if (static_cast<size_t> (imports_.dimension_0 ()) != totalImportPackets ||
                host_imports_.dimension_0 () != imports_.dimension_0 ()) {
              Kokkos::realloc (imports_, totalImportPackets);
              // This is doTransferOld, so we need the host version of
              // imports_.  We'll wrap its data in a Teuchos::ArrayView
              // and use the backwards compatibility interface.
              host_imports_ = Kokkos::create_mirror_view (imports_);
            }
            Teuchos::ArrayView<packet_type> hostImports (host_imports_.ptr_on_device (),
                                                         host_imports_.dimension_0 ());
            distor.doPostsAndWaits (exports_old_().getConst(),
                                    numExportPacketsPerLID,
                                    hostImports,
                                    numImportPacketsPerLID);
          }
          else {
            Teuchos::ArrayView<packet_type> hostImports (host_imports_.ptr_on_device (),
                                                         host_imports_.dimension_0 ());
            distor.doPostsAndWaits (exports_old_ ().getConst (),
                                    constantNumPackets,
                                    hostImports);
          }
        }
        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::ArrayView<size_t> numImportPacketsPerLID (host_numImportPacketsPerLID_.ptr_on_device (), host_numImportPacketsPerLID_.dimension_0 ());
          Teuchos::ArrayView<packet_type> hostImports (host_imports_.ptr_on_device (),
                                                       host_imports_.dimension_0 ());
          unpackAndCombine (remoteLIDs, hostImports, numImportPacketsPerLID,
                            constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    this->releaseViews ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransferNew (const SrcDistObject& src,
                 CombineMode CM,
                 size_t numSameIDs,
                 const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs_,
                 Distributor &distor,
                 ReverseOption revOp)
  {
    using Kokkos::Compat::getArrayView;
    using Kokkos::Compat::getConstArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::Compat::create_const_view;

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    // Convert arguments to Kokkos::View's (involves deep copy to device)
    //
    // FIXME (mfh 17 Feb 2014) getKokkosViewDeepCopy _always_ does a
    // deep copy.  It has to, since it returns a managed Kokkos::View,
    // but the input Teuchos::ArrayView is unmanaged.  One way to fix
    // this would be to change the interface by replacing the
    // Teuchos::ArrayView inputs with Kokkos::View inputs.  This is
    // the approach taken by Kokkos::DistObjectKA.  Some points:
    //
    //   1. It's perfectly OK to change the interface of doTransfer.
    //      It should take Kokkos::View by default, and convert to
    //      Teuchos::Array{RCP,View} if needed internally.
    //   2. Recall that Teuchos::ArrayView is an unmanaged view.
    //   3. If DistObject ever gets a nonblocking interface, that
    //      interface should take managed views.
    typedef Kokkos::View<const LocalOrdinal*, execution_space> lo_const_view_type;
    lo_const_view_type permuteToLIDs =
      getKokkosViewDeepCopy<execution_space> (permuteToLIDs_);
    lo_const_view_type permuteFromLIDs =
      getKokkosViewDeepCopy<execution_space> (permuteFromLIDs_);
    lo_const_view_type remoteLIDs =
      getKokkosViewDeepCopy<execution_space> (remoteLIDs_);
    lo_const_view_type exportLIDs =
      getKokkosViewDeepCopy<execution_space> (exportLIDs_);

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

    // FIXME (mfh 17 Feb 2014) We're assuming that MPI can read device
    // memory (that's even pre-UVM), so there is no need to create
    // host views of the source object's data.

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

    // FIXME (mfh 17 Feb 2014) We're assuming that MPI can read device
    // memory (that's even pre-UVM), so there is no need to create
    // host views of the target object's data.

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

      // FIXME (mfh 17 Feb 2014) Nobody implements DistObject
      // subclasses but Tpetra developers anyway, so don't bother with
      // this "new" business.  Just write the interface you want.

      // There is at least one GID to copy or permute.
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

        // FIXME (mfh 17 Feb 2014) Don't "realloc" unless you really
        // need to.  That will avoid a bit of time for reinitializing
        // the Views' data.
        Kokkos::realloc (numExportPacketsPerLID_, exportLIDs.size ());
        host_numExportPacketsPerLID_ = Kokkos::create_mirror_view (numExportPacketsPerLID_);
        Kokkos::realloc (numImportPacketsPerLID_, remoteLIDs.size ());
        host_numImportPacketsPerLID_ = Kokkos::create_mirror_view (numImportPacketsPerLID_);
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
        packAndPrepareNew (src, exportLIDs, exports_, numExportPacketsPerLID_,
                        constantNumPackets, distor);
      }
    }

    // FIXME (mfh 17 Feb 2014) See comments above; there is no need to
    // create host views of the source object's data.

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
        if (static_cast<size_t> (imports_.dimension_0 ()) != rbufLen) {
          // FIXME (mfh 28 Mar 2016) This helps fix #227.
          execution_space::fence ();

          Kokkos::realloc (imports_, rbufLen);
          // This is doTransferNew, so we don't need to allocate the
          // host version of imports_.
        }
      }

      // FIXME (mfh 17 Feb 2014) Why do we need mirror views?
      // Furthermore, if we do need to do this, we should only do it
      // _once_, since the arrays are constant (they come from the
      // Import / Export object, which is constant).

      TEUCHOS_TEST_FOR_EXCEPTION
        (host_numExportPacketsPerLID_.dimension_0 () !=
         numExportPacketsPerLID_.dimension_0 (), std::logic_error, "Tpetra::"
         "DistObject::doTransferNew: host_numExportPacketsPerLID_.dimension_0()"
         " = " << host_numExportPacketsPerLID_.dimension_0 () << " != "
         "numExportPacketsPerLID_.dimension_0() = " <<
         numExportPacketsPerLID_.dimension_0() << ".  Please report this bug "
         "to the Tpetra developers.");

      // FIXME (mfh 28 Mar 2016) This helps fix #227.
      execution_space::fence ();

      // Copy numExportPacketsPerLID to host
      Kokkos::deep_copy (host_numExportPacketsPerLID_, numExportPacketsPerLID_);

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

      // FIXME (mfh 17 Feb 2014) Distributor doesn't actually inspect
      // the contents of the "exports" or "imports" arrays, other than
      // to do a deep copy in the (should be technically unnecessary,
      // but isn't for some odd reason) "self-message" case.
      // Distributor thus doesn't need host views; it could do just
      // fine with device views, assuming that MPI knows how to read
      // device memory (which doesn't even require UVM).

      if (needCommunication) {
        if (revOp == DoReverse) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            distor.doReversePostsAndWaits (create_const_view (host_numExportPacketsPerLID_),
                                           1,
                                           host_numImportPacketsPerLID_);
            size_t totalImportPackets = 0;
            // FIXME (mfh 17 Feb 2014) This would be a good place for
            // a Kokkos reduction.  numImportPacketsPerLID_ has as
            // many entries as the number of LIDs on the calling
            // process.
            for (view_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
              totalImportPackets += host_numImportPacketsPerLID_[i];
            }
            if (static_cast<size_t> (imports_.dimension_0 ()) != totalImportPackets) {
              Kokkos::realloc (imports_, totalImportPackets);
              // This is doTransferNew, so we don't need to allocate the
              // host version of imports_.
            }
            distor.doReversePostsAndWaits (create_const_view (exports_),
                                           getArrayView (host_numExportPacketsPerLID_),
                                           imports_,
                                           getArrayView (host_numImportPacketsPerLID_));
          }
          else {
            distor.doReversePostsAndWaits (create_const_view (exports_),
                                           constantNumPackets,
                                           imports_);
          }
        }
        else { // revOp == DoForward
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            distor.doPostsAndWaits (create_const_view (host_numExportPacketsPerLID_), 1,
                                    host_numImportPacketsPerLID_);
            size_t totalImportPackets = 0;
            // FIXME (mfh 17 Feb 2014) This would be a good place for
            // a Kokkos reduction.  numImportPacketsPerLID_ has as
            // many entries as the number of LIDs on the calling
            // process.
            for (view_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
              totalImportPackets += host_numImportPacketsPerLID_[i];
            }

            if (static_cast<size_t> (imports_.dimension_0 ()) != totalImportPackets) {
              Kokkos::realloc (imports_, totalImportPackets);
              // This is doTransferNew, so we don't need to allocate
              // the host version of imports_.
            }
            distor.doPostsAndWaits (create_const_view (exports_),
                                    getArrayView (host_numExportPacketsPerLID_),
                                    imports_,
                                    getArrayView (host_numImportPacketsPerLID_));
          }
          else {
            distor.doPostsAndWaits (create_const_view (exports_),
                                    constantNumPackets,
                                    imports_);
          }
        }

        // FIXME (mfh 17 Feb 2014) This array should just live on the
        // device and stay there.  There is no need for a host view,
        // as long as unpackAndCombine(new) knows what to do with a
        // device view.
        //
        // Copy numImportPacketsPerLID to device
        Kokkos::deep_copy (numImportPacketsPerLID_, host_numImportPacketsPerLID_);

        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          unpackAndCombineNew (remoteLIDs, imports_, numImportPacketsPerLID_,
                            constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    // FIXME (mfh 17 Dec(??? probably Feb) 2014) We don't have to call
    // releaseViews() on the destination object any more, since MPI
    // knows how to read device memory.

    this->releaseViews ();
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
