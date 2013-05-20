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

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Distributor.hpp"

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_DistObject_decl.hpp"
#endif // DOXYGEN_USE_ONLY


namespace Tpetra {
  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  DistObject (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map)
    : map_ (map)
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  DistObject (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source)
    : map_ (source.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::~DistObject()
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "Tpetra::DistObject<"
       << TypeNameTraits<Packet>::name ()
       << ", " << TypeNameTraits<LocalOrdinal>::name ()
       << ", " << TypeNameTraits<GlobalOrdinal>::name ()
       << ", " << TypeNameTraits<Node>::name ()
       << ">";
    return os.str ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::rcpFromRef;
    using std::endl;

    const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ?
      Teuchos::VERB_LOW : verbLevel;

    if (vl != Teuchos::VERB_NONE) {
      out << this->description () << endl;
      Teuchos::OSTab tab (rcpFromRef (out));
      out << "Export buffer size (in packets): " << exports_.size() << endl
          << "Import buffer size (in packets): " << imports_.size() << endl
          << "Map over which this object is distributed:" << endl;
      map_->describe (out, vl);
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Tpetra::DistObject::removeEmptyProcessesInPlace: Not implemented");
  }

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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(*getMap() != *importer.getTargetMap(),
      std::invalid_argument, "doImport: The target DistObject's Map is not "
      "identical to the Import's target Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(*A.getMap() != *importer.getSourceMap(),
      std::invalid_argument, "doImport: The source DistObject's Map is not "
      "identical to the Import's source Map.");
    size_t numSameIDs = importer.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    const view_type exportLIDs      = importer.getExportLIDs();
    const view_type remoteLIDs      = importer.getRemoteLIDs();
    const view_type permuteToLIDs   = importer.getPermuteToLIDs();
    const view_type permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs,
                      remoteLIDs, exportLIDs, importer.getDistributor (),
                      DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(   *getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doExport: The target DistObject's Map is not identical to the Export's target Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doExport: The source DistObject's Map is not identical to the Export's source Map.");
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getExportLIDs();
    view_type remoteLIDs      = exporter.getRemoteLIDs();
    view_type permuteToLIDs   = exporter.getPermuteToLIDs();
    view_type permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(  * getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doImport (with Export): The target DistObject's Map is not identical to the Export's source Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doImport (with Export): The source DistObject's Map is not identical to the Export's target Map.");
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getRemoteLIDs();
    view_type remoteLIDs      = exporter.getExportLIDs();
    view_type permuteToLIDs   = exporter.getPermuteFromLIDs();
    view_type permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( *getMap() != *importer.getSourceMap(),
      std::invalid_argument, "doExport (with Import): The target object's Map "
      "is not identical to the Import's source Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *importer.getTargetMap(),
      std::invalid_argument, "doExport (with Import): The source object's Map "
      "is not identical to the Import's target Map.");
    size_t numSameIDs = importer.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = importer.getRemoteLIDs();
    view_type remoteLIDs      = importer.getExportLIDs();
    view_type permuteToLIDs   = importer.getPermuteFromLIDs();
    view_type permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, importer.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::isDistributed() const {
    return map_->isDistributed ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  constantNumberOfPackets () const {
    return 0; // default implementation; subclasses may override
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doTransfer (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source,
              CombineMode CM,
              size_t numSameIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
              Distributor &distor,
              ReverseOption revOp)
  {
    using Teuchos::as;
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    TEUCHOS_TEST_FOR_EXCEPTION( ! checkSizes(source), std::invalid_argument,
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");
    Kokkos::ReadWriteOption rwo = Kokkos::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      const size_t numIDsToWrite = numSameIDs +
        as<size_t> (permuteToLIDs.size ()) +
        as<size_t> (remoteLIDs.size ());
      if (numIDsToWrite == this->getMap ()->getNodeNumElements ()) {
        // We're overwriting all of our local data in the destination
        // object, so a write-only view suffices.
        //
        // FIXME (mfh 10 Apr 2012) This doesn't make sense for a
        // CrsMatrix with a dynamic graph.  INSERT mode could mean
        // that we're adding new entries to the object, but we don't
        // want to get rid of the old ones.
        rwo = Kokkos::WriteOnly;
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
    source.createViews();

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
    this->createViewsNonConst(rwo);

    if (numSameIDs + permuteToLIDs.size()) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
      // There is at least one GID to copy or permute.
      copyAndPermute (source, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }
    // The method may return zero even if the implementation actually
    // does have a constant number of packets per LID.  However, if it
    // returns nonzero, we may use this information to avoid
    // (re)allocating num{Ex,Im}portPacketsPerLID_.  packAndPrepare()
    // will set this to its final value.
    size_t constantNumPackets = this->constantNumberOfPackets ();

    if (constantNumPackets == 0) {
      numExportPacketsPerLID_.resize (exportLIDs.size ());
      numImportPacketsPerLID_.resize (remoteLIDs.size ());
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
      packAndPrepare (source, exportLIDs, exports_, numExportPacketsPerLID_(),
                      constantNumPackets, distor);
    }

    // We don't need the source's data anymore, so it can let go of
    // its views.  On a discrete accelerator, this frees host memory,
    // since device memory has the "master" version of the data.
    source.releaseViews();

    if (constantNumPackets != 0) {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = remoteLIDs.size() * constantNumPackets;
      if (as<size_t> (imports_.size()) != rbufLen) {
        imports_.resize (rbufLen);
      }
    }
    if ((isDistributed() && revOp == DoReverse) ||
        (source.isDistributed() && revOp == DoForward)) {
      // call one of the doPostsAndWaits functions
      if (revOp == DoReverse) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doReversePostsAndWaits (numExportPacketsPerLID_().getConst(), 1,
                                         numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doReversePostsAndWaits (exports_().getConst(),
                                         numExportPacketsPerLID_(),
                                         imports_(),
                                         numImportPacketsPerLID_());
        }
        else {
          distor.doReversePostsAndWaits (exports_().getConst(),
                                         constantNumPackets,
                                         imports_());
        }
      }
      else { // revOp == DoForward
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doPostsAndWaits (numExportPacketsPerLID_().getConst(), 1,
                                  numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doPostsAndWaits (exports_().getConst(),
                                  numExportPacketsPerLID_(),
                                  imports_(),
                                  numImportPacketsPerLID_());
        }
        else {
          distor.doPostsAndWaits (exports_().getConst(),
                                  constantNumPackets,
                                  imports_());
        }
      }
      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        unpackAndCombine (remoteLIDs, imports_(), numImportPacketsPerLID_(),
                          constantNumPackets, distor, CM);
      }
    }
    this->releaseViews();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::print (std::ostream &os) const
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using std::endl;

    RCP<FancyOStream> out = getFancyOStream (rcpFromRef (os));
    this->describe (*out, Teuchos::VERB_DEFAULT);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::createViews () const
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  createViewsNonConst (Kokkos::ReadWriteOption /*rwo*/)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  releaseViews () const
  {}

#define TPETRA_DISTOBJECT_INSTANT(SCALAR, LO, GO, NODE) \
  \
  template class DistObject< SCALAR , LO , GO , NODE >;

 // The "SLGN" stuff above doesn't work for Packet=char.
#define TPETRA_DISTOBJECT_INSTANT_CHAR(LO, GO, NODE) \
  \
  template class DistObject< char , LO , GO , NODE >;


} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_DEF_HPP */
