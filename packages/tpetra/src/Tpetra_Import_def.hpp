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

#ifndef TPETRA_IMPORT_DEF_HPP
#define TPETRA_IMPORT_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
#  include <Tpetra_Import_decl.hpp>
#endif // DOXYGEN_USE_ONLY

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Export.hpp>
#include <Teuchos_as.hpp>

namespace {
  // Default value of Import's "Debug" parameter.
  const bool tpetraImportDebugDefault = false;
} // namespace (anonymous)

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    bool debug = tpetraImportDebugDefault;
    if (! plist.is_null ()) {
      try {
        debug = plist->get<bool> ("Debug");
      } catch (Teuchos::Exceptions::InvalidParameter&) {}
    }
    debug_ = debug;
    ImportData_->distributor_.setParameterList (plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  init (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
        const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
        const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using Teuchos::rcp;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    // Read "Debug" parameter from the input ParameterList.
    bool debug = tpetraImportDebugDefault;
    if (! plist.is_null ()) {
      try {
        debug = plist->get<bool> ("Debug");
      } catch (Teuchos::Exceptions::InvalidParameter&) {}
    }
    debug_ = debug;

    if (! out_.is_null ()) {
      out_->pushTab ();
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Import ctor" << endl;
      *out_ << os.str ();
    }
    ImportData_ = rcp (new data_type (source, target, out_, plist));
    setupSamePermuteRemote ();
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Import ctor: "
         << "setupSamePermuteRemote done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupExport ();
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Import ctor: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
    remoteGIDs_ = null; // Don't need this anymore.
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraImportDebugDefault)
  {
    init (source, target, Teuchos::null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const RCP<Teuchos::FancyOStream>& out) :
    out_ (out),
    debug_ (tpetraImportDebugDefault)
  {
    init (source, target, Teuchos::null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraImportDebugDefault)
  {
    init (source, target, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (out),
    debug_ (tpetraImportDebugDefault)
  {
    init (source, target, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Import<LocalOrdinal,GlobalOrdinal,Node>& rhs)
    : ImportData_ (rhs.ImportData_)
    , out_ (rhs.out_)
    , debug_ (rhs.debug_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter)
    : out_ (exporter.out_)
    , debug_ (exporter.debug_)
  {
    if (! exporter.ExportData_.is_null ()) {
      ImportData_ = exporter.ExportData_->reverseClone ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::~Import()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumSameIDs() const {
    return ImportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumPermuteIDs() const {
    return ImportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ImportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ImportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ImportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const int>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportPIDs() const {
    return ImportData_->exportPIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>&
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  operator= (const Import<LocalOrdinal,GlobalOrdinal,Node>& rhs) {
    if (&rhs != this) {
      ImportData_ = rhs.ImportData_;
    }
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::
  print (std::ostream& os) const
  {
    using Teuchos::Comm;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using Teuchos::toString;
    using std::endl;

    RCP<const Comm<int> > comm = getSourceMap ()->getComm ();
    const int myImageID = comm->getRank ();
    const int numImages = comm->getSize ();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        os << endl;
        if (myImageID == 0) { // I'm the root node (only output this info once)
          os << "Import Data Members:" << endl;
        }
        os << "Image ID       : " << myImageID << endl;

        os << "permuteFromLIDs: " << toString (getPermuteFromLIDs ()) << endl;
        os << "permuteToLIDs  : " << toString (getPermuteToLIDs ()) << endl;
        os << "remoteLIDs     : " << toString (getRemoteLIDs ()) << endl;
        os << "exportLIDs     : " << toString (getExportLIDs ()) << endl;
        os << "exportPIDs     : " << toString (getExportPIDs ()) << endl;

        os << "numSameIDs     : " << getNumSameIDs () << endl;
        os << "numPermuteIDs  : " << getNumPermuteIDs () << endl;
        os << "numRemoteIDs   : " << getNumRemoteIDs () << endl;
        os << "numExportIDs   : " << getNumExportIDs () << endl;
      }
      // A few global barriers give output a chance to complete.
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }

    const bool printMaps = false;
    if (printMaps) {
      if (myImageID == 0) {
        os << endl << endl << "Source Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getSourceMap();
      comm->barrier();

      if (myImageID == 0) {
        os << endl << endl << "Target Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getTargetMap();
      comm->barrier();
    }

    // It's also helpful for debugging to print the Distributor
    // object.  Epetra_Import::Print() does this, so we can do a
    // side-by-side comparison.
    if (myImageID == 0) {
      os << endl << endl << "Distributor:" << endl << std::flush;
    }
    comm->barrier();
    getDistributor().describe (*(getFancyOStream (rcpFromRef (os))),
                               Teuchos::VERB_EXTREME);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setupSamePermuteRemote()
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::null;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;
    const Map<LO,GO,Node>& source = * (getSourceMap ());
    const Map<LO,GO,Node>& target = * (getTargetMap ());
    ArrayView<const GO> sourceGIDs = source.getNodeElementList ();
    ArrayView<const GO> targetGIDs = target.getNodeElementList ();

#ifdef HAVE_TPETRA_DEBUG
    ArrayView<const GO> rawSrcGids = sourceGIDs;
    ArrayView<const GO> rawTgtGids = targetGIDs;
#else
    const GO* const rawSrcGids = sourceGIDs.getRawPtr ();
    const GO* const rawTgtGids = targetGIDs.getRawPtr ();
#endif // HAVE_TPETRA_DEBUG
    const size_type numSrcGids = sourceGIDs.size ();
    const size_type numTgtGids = targetGIDs.size ();
    const size_type numGids = std::min (numSrcGids, numTgtGids);

    // Compute numSameIDs_: the number of initial GIDs that are the
    // same (and occur in the same order) in both Maps.  The point of
    // numSameIDs_ is for the common case of an Import where all the
    // overlapping GIDs are at the end of the target Map, but
    // otherwise the source and target Maps are the same.  This allows
    // a fast contiguous copy for the initial "same IDs."
    size_type numSameGids = 0;
    for ( ; numSameGids < numGids && rawSrcGids[numSameGids] == rawTgtGids[numSameGids]; ++numSameGids)
      {} // third clause of 'for' does everything
    ImportData_->numSameIDs_ = numSameGids;

    // Compute permuteToLIDs_, permuteFromLIDs_, remoteGIDs_, and
    // remoteLIDs_.  The first two arrays are IDs to be permuted, and
    // the latter two arrays are IDs to be received ("imported"),
    // called "remote" IDs.
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  IDs to receive are in the
    // target Map, but not the source Map.

    remoteGIDs_ = rcp (new Array<GO> ());
    Array<GO>& remoteGids = *remoteGIDs_;
    Array<LO>& permuteToLIDs = ImportData_->permuteToLIDs_;
    Array<LO>& permuteFromLIDs = ImportData_->permuteFromLIDs_;
    Array<LO>& remoteLIDs = ImportData_->remoteLIDs_;
    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO numTgtLids = as<LO> (numTgtGids);
    // Iterate over the target Map's LIDs, since we only need to do
    // GID -> LID lookups for the source Map.
    for (LO tgtLid = numSameGids; tgtLid < numTgtLids; ++tgtLid) {
      const GO curTargetGid = rawTgtGids[tgtLid];
      // getLocalElement() returns LINVALID if the GID isn't in the source Map.
      // This saves us a lookup (which isNodeGlobalElement() would do).
      const LO srcLid = source.getLocalElement (curTargetGid);
      if (srcLid != LINVALID) { // if source.isNodeGlobalElement (curTargetGid)
        permuteToLIDs.push_back (tgtLid);
        permuteFromLIDs.push_back (srcLid);
      } else {
        remoteGids.push_back (curTargetGid);
        remoteLIDs.push_back (tgtLid);
      }
    }

    TPETRA_ABUSE_WARNING(
      getNumRemoteIDs() > 0 && ! source.isDistributed(),
      std::runtime_error,
      "::setupSamePermuteRemote(): Target has remote LIDs but Source is not "
      "distributed globally." << std::endl
      << "Importing to a submap of the target map.");
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupExport() {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Array<int>::difference_type size_type;
    const Map<LO, GO, Node> & source = *getSourceMap ();

    Teuchos::OSTab tab (out_);

    if (debug_) {
      std::ostringstream os;
      const int myRank = source.getComm ()->getRank ();
      os << myRank << ": Import::setupExport" << endl;
      *out_ << os.str ();
    }

    // For each entry remoteGIDs[i], remoteProcIDs[i] will contain
    // the process ID of the process that owns that GID.
    ArrayView<GO> remoteGIDs = (*remoteGIDs_) ();
    Array<int> remoteProcIDs (remoteGIDs.size ());

    // lookup == IDNotPresent means that the source Map wasn't able to
    // figure out to which processes one or more of the GIDs in the
    // given list of remoteGIDs belong.
    //
    // The previous abuse warning said "The target Map has GIDs not
    // found in the source Map."  This statement could be confusing,
    // because it doesn't refer to ownership by the current process,
    // but rather to ownership by _any_ process participating in the
    // Map.  (It could not possibly refer to ownership by the current
    // process, since remoteGIDs is exactly the list of GIDs owned by
    // the target Map but not owned by the source Map.  It was
    // constructed that way by setupSamePermuteRemote().)
    //
    // What this statement means is that the source and target Maps
    // don't contain the same set of GIDs globally (over all
    // processes).  That is, there is at least one GID owned by some
    // process in the target Map, which is not owned by _any_ process
    // in the source Map.
    const LookupStatus lookup =
      source.getRemoteIndexList (remoteGIDs, remoteProcIDs ());
    if (debug_) {
      std::ostringstream os;
      const int myRank = source.getComm ()->getRank ();
      os << myRank << ": Import::setupExport: finished lookup" << endl;
      *out_ << os.str ();
    }

    TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
      "::setupExport(): the source Map wasn't able to figure out which process "
      "owns one or more of the GIDs in the list of remote GIDs.  This probably "
      "means that there is at least one GID owned by some process in the target"
      " Map which is not owned by any process in the source Map.  (That is, the"
      " source and target Maps do not contain the same set of GIDs globally.)");

    // Ignore remote GIDs that aren't owned by any process in the
    // source Map.  getRemoteIndexList() gives each of these a process
    // ID of -1.
    if (lookup == IDNotPresent) {
      const size_type numInvalidRemote =
        std::count_if (remoteProcIDs.begin (), remoteProcIDs.end (),
                       std::bind1st (std::equal_to<int> (), -1));
      // If all of them are invalid, we can delete the whole array.
      const size_type totalNumRemote = getNumRemoteIDs ();
      if (numInvalidRemote == totalNumRemote) {
        // all remotes are invalid; we have no remotes; we can delete the remotes
        remoteProcIDs.clear ();
        (*remoteGIDs_).clear (); // This invalidates the view remoteGIDs
        ImportData_->remoteLIDs_.clear();
      }
      else {
        // Some remotes are valid; we need to keep the valid ones.
        // Pack and resize remoteProcIDs, remoteGIDs_, and
        // remoteLIDs_.
        size_type numValidRemote = 0;
#ifdef HAVE_TPETRA_DEBUG
        ArrayView<GlobalOrdinal> remoteGIDsPtr = remoteGIDs;
#else
        GlobalOrdinal* const remoteGIDsPtr = remoteGIDs.getRawPtr ();
#endif // HAVE_TPETRA_DEBUG
        for (size_type r = 0; r < totalNumRemote; ++r) {
          // Pack in all the valid remote PIDs and GIDs.
          if (remoteProcIDs[r] != -1) {
            remoteProcIDs[numValidRemote] = remoteProcIDs[r];
            remoteGIDsPtr[numValidRemote] = remoteGIDsPtr[r];
            ImportData_->remoteLIDs_[numValidRemote] = ImportData_->remoteLIDs_[r];
            ++numValidRemote;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          numValidRemote != totalNumRemote - numInvalidRemote, std::logic_error,
          "Tpetra::Import::setupExport(): After removing invalid remote GIDs and"
          " packing the valid remote GIDs, numValidRemote = " << numValidRemote
          << " != totalNumRemote - numInvalidRemote = "
          << totalNumRemote - numInvalidRemote
          << ".  Please report this bug to the Tpetra developers.");

        remoteProcIDs.resize (numValidRemote);
        (*remoteGIDs_).resize (numValidRemote);
        ImportData_->remoteLIDs_.resize (numValidRemote);
      }
      // Revalidate the view after clear or resize.
      remoteGIDs = (*remoteGIDs_)();
    }

    // Sort remoteProcIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs_ and remoteLIDs_.  This ensures that
    // remoteProcIDs[i], remoteGIDs_[i], and remoteLIDs_[i] all refer
    // to the same thing.
    sort3 (remoteProcIDs.begin(),
           remoteProcIDs.end(),
           remoteGIDs.begin(),
           ImportData_->remoteLIDs_.begin());

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDs and remoteProcIDs_ are input;
    // exportGIDs and exportProcIDs_ are output arrays which are
    // allocated by createFromRecvs().
    Array<GO> exportGIDs;
    ImportData_->distributor_.createFromRecvs (remoteGIDs ().getConst (),
                                               remoteProcIDs, exportGIDs,
                                               ImportData_->exportPIDs_);

    // Find the LIDs corresponding to the (outgoing) GIDs in
    // exportGIDs.  For sparse matrix-vector multiply, this tells the
    // calling process how to index into the source vector to get the
    // elements which it needs to send.
    const size_type numExportIDs = exportGIDs.size ();
    if (numExportIDs > 0) {
      ImportData_->exportLIDs_.resize(numExportIDs);

      ArrayView<const GO> expGIDs = exportGIDs ();
      ArrayView<LO> expLIDs = ImportData_->exportLIDs_ ();
      for (size_type k = 0; k < numExportIDs; ++k) {
        expLIDs[k] = source.getLocalElement (expGIDs[k]);
      }
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source.getComm ()->getRank ();
      os << myRank << ": Import::setupExport: done" << endl;
      *out_ << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setUnionNaiveImpl (const Import<LocalOrdinal, GlobalOrdinal, Node>& rhs) const
  {
    // mfh 22 Apr 2013: We plan to optimize this in the future, but
    // provide an unoptimized implementation for now to allow
    // development to proceed.
    typedef Tpetra::global_size_t GST;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Import<LO, GO, Node> import_type;

    RCP<const map_type> srcMap = this->getSourceMap ();
    RCP<const map_type> tgtMap1 = this->getTargetMap ();
    RCP<const map_type> tgtMap2 = rhs.getTargetMap ();
    RCP<const Comm<int> > comm = srcMap->getComm ();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! srcMap->isSameAs (* (rhs.getSourceMap ())), std::invalid_argument,
      "Tpetra::Import::setUnion: The source Map of the input Import must be the "
      "same as (in the sense of Map::isSameAs) the source Map of this Import.");
#endif // HAVE_TPETRA_DEBUG

    ArrayView<const GO> gids1 = tgtMap1->getNodeElementList ();
    ArrayView<const GO> gids2 = tgtMap2->getNodeElementList ();

    // Union needs sorted sequences.
    // Only copy a view if it's not sorted.
    Array<GO> gids1Copy, gids2Copy;
    if (! SortDetails::isAlreadySorted (gids1.begin (), gids1.end ())) {
      gids1Copy.reserve (gids1.size ());
      gids1Copy.assign (gids1.begin (), gids1.end ());
      std::sort (gids1Copy.begin (), gids1Copy.end ());
      gids1 = gids1Copy ();
    }
    if (! SortDetails::isAlreadySorted (gids2.begin (), gids2.end ())) {
      gids2Copy.reserve (gids2.size ());
      gids2Copy.assign (gids2.begin (), gids2.end ());
      std::sort (gids2Copy.begin (), gids2Copy.end ());
      gids2 = gids2Copy ();
    }
    // Compute union.  Reservation is only a reasonable guess.
    Array<GO> gidsUnion;
    gidsUnion.reserve (std::max (gids1.size (), gids2.size ()));
    std::set_union (gids1.begin (), gids1.end (), gids2.begin (), gids2.end (),
                    std::back_inserter (gidsUnion));
    // Find the index base, which must also be the Map's global min
    // GID.  This is inefficient, because we could have done that
    // while computing the set union above.
    typename Array<GO>::const_iterator minPos =
      std::min_element (gidsUnion.begin (), gidsUnion.end ());
    const GO myMinGid = (minPos == gidsUnion.end ()) ? 0 : *minPos;
    GO globalMinGid = 0;
    reduceAll<int, GO> (*comm, REDUCE_MIN, myMinGid, outArg (globalMinGid));
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const map_type> unionMap =
      rcp (new map_type (INVALID, gidsUnion (), globalMinGid,
                         comm, srcMap->getNode ()));
    return rcp (new import_type (srcMap, unionMap));
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setUnion (const Import<LocalOrdinal, GlobalOrdinal, Node>& rhs) const
  {
    // mfh 22 Apr 2013: We plan to optimize this in the future, but
    // provide an unoptimized implementation for now to allow
    // development to proceed.
    return setUnionNaiveImpl (rhs);

//     using Teuchos::Array;
//     using Teuchos::ArrayView;
//     using Teuchos::RCP;
//     using Teuchos::rcp;
//     typedef Import<class LocalOrdinal, class GlobalOrdinal, class Node> import_type;
//     typedef Map<class LocalOrdinal, class GlobalOrdinal, class Node> map_type;

//     const map_type& srcMap = * (this->getSourceMap ());
//     const map_type& tgtMap1 = * (this->getTargetMap ());
//     const map_type& tgtMap2 = * (rhs.getTargetMap ());
// #ifdef HAVE_TPETRA_DEBUG
//     TEUCHOS_TEST_FOR_EXCEPTION(
//       srcMap->isSameAs (* (rhs.getSourceMap ())), std::invalid_argument,
//       "Tpetra::Import::union: The source Map of the input Import must be the "
//       "same as (in the sense of Map::isSameAs) the source Map of this Import.");
// #endif // HAVE_TPETRA_DEBUG

//     ArrayView<const GO> gids1 = tgtMap1.getNodeElementList ();
//     ArrayView<const GO> gids2 = tgtMap2.getNodeElementList ();

//     ArrayView<const LO> permuteToLIDs1 = this->getPermuteToLIDs ();
//     ArrayView<const LO> permuteToLIDs2 = rhs.getPermuteToLIDs ();
//     const size_type numPermuteToLIDs1 = permuteToLIDs1.size ();
//     const size_type numPermuteToLIDs2 = permuteToLIDs2.size ();
//     ArrayView<const LO> remoteLIDs1 = this->getRemoteLIDs ();
//     ArrayView<const LO> remoteLIDs2 = rhs.getRemoteLIDs ();
//     const size_type numRemoteLIDs1 = remoteLIDs1.size ();
//     const size_type numRemoteLIDs2 = remoteLIDs2.size ();

//     // We could get the remote PIDs from the two Imports' Distributor
//     // instances.  This would help us bypass much of the work of
//     // createFromRecvs().

//     // The union has to be at least as big as the max of the two
//     // sizes.  This is just an optimization so we don't have to get
//     // the size exactly right here.
//     Array<GO> gidsUnion;
//     gidsUnion.reserve (std::max (gids1.size (), gids2.size ()));

//     // Get the "same IDs" for each Map: that is, those GIDs that are
//     // initially the same and in the same order as in the source Map.
//     // We'll add the larger of these two sets to the union Map.
//     const size_t numSame1 = tgtMap1.getNumSameIDs ();
//     const size_t numSame2 = tgtMap2.getNumSameIDs ();
//     const size_t numSameUnion = std::max (numSame1, numSame2);
//     ArrayView<const GO> gids1Same = gids1.view (0, numSame1);
//     ArrayView<const GO> gids2Same = gids2.view (0, numSame2);

//     // Add the larger of the two "same" lists to the union GIDs.
//     // Don't sort this part of the union list; we want to keep it in
//     // the same order so that max(numSame1,numSame2) is the new number
//     // of same IDs as in the source Map.
//     if (numSame1 >= numSame2) {
//       gidsUnion.assign (gids1Same.begin (), gids1Same.end ());
//     } else {
//       gidsUnion.assign (gids2Same.begin (), gids2Same.end ());
//     }

//     // GIDs that remain from each Map after the "same" lists.  These
//     // are either "permute-to" GIDs (in both source and target Maps,
//     // but in a different order) or "remote" GIDs (in target Map, not
//     // in source Map, so they must be received).
//     ArrayView<const GO> gids1Diff = gids1.view (numSame1, gids1.size () - numSame1);
//     ArrayView<const GO> gids2Diff = gids2.view (numSame2, gids2.size () - numSame2);

//     // Now we have to figure out the permute-to and remote IDs in the
//     // union Map.  We could just imitate setupSamePermuteRemote(),
//     // except going through both target Maps instead of just one.
//     // However, we already have permuteToLIDs and remoteLIDs for each
//     // of the target Maps, because we have both of their Import
//     // objects.  A straightforward setupSamePermuteRemote() analogue
//     // would require a lot of GID -> LID lookups.  LID -> GID lookups
//     // are much cheaper.

//     // We must convert the permute-to LIDs to GIDs, since their LIDs
//     // might be different in the union Map.
//     Array<GO> permuteToGIDs1 (numPermuteToLIDs1);
//     for (size_type k = 0; k < numPermuteToLIDs1; ++k) {
//       permuteToGIDs1[k] = getGlobalElement (permuteToLIDs1[k]);
//     }
//     Array<GO> permuteToGIDs2 (numPermuteToLIDs2);
//     for (size_type k = 0; k < numPermuteToLIDs2; ++k) {
//       permuteToGIDs2[k] = getGlobalElement (permuteToLIDs2[k]);
//     }

//     // Sort the two permute-to GID lists and compute their union.
//     std::sort (permuteToGIDs1.begin (), permuteToGIDs1.end ());
//     std::sort (permuteToGIDs2.begin (), permuteToGIDs2.end ());
//     std::set_union (permuteToGIDs1.begin (), permuteToGIDs1.end (),
//                     permuteToGIDs2.begin (), permuteToGIDs2.end (),
//                     std::back_inserter (gidsUnion));
//     ArrayView<const GO> gidsPermuteToUnion =
//       gidsUnion.view (numSameUnion, gidsUnion.size () - numSameUnion);
//     const size_type numPermuteToUnion = gidsPermuteToUnion.size ();

//     // We may actually compute permuteToLIDs in the union Map now,
//     // since we have all the num-same and permute-to GIDs in the union
//     // Map in their correct order.  While we're at it, we can compute
//     // permuteFromLIDs (LIDs in the source Map).
//     const size_type numPermuteToUnion = gidsPermuteToUnion.size ();
//     Array<LO> lidsPermuteToUnion (numPermuteToUnion);
//     Array<LO> lidsPermuteFromUnion (numPermuteToUnion);
//     for (size_type k = 0; k < numPermuteToUnion; ++k) {
//       lidsPermuteToUnion[k] = as<LO> (numSameUnion) + as<LO> (k);
//       lidsPermuteFromUnion[k] = srcMap.getLocalElement (gidsPermuteToUnion[k]);
//     }

//     // At this point, we have the following:
//     //   - All same-as and permute-to GIDs in the union target Map
//     //   - Permute-to LIDs (in the union target Map) for the union Import
//     //   - Permute-from LIDs (in the source Map) for the union Import

//     // Compute remote IDs (in source Map, not in the union target Map
//     // -- that is, not in either of the two target Maps).  We already
//     // have the remote LIDs from the two Imports, and LID -> GID
//     // conversion is cheap.

//     Array<GO> remoteGIDs1 (numRemoteLIDs1);
//     for (size_type k = 0; k < numRemoteLIDs1; ++k) {
//       remoteGIDs1[k] = getGlobalElement (remoteLIDs1[k]);
//     }
//     Array<GO> remoteGIDs2 (numRemoteLIDs2);
//     for (size_type k = 0; k < numRemoteLIDs2; ++k) {
//       remoteGIDs2[k] = getGlobalElement (remoteLIDs2[k]);
//     }

//     // Sort the two remote GID lists and compute their union.
//     std::sort (remoteGIDs1.begin (), remoteGIDs1.end ());
//     std::sort (remoteGIDs2.begin (), remoteGIDs2.end ());
//     std::set_union (remoteGIDs1.begin (), remoteGIDs1.end (),
//                     remoteGIDs2.begin (), remoteGIDs2.end (),
//                     std::back_inserter (gidsUnion));
//     const size_type numNonRemoteUnion = numSameUnion + numPermuteUnion;
//     const size_type numRemoteUnion =
//       gidsUnion.size () - numNonRemoteUnion;
//     ArrayView<const GO> gidsRemoteUnion =
//       gidsUnion.view (numNonRemoteUnion, numRemoteUnion);

//     // Compute remote LIDs (in the union target Map).
//     Array<LO> lidsRemoteUnion (numRemoteUnion);
//     for (size_type k = 0; k < numRemoteUnion; ++k) {
//       lidsRemoteUnion[k] = as<LO> (numNonRemoteUnion) + as<LO> (k);
//     }

//     // At this point, we have the following:
//     //   - All same-as and permute-to GIDs in the union target Map
//     //   - Permute-to LIDs (in the union target Map) for the union Import
//     //   - Permute-from LIDs (in the source Map) for the union Import
//     //   - Remote GIDs and their LIDs in the union target Map

//     // Now we need:
//     //   - Remote PIDs (in the source Map)
//     //   - To pack (eliminate invalid) and resize remote PID, GIDs,
//     //     and LIDs
//     //   - To sort remote PIDs, and apply the resulting permutation to
//     //     remote GIDs and LIDs
//     //   - To create the Distributor




//     // Same IDs and permute IDs are disjoint in each Map.
//     // permuteToLIDs (from each the target Map) occur in order, but
//     // not necessarily the same order as in the whole element list,
//     // because they might be interspersed with remote IDs.


//     Array<GO> gids1DiffCopy (gids1Diff.begin (), gids1Diff.end ());
//     Array<GO> gids2DiffCopy (gids2Diff.begin (), gids2Diff.end ());

//     std::sort (gids1DiffCopy.begin (), gids1DiffCopy.end ());
//     std::sort (gids2DiffCopy.begin (), gids2DiffCopy.end ());
//     Array<GO> gidsDiffUnion;
//     gidsUnion.reserve (gids1DiffCopy.size ());
//     std::set_union (gids1DiffCopy.begin (), gids1DiffCopy.end (),
//                     gids2DiffCopy.begin (), gids2DiffCopy.end (),
//                     std::back_inserter (gidsDiffUnion));




//     RCP<const map_type> unionMap = this->getTargetMap ()->union (rhs.getTargetMap ());
//     return rcp (new import_type (srcMap, unionMap));
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setUnionImpl (const Import<LocalOrdinal, GlobalOrdinal, Node>& rhs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Optimized "
      "implementation of setUnion() has not yet been implemented.  "
      "setUnion() should call setUnionNaiveImpl() for now.  "
      "Please report this bugt to the Tpetra developers.");
    return Teuchos::null; // forestall "no return value" warning

//     using Teuchos::Array;
//     using Teuchos::ArrayView;
//     using Teuchos::RCP;
//     using Teuchos::rcp;
//     typedef Import<class LocalOrdinal, class GlobalOrdinal, class Node> import_type;
//     typedef Map<class LocalOrdinal, class GlobalOrdinal, class Node> map_type;

//     const map_type& srcMap = * (this->getSourceMap ());
//     const map_type& tgtMap1 = * (this->getTargetMap ());
//     const map_type& tgtMap2 = * (rhs.getTargetMap ());
// #ifdef HAVE_TPETRA_DEBUG
//     TEUCHOS_TEST_FOR_EXCEPTION(
//       srcMap->isSameAs (* (rhs.getSourceMap ())), std::invalid_argument,
//       "Tpetra::Import::union: The source Map of the input Import must be the "
//       "same as (in the sense of Map::isSameAs) the source Map of this Import.");
// #endif // HAVE_TPETRA_DEBUG

//     ArrayView<const GO> gids1 = tgtMap1.getNodeElementList ();
//     ArrayView<const GO> gids2 = tgtMap2.getNodeElementList ();

//     ArrayView<const LO> permuteToLIDs1 = this->getPermuteToLIDs ();
//     ArrayView<const LO> permuteToLIDs2 = rhs.getPermuteToLIDs ();
//     const size_type numPermuteToLIDs1 = permuteToLIDs1.size ();
//     const size_type numPermuteToLIDs2 = permuteToLIDs2.size ();
//     ArrayView<const LO> remoteLIDs1 = this->getRemoteLIDs ();
//     ArrayView<const LO> remoteLIDs2 = rhs.getRemoteLIDs ();
//     const size_type numRemoteLIDs1 = remoteLIDs1.size ();
//     const size_type numRemoteLIDs2 = remoteLIDs2.size ();

//     // We could get the remote PIDs from the two Imports' Distributor
//     // instances.  This would help us bypass much of the work of
//     // createFromRecvs().

//     // The union has to be at least as big as the max of the two
//     // sizes.  This is just an optimization so we don't have to get
//     // the size exactly right here.
//     Array<GO> gidsUnion;
//     gidsUnion.reserve (std::max (gids1.size (), gids2.size ()));

//     // Get the "same IDs" for each Map: that is, those GIDs that are
//     // initially the same and in the same order as in the source Map.
//     // We'll add the larger of these two sets to the union Map.
//     const size_t numSame1 = tgtMap1.getNumSameIDs ();
//     const size_t numSame2 = tgtMap2.getNumSameIDs ();
//     const size_t numSameUnion = std::max (numSame1, numSame2);
//     ArrayView<const GO> gids1Same = gids1.view (0, numSame1);
//     ArrayView<const GO> gids2Same = gids2.view (0, numSame2);

//     // Add the larger of the two "same" lists to the union GIDs.
//     // Don't sort this part of the union list; we want to keep it in
//     // the same order so that max(numSame1,numSame2) is the new number
//     // of same IDs as in the source Map.
//     if (numSame1 >= numSame2) {
//       gidsUnion.assign (gids1Same.begin (), gids1Same.end ());
//     } else {
//       gidsUnion.assign (gids2Same.begin (), gids2Same.end ());
//     }

//     // GIDs that remain from each Map after the "same" lists.  These
//     // are either "permute-to" GIDs (in both source and target Maps,
//     // but in a different order) or "remote" GIDs (in target Map, not
//     // in source Map, so they must be received).
//     ArrayView<const GO> gids1Diff = gids1.view (numSame1, gids1.size () - numSame1);
//     ArrayView<const GO> gids2Diff = gids2.view (numSame2, gids2.size () - numSame2);

//     // Now we have to figure out the permute-to and remote IDs in the
//     // union Map.  We could just imitate setupSamePermuteRemote(),
//     // except going through both target Maps instead of just one.
//     // However, we already have permuteToLIDs and remoteLIDs for each
//     // of the target Maps, because we have both of their Import
//     // objects.  A straightforward setupSamePermuteRemote() analogue
//     // would require a lot of GID -> LID lookups.  LID -> GID lookups
//     // are much cheaper.

//     // We must convert the permute-to LIDs to GIDs, since their LIDs
//     // might be different in the union Map.
//     Array<GO> permuteToGIDs1 (numPermuteToLIDs1);
//     for (size_type k = 0; k < numPermuteToLIDs1; ++k) {
//       permuteToGIDs1[k] = getGlobalElement (permuteToLIDs1[k]);
//     }
//     Array<GO> permuteToGIDs2 (numPermuteToLIDs2);
//     for (size_type k = 0; k < numPermuteToLIDs2; ++k) {
//       permuteToGIDs2[k] = getGlobalElement (permuteToLIDs2[k]);
//     }

//     // Sort the two permute-to GID lists and compute their union.
//     std::sort (permuteToGIDs1.begin (), permuteToGIDs1.end ());
//     std::sort (permuteToGIDs2.begin (), permuteToGIDs2.end ());
//     std::set_union (permuteToGIDs1.begin (), permuteToGIDs1.end (),
//                     permuteToGIDs2.begin (), permuteToGIDs2.end (),
//                     std::back_inserter (gidsUnion));
//     ArrayView<const GO> gidsPermuteToUnion =
//       gidsUnion.view (numSameUnion, gidsUnion.size () - numSameUnion);
//     const size_type numPermuteToUnion = gidsPermuteToUnion.size ();

//     // We may actually compute permuteToLIDs in the union Map now,
//     // since we have all the num-same and permute-to GIDs in the union
//     // Map in their correct order.  While we're at it, we can compute
//     // permuteFromLIDs (LIDs in the source Map).
//     const size_type numPermuteToUnion = gidsPermuteToUnion.size ();
//     Array<LO> lidsPermuteToUnion (numPermuteToUnion);
//     Array<LO> lidsPermuteFromUnion (numPermuteToUnion);
//     for (size_type k = 0; k < numPermuteToUnion; ++k) {
//       lidsPermuteToUnion[k] = as<LO> (numSameUnion) + as<LO> (k);
//       lidsPermuteFromUnion[k] = srcMap.getLocalElement (gidsPermuteToUnion[k]);
//     }

//     // At this point, we have the following:
//     //   - All same-as and permute-to GIDs in the union target Map
//     //   - Permute-to LIDs (in the union target Map) for the union Import
//     //   - Permute-from LIDs (in the source Map) for the union Import

//     // Compute remote IDs (in source Map, not in the union target Map
//     // -- that is, not in either of the two target Maps).  We already
//     // have the remote LIDs from the two Imports, and LID -> GID
//     // conversion is cheap.

//     Array<GO> remoteGIDs1 (numRemoteLIDs1);
//     for (size_type k = 0; k < numRemoteLIDs1; ++k) {
//       remoteGIDs1[k] = getGlobalElement (remoteLIDs1[k]);
//     }
//     Array<GO> remoteGIDs2 (numRemoteLIDs2);
//     for (size_type k = 0; k < numRemoteLIDs2; ++k) {
//       remoteGIDs2[k] = getGlobalElement (remoteLIDs2[k]);
//     }

//     // Sort the two remote GID lists and compute their union.
//     std::sort (remoteGIDs1.begin (), remoteGIDs1.end ());
//     std::sort (remoteGIDs2.begin (), remoteGIDs2.end ());
//     std::set_union (remoteGIDs1.begin (), remoteGIDs1.end (),
//                     remoteGIDs2.begin (), remoteGIDs2.end (),
//                     std::back_inserter (gidsUnion));
//     const size_type numNonRemoteUnion = numSameUnion + numPermuteUnion;
//     const size_type numRemoteUnion =
//       gidsUnion.size () - numNonRemoteUnion;
//     ArrayView<const GO> gidsRemoteUnion =
//       gidsUnion.view (numNonRemoteUnion, numRemoteUnion);

//     // Compute remote LIDs (in the union target Map).
//     Array<LO> lidsRemoteUnion (numRemoteUnion);
//     for (size_type k = 0; k < numRemoteUnion; ++k) {
//       lidsRemoteUnion[k] = as<LO> (numNonRemoteUnion) + as<LO> (k);
//     }

//     // At this point, we have the following:
//     //   - All same-as and permute-to GIDs in the union target Map
//     //   - Permute-to LIDs (in the union target Map) for the union Import
//     //   - Permute-from LIDs (in the source Map) for the union Import
//     //   - Remote GIDs and their LIDs in the union target Map

//     // Now we need:
//     //   - Remote PIDs (in the source Map)
//     //   - To pack (eliminate invalid) and resize remote PID, GIDs,
//     //     and LIDs
//     //   - To sort remote PIDs, and apply the resulting permutation to
//     //     remote GIDs and LIDs
//     //   - To create the Distributor

  }

} // namespace Tpetra

// Explicit instantiation macro.
// Only invoke this when in the Tpetra namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_IMPORT_INSTANT(LO, GO, NODE) \
  \
  template class Import< LO , GO , NODE >;

#endif // TPETRA_IMPORT_DEF_HPP
