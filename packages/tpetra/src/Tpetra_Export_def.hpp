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

#ifndef TPETRA_EXPORT_DEF_HPP
#define TPETRA_EXPORT_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
#  include <Tpetra_Export_decl.hpp>
#endif // DOXYGEN_USE_ONLY

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Import.hpp>
#include <Teuchos_as.hpp>

namespace {
  // Default value of Export's "Debug" parameter.
  const bool tpetraExportDebugDefault = false;
} // namespace (anonymous)

namespace Tpetra {
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    bool debug = tpetraExportDebugDefault;
    if (! plist.is_null ()) {
      try {
        debug = plist->get<bool> ("Debug");
      } catch (Teuchos::Exceptions::InvalidParameter&) {}
    }
    debug_ = debug;
    ExportData_->distributor_.setParameterList (plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraExportDebugDefault)
  {
    using Teuchos::rcp;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    if (! out_.is_null ()) {
      out_->pushTab ();
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor" << endl;
      *out_ << os.str ();
    }
    ExportData_ = rcp (new data_type (source, target, out_));
    Teuchos::Array<GlobalOrdinal> exportGIDs;
    setupSamePermuteExport (exportGIDs);
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: "
         << "setupSamePermuteExport done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupRemote (exportGIDs);
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const RCP<Teuchos::FancyOStream>& out) :
    out_ (out),
    debug_ (tpetraExportDebugDefault)
  {
    using Teuchos::rcp;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    if (! out_.is_null ()) {
      out_->pushTab ();
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor" << endl;
      *out_ << os.str ();
    }
    ExportData_ = rcp (new data_type (source, target, out));
    Teuchos::Array<GlobalOrdinal> exportGIDs;
    setupSamePermuteExport (exportGIDs);
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: "
         << "setupSamePermuteExport done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupRemote (exportGIDs);
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraExportDebugDefault)
  {
    using Teuchos::rcp;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    // Read "Debug" parameter from the input ParameterList.
    bool debug = tpetraExportDebugDefault;
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
      os << myRank << ": Export ctor" << endl;
      *out_ << os.str ();
    }
    ExportData_ = rcp (new data_type (source, target, out_, plist));
    Teuchos::Array<GlobalOrdinal> exportGIDs;
    setupSamePermuteExport (exportGIDs);
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: "
         << "setupSamePermuteExport done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupRemote (exportGIDs);
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraExportDebugDefault)
  {
    using Teuchos::rcp;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    // Read "Debug" parameter from the input ParameterList.
    bool debug = tpetraExportDebugDefault;
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
      os << myRank << ": Export ctor" << endl;
      *out_ << os.str ();
    }
    ExportData_ = rcp (new data_type (source, target, out, plist));
    Teuchos::Array<GlobalOrdinal> exportGIDs;
    setupSamePermuteExport (exportGIDs);
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: "
         << "setupSamePermuteExport done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupRemote (exportGIDs);
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs)
    : ExportData_ (rhs.ExportData_),
      out_ (rhs.out_),
      debug_ (rhs.debug_)
  {
    using std::endl;

    if (! out_.is_null ()) {
      out_->pushTab ();
    }
    if (debug_) {
      std::ostringstream os;
      const int myRank = getSourceMap ()->getComm ()->getRank ();
      os << myRank << ": Export copy ctor (done)" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Import<LocalOrdinal,GlobalOrdinal,Node>& importer)
    : out_ (importer.out_)
    , debug_ (importer.debug_)
  {
    if(!importer.ImportData_.is_null())  ExportData_ = importer.ImportData_->reverseClone();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::~Export()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumSameIDs() const {
    return ExportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumPermuteIDs() const {
    return ExportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ExportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ExportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ExportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ExportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ExportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ExportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const int>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getExportPIDs() const {
    return ExportData_->exportPIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ExportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ExportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ExportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>&
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  operator= (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs) {
    if (&rhs != this) {
      ExportData_ = rhs.ExportData_;
    }
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Export<LocalOrdinal,GlobalOrdinal,Node>::
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
          os << "Export Data Members:" << endl;
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

    // It's also helpful for debugging to print the Distributor
    // object.  Epetra_Export::Print() does this, so we can do a
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
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  setupSamePermuteExport (Teuchos::Array<GlobalOrdinal>& exportGIDs)
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
    // numSameIDs_ is for the common case of an Export where all the
    // overlapping GIDs are at the end of the source Map, but
    // otherwise the source and target Maps are the same.  This allows
    // a fast contiguous copy for the initial "same IDs."
    size_type numSameGids = 0;
    for ( ; numSameGids < numGids && rawSrcGids[numSameGids] == rawTgtGids[numSameGids]; ++numSameGids)
      {} // third clause of 'for' does everything
    ExportData_->numSameIDs_ = numSameGids;

    // Compute permuteToLIDs_, permuteFromLIDs_, exportGIDs, and
    // exportLIDs_.  The first two arrays are IDs to be permuted, and
    // the latter two arrays are IDs to sent out ("exported"), called
    // "export" IDs.
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  IDs to send are in the
    // source Map, but not in the target Map.

    exportGIDs.resize (0);
    Array<LO>& permuteToLIDs = ExportData_->permuteToLIDs_;
    Array<LO>& permuteFromLIDs = ExportData_->permuteFromLIDs_;
    Array<LO>& exportLIDs = ExportData_->exportLIDs_;
    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO numSrcLids = as<LO> (numSrcGids);
    // Iterate over the source Map's LIDs, since we only need to do
    // GID -> LID lookups for the target Map.
    for (LO srcLid = numSameGids; srcLid < numSrcLids; ++srcLid) {
      const GO curSrcGid = rawSrcGids[srcLid];
      // getLocalElement() returns LINVALID if the GID isn't in the target Map.
      // This saves us a lookup (which isNodeGlobalElement() would do).
      const LO tgtLid = target.getLocalElement (curSrcGid);
      if (tgtLid != LINVALID) { // if target.isNodeGlobalElement (curSrcGid)
        permuteToLIDs.push_back (tgtLid);
        permuteFromLIDs.push_back (srcLid);
      } else {
        exportGIDs.push_back (curSrcGid);
        exportLIDs.push_back (srcLid);
      }
    }

    // exportLIDs_ is the list of this process' LIDs that it has to
    // send out.  Since this is an Export, and therefore the target
    // Map is nonoverlapping, we know that each export LID only needs
    // to be sent to one process.  However, the source Map may be
    // overlapping, so multiple processes might send to the same LID
    // on a receiving process.

    TPETRA_ABUSE_WARNING(
      getNumExportIDs() > 0 && ! source.isDistributed(),
      std::runtime_error,
      "::setupSamePermuteExport(): Source has export LIDs but Source is not "
      "distributed globally." << std::endl
      << "Exporting to a submap of the target map.");

    // Compute exportPIDs_ ("outgoing" process IDs).
    //
    // For each GID in exportGIDs (GIDs to which this process must
    // send), find its corresponding owning process (a.k.a. "image")
    // ID in the target Map.  Store these process IDs in
    // exportPIDs_.  These are the process IDs to which the Export
    // needs to send data.
    //
    // We only need to do this if the source Map is distributed;
    // otherwise, the Export doesn't have to perform any
    // communication.
    if (source.isDistributed ()) {
      ExportData_->exportPIDs_.resize(exportGIDs.size ());
      // This call will assign any GID in the target Map with no
      // corresponding process ID a fake process ID of -1.  We'll use
      // this below to remove exports for processses that don't exist.
      const LookupStatus lookup =
        target.getRemoteIndexList (exportGIDs(),
                                   ExportData_->exportPIDs_ ());
      TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
        "::setupSamePermuteExport(): The source Map has GIDs not found "
        "in the target Map.");

      // Get rid of process IDs not in the target Map.  This prevents
      // exporting to GIDs which don't belong to any process in the
      // target Map.
      typedef typename ArrayRCP<int>::difference_type size_type;
      if (lookup == IDNotPresent) {
        const size_type numInvalidExports =
          std::count_if (ExportData_->exportPIDs_().begin(),
                         ExportData_->exportPIDs_().end(),
                         std::bind1st (std::equal_to<int>(), -1));

        // count number of valid and total number of exports
        const size_type totalNumExports = ExportData_->exportPIDs_.size();
        if (numInvalidExports == totalNumExports) {
          // all exports are invalid; we have no exports; we can delete all exports
          exportGIDs.resize(0);
          ExportData_->exportLIDs_.resize(0);
          ExportData_->exportPIDs_.resize(0);
        }
        else {
          // some exports are valid; we need to keep the valid exports
          // pack and resize
          size_type numValidExports = 0;
          for (size_type e = 0; e < totalNumExports; ++e) {
            if (ExportData_->exportPIDs_[e] != -1) {
              exportGIDs[numValidExports]               = exportGIDs[e];
              ExportData_->exportLIDs_[numValidExports] = ExportData_->exportLIDs_[e];
              ExportData_->exportPIDs_[numValidExports] = ExportData_->exportPIDs_[e];
              ++numValidExports;
            }
          }
          exportGIDs.resize (numValidExports);
          ExportData_->exportLIDs_.resize (numValidExports);
          ExportData_->exportPIDs_.resize (numValidExports);
        }
      }
    }
  } // setupSamePermuteExport()

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::setupRemote(Teuchos::Array<GlobalOrdinal> & exportGIDs)
  {
    using std::endl;
    const Map<LocalOrdinal,GlobalOrdinal,Node>& target = * (getTargetMap ());
    const int myRank = target.getComm ()->getRank ();

    if (! out_.is_null ()) {
      out_->pushTab ();
    }
    if (debug_) {
      std::ostringstream os;
      os << myRank << ": Export::setupRemote" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->pushTab ();
    }

    // Sort exportPIDs_ in ascending order, and apply the same
    // permutation to exportGIDs_ and exportLIDs_.  This ensures that
    // exportPIDs_[i], exportGIDs_[i], and exportLIDs_[i] all
    // refer to the same thing.
    sort3 (ExportData_->exportPIDs_.begin(),
           ExportData_->exportPIDs_.end(),
           exportGIDs.begin(),
           ExportData_->exportLIDs_.begin());

    if (debug_) {
      std::ostringstream os;
      os << myRank << ": Export::setupRemote: Calling createFromSends" << endl;
      *out_ << os.str ();
    }

    // Construct the list of entries that calling image needs to send
    // as a result of everyone asking for what it needs to receive.
    //
    // mfh 05 Jan 2012: I understand the above comment as follows:
    // Construct the communication plan from the list of image IDs to
    // which we need to send.
    size_t numRemoteIDs;
    numRemoteIDs = ExportData_->distributor_.createFromSends (ExportData_->exportPIDs_ ());

    if (debug_) {
      std::ostringstream os;
      os << myRank << ": Export::setupRemote: Calling doPostsAndWaits" << endl;
      *out_ << os.str ();
    }

    // Use the communication plan with ExportGIDs to find out who is
    // sending to us and get the proper ordering of GIDs for incoming
    // remote entries (these will be converted to LIDs when done).
    Array<GlobalOrdinal> remoteGIDs (numRemoteIDs);
    ExportData_->distributor_.doPostsAndWaits (exportGIDs().getConst (), 1, remoteGIDs());

    // Remote (incoming) IDs come in as GIDs; convert to LIDs.  LIDs
    // tell this process where to store the incoming remote data.
    ExportData_->remoteLIDs_.resize (numRemoteIDs);
    {
      typename Array<GlobalOrdinal>::const_iterator i = remoteGIDs.begin();
      typename Array<LocalOrdinal>::iterator        j = ExportData_->remoteLIDs_.begin();
      while (i != remoteGIDs.end()) {
        *j++ = target.getLocalElement(*i++);
      }
    }

    if (! out_.is_null ()) {
      out_->popTab ();
    }
    if (debug_) {
      std::ostringstream os;
      os << myRank << ": Export::setupRemote: done" << endl;
      *out_ << os.str ();
    }
    if (! out_.is_null ()) {
      out_->popTab ();
    }
  }
} // namespace Tpetra

// Explicit instantiation macro.
// Only invoke this when in the Tpetra namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_EXPORT_INSTANT(LO, GO, NODE) \
  \
  template class Export< LO , GO , NODE >;

#endif // TPETRA_EXPORT_DEF_HPP
