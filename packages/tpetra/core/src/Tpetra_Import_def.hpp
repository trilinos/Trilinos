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

#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Import_Util.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Teuchos_as.hpp>
#ifdef HAVE_TPETRA_MMM_TIMINGS
#include <Teuchos_TimeMonitor.hpp>
#endif


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
  init (const Teuchos::RCP<const map_type>& source,
        const Teuchos::RCP<const map_type>& target,
        bool useRemotePIDs,
        Teuchos::Array<int> & remotePIDs,
        const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::Ptr;
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

    Array<GlobalOrdinal> remoteGIDs;
    setupSamePermuteRemote (remoteGIDs);
    if (debug_) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Import ctor: "
         << "setupSamePermuteRemote done" << endl;
      *out_ << os.str ();
    }
    if (source->isDistributed ()) {
      setupExport (remoteGIDs,useRemotePIDs,remotePIDs);
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
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraImportDebugDefault)
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, Teuchos::null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out) :
    out_ (out),
    debug_ (tpetraImportDebugDefault)
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, Teuchos::null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr))),
    debug_ (tpetraImportDebugDefault)
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    out_ (out),
    debug_ (tpetraImportDebugDefault)
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          Teuchos::Array<int> & remotePIDs) :
    debug_ (tpetraImportDebugDefault)
  {
    init (source, target, true, remotePIDs, Teuchos::null);
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

  // This is the "createExpert" version of the constructor

   template <class LocalOrdinal, class GlobalOrdinal, class Node>
   Import<LocalOrdinal,GlobalOrdinal,Node>::
   Import(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& target,
          Teuchos::Array<int>& userRemotePIDs,
          Teuchos::Array<GlobalOrdinal>& remoteGIDs,
          const Teuchos::ArrayView<const LocalOrdinal> & userExportLIDs,
          const Teuchos::ArrayView<const int> & userExportPIDs,
          const bool useRemotePIDGID,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<Teuchos::FancyOStream>& out) :
     out_ (out.is_null () ?
           Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) : out)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::rcp;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Teuchos::Array<int>::size_type size_type;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    out_->pushTab ();

    ArrayView<const GO> sourceGIDs = source->getNodeElementList ();
    ArrayView<const GO> targetGIDs = target->getNodeElementList ();
    const size_type numSrcGids = sourceGIDs.size ();
    const size_type numTgtGids = targetGIDs.size ();
    const size_type numGids = std::min (numSrcGids, numTgtGids);

    size_type numSameGids = 0;
    for ( ; numSameGids < numGids && sourceGIDs[numSameGids] == targetGIDs[numSameGids]; ++numSameGids)
      {}

    // Read "Debug" parameter from the input ParameterList.
    bool debug = tpetraImportDebugDefault;
    if (! plist.is_null ()) {
      try {
        debug = plist->get<bool> ("Debug");
      } catch (Teuchos::Exceptions::InvalidParameter&) {}
    }
    debug_ = debug;

    if (debug_ && ! out_.is_null ()) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": constructExpert " << std::endl;
      *out_ << os.str ();
    }
    ImportData_ = rcp (new data_type (source, target, out_, plist));
    ImportData_->numSameIDs_ = numSameGids;

    Array<LO>& permuteToLIDs = ImportData_->permuteToLIDs_;
    Array<LO>& permuteFromLIDs = ImportData_->permuteFromLIDs_;
    Array<LO>& remoteLIDs = ImportData_->remoteLIDs_;
    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO numTgtLids = as<LO> (numTgtGids);

    if(!useRemotePIDGID) {
      remoteGIDs.clear();
      remoteLIDs.clear();
    }

    for (LO tgtLid = numSameGids; tgtLid < numTgtLids; ++tgtLid) {
      const GO curTargetGid = targetGIDs[tgtLid];
      // getLocalElement() returns LINVALID if the GID isn't in the source Map.
      const LO srcLid = source->getLocalElement (curTargetGid);
      if (srcLid != LINVALID) {
        permuteToLIDs.push_back (tgtLid);
        permuteFromLIDs.push_back (srcLid);
      } else {
        if(!useRemotePIDGID) {
          remoteGIDs.push_back (curTargetGid);
          remoteLIDs.push_back (tgtLid);
        }
      }
    }

    TPETRA_ABUSE_WARNING(
      getNumRemoteIDs() > 0 && ! source->isDistributed(),
      std::runtime_error,
      "::constructExpert(): Target has remote LIDs but Source is not "
      "distributed globally." << std::endl
      << "Importing to a submap of the target map.");

    Array<int> remotePIDs;
    remotePIDs.resize (remoteGIDs.size (),0);
    LookupStatus lookup = AllIDsPresent;

    ArrayView<GO> remoteGIDsView = remoteGIDs ();
    lookup = source->getRemoteIndexList (remoteGIDsView, remotePIDs ());
    remoteGIDsView = remoteGIDs ();

    Array<int>& remoteProcIDs = (useRemotePIDGID) ? userRemotePIDs : remotePIDs;

    TEUCHOS_TEST_FOR_EXCEPTION( lookup == IDNotPresent, std::runtime_error,
      "Import::Import createExpert: the source Map wasn't able to figure out which process "
      "owns one or more of the GIDs in the list of remote GIDs.  This probably "
      "means that there is at least one GID owned by some process in the target"
      " Map which is not owned by any process in the source Map.  (That is, the"
      " source and target Maps do not contain the same set of GIDs globally.)");

    // Sort remoteProcIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs and remoteLIDs_.  This ensures that
    // remoteProcIDs[i], remoteGIDs[i], and remoteLIDs_[i] all refer
    // to the same thing.

    TEUCHOS_TEST_FOR_EXCEPTION( !(remoteProcIDs.size() == remoteGIDsView.size() &&remoteGIDsView.size() == remoteLIDs.size()), std::runtime_error,
                               "Import::Import createExpert version: Size miss match on RemoteProcIDs, remoteGIDsView and remoteLIDs Array's to sort3. This will produce produce an error, aborting ");

    sort3 (remoteProcIDs.begin (),
           remoteProcIDs.end (),
           remoteGIDsView.begin (),
           remoteLIDs.begin ());

    ImportData_->remoteLIDs_ = remoteLIDs;
    ImportData_->distributor_ =  Distributor (source->getComm(),this->out_);
    ImportData_->exportPIDs_ = Teuchos::Array<int>(userExportPIDs.size(),0);
    ImportData_->exportLIDs_ = Teuchos::Array<int>(userExportPIDs.size(),0);

    bool locallyComplete = true;
    for(size_type i=0; i<userExportPIDs.size(); i++)  {
      if (userExportPIDs[i] == -1) {
        locallyComplete = false;
      }
      ImportData_->exportPIDs_[i] = userExportPIDs[i];
      ImportData_->exportLIDs_[i] = userExportLIDs[i];
    }
    ImportData_->isLocallyComplete_ = locallyComplete;

    ImportData_->distributor_.createFromSendsAndRecvs(ImportData_->exportPIDs_,remoteProcIDs);

  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target,
          const size_t numSameIDs,
          Teuchos::Array<LocalOrdinal>& permuteToLIDs,
          Teuchos::Array<LocalOrdinal>& permuteFromLIDs,
          Teuchos::Array<LocalOrdinal>& remoteLIDs,
          Teuchos::Array<LocalOrdinal>& exportLIDs,
          Teuchos::Array<int>& exportPIDs,
          Distributor& distributor,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
     out_ (out.is_null () ? Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) : out),
    debug_ (tpetraImportDebugDefault)
  {
    using Teuchos::null;
    using Teuchos::Ptr;
    using Teuchos::rcp;
    using std::cerr;
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
    if (debug_ && ! out_.is_null ()) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Import expert ctor" << endl;
      *out_ << os.str ();
    }
    ImportData_ = rcp (new data_type (source, target, out_, plist));

    bool locallyComplete = true;
    for (Teuchos::Array<int>::size_type i = 0; i < exportPIDs.size (); ++i) {
      if (exportPIDs[i] == -1) {
        locallyComplete = false;
      }
    }
    ImportData_->isLocallyComplete_ = locallyComplete;

    ImportData_->numSameIDs_ = numSameIDs;
    ImportData_->permuteToLIDs_.swap (permuteToLIDs);
    ImportData_->permuteFromLIDs_.swap (permuteFromLIDs);
    ImportData_->remoteLIDs_.swap (remoteLIDs);
    ImportData_->distributor_.swap (distributor);
    ImportData_->exportLIDs_.swap (exportLIDs);
    ImportData_->exportPIDs_.swap (exportPIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& sourceMap,
          const GlobalOrdinal targetMapGlobalIndices[],
          const int /* processFromWhichToReceiveEachOfThoseGlobalIndices */ [],
          const std::size_t numTargetMapGlobalIndices,
          const bool /* mayReorderTargetMapIndicesLocally */ ,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<Teuchos::FancyOStream>& out) :
    Import (sourceMap,
            Teuchos::rcp (new map_type (Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid (),
                                        targetMapGlobalIndices,
                                        numTargetMapGlobalIndices,
                                        sourceMap->getIndexBase (),
                                        sourceMap->getComm ())),
            plist,
            out)
  {
    // stub implementation for now; just dispatch to the usual Import ctor
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
  Teuchos::ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ImportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ImportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ImportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const int>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportPIDs() const {
    return ImportData_->exportPIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename Import<LocalOrdinal,GlobalOrdinal,Node>::map_type>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename Import<LocalOrdinal,GlobalOrdinal,Node>::map_type>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Import<LocalOrdinal,GlobalOrdinal,Node>::isLocallyComplete () const {
    return ImportData_->isLocallyComplete_;
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
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    // Call the base class' method.  It does all the work.
    this->describeImpl (out, "Tpetra::Import", verbLevel);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::
  print (std::ostream& os) const
  {
    auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));
    // "Print" traditionally meant "everything."
    this->describe (*out, Teuchos::VERB_EXTREME);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setupSamePermuteRemote (Teuchos::Array<GlobalOrdinal>& remoteGIDs)
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
    const map_type& source = * (getSourceMap ());
    const map_type& target = * (getTargetMap ());
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

    // Compute permuteToLIDs_, permuteFromLIDs_, remoteGIDs, and
    // remoteLIDs_.  The first two arrays are IDs to be permuted, and
    // the latter two arrays are IDs to be received ("imported"),
    // called "remote" IDs.
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  IDs to receive are in the
    // target Map, but not the source Map.

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
        remoteGIDs.push_back (curTargetGid);
        remoteLIDs.push_back (tgtLid);
      }
    }

    if (remoteLIDs.size () != 0 && ! source.isDistributed ()) {
      // This Import has remote LIDs, meaning that the target Map has
      // entries on this process that are not in the source Map on
      // this process.  However, the source Map is not distributed
      // globally.  This implies that this Import is not locally
      // complete on this process.
      ImportData_->isLocallyComplete_ = false;
      // mfh 12 Sep 2016: I disagree that this is "abuse"; it may be
      // correct behavior, depending on the circumstances.
      TPETRA_ABUSE_WARNING
        (true, std::runtime_error, "::setupSamePermuteRemote(): Target has "
         "remote LIDs but Source is not distributed globally.  Importing to a "
         "submap of the target map.");
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::
  setupExport (Teuchos::Array<GlobalOrdinal>& remoteGIDs,
               bool useRemotePIDs,
               Teuchos::Array<int>& userRemotePIDs)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Array<int>::difference_type size_type;
    const char tfecfFuncName[] = "setupExport: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (getSourceMap ().is_null (), std::logic_error, "Source Map is null.  "
       "Please report this bug to the Tpetra developers.");
    const map_type& source = * (getSourceMap ());

    Teuchos::OSTab tab (out_);

    // if (debug_ && ! out_.is_null ()) {
    //   std::ostringstream os;
    //   const int myRank = source.getComm ()->getRank ();
    //   os << myRank << ": Import::setupExport:" << endl;
    // }

    // Sanity checks
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! useRemotePIDs && (userRemotePIDs.size() > 0), std::invalid_argument,
       "remotePIDs are non-empty but their use has not been requested.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (userRemotePIDs.size () > 0 && remoteGIDs.size () != userRemotePIDs.size (),
       std::invalid_argument, "remotePIDs must either be of size zero or match "
       "the size of remoteGIDs.");

    // For each entry remoteGIDs[i], remoteProcIDs[i] will contain
    // the process ID of the process that owns that GID.
    ArrayView<GO> remoteGIDsView = remoteGIDs ();
    ArrayView<int> remoteProcIDsView;

    // lookup == IDNotPresent means that the source Map wasn't able to
    // figure out to which processes one or more of the GIDs in the
    // given list of remote GIDs belong.
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
    Array<int> newRemotePIDs;
    LookupStatus lookup = AllIDsPresent;

    if (! useRemotePIDs) {
      newRemotePIDs.resize (remoteGIDsView.size ());
      if (debug_ && ! out_.is_null ()) {
        std::ostringstream os;
        const int myRank = source.getComm ()->getRank ();
        os << myRank << ": Import::setupExport: about to call "
          "getRemoteIndexList on source Map" << endl;
        *out_ << os.str ();
      }
      lookup = source.getRemoteIndexList (remoteGIDsView, newRemotePIDs ());
    }
    Array<int>& remoteProcIDs = useRemotePIDs ? userRemotePIDs : newRemotePIDs;

    if (lookup == IDNotPresent) {
      // There is at least one GID owned by the calling process in the
      // target Map, which is not owned by any process in the source
      // Map.
      ImportData_->isLocallyComplete_ = false;

      // mfh 12 Sep 2016: I disagree that this is "abuse"; it may be
      // correct behavior, depending on the circumstances.
      TPETRA_ABUSE_WARNING
        (true, std::runtime_error, "::setupExport(): the source Map wasn't "
         "able to figure out which process owns one or more of the GIDs in the "
         "list of remote GIDs.  This probably means that there is at least one "
         "GID owned by some process in the target Map which is not owned by any"
         " process in the source Map.  (That is, the source and target Maps do "
         "not contain the same set of GIDs globally.)");

      // Ignore remote GIDs that aren't owned by any process in the
      // source Map.  getRemoteIndexList() gives each of these a
      // process ID of -1.

      const size_type numInvalidRemote =
        std::count_if (remoteProcIDs.begin (), remoteProcIDs.end (),
                       std::bind1st (std::equal_to<int> (), -1));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numInvalidRemote == 0, std::logic_error, "Calling getRemoteIndexList "
         "on the source Map returned IDNotPresent, but none of the returned "
         "\"remote\" process ranks are -1.  Please report this bug to the "
         "Tpetra developers.");

      // If all of them are invalid, we can delete the whole array.
      const size_type totalNumRemote = getNumRemoteIDs ();
      if (numInvalidRemote == totalNumRemote) {
        // all remotes are invalid; we have no remotes; we can delete the remotes
        remoteProcIDs.clear ();
        remoteGIDs.clear (); // This invalidates the view remoteGIDsView
        ImportData_->remoteLIDs_.clear();
      }
      else {
        // Some remotes are valid; we need to keep the valid ones.
        // Pack and resize remoteProcIDs, remoteGIDs, and remoteLIDs_.
        size_type numValidRemote = 0;
#ifdef HAVE_TPETRA_DEBUG
        ArrayView<GlobalOrdinal> remoteGIDsPtr = remoteGIDsView;
#else
        GlobalOrdinal* const remoteGIDsPtr = remoteGIDsView.getRawPtr ();
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
        remoteGIDs.resize (numValidRemote);
        ImportData_->remoteLIDs_.resize (numValidRemote);
      }
      // Revalidate the view after clear or resize.
      remoteGIDsView = remoteGIDs ();
    }

    // Sort remoteProcIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs and remoteLIDs_.  This ensures that
    // remoteProcIDs[i], remoteGIDs[i], and remoteLIDs_[i] all refer
    // to the same thing.

    sort3 (remoteProcIDs.begin (),
           remoteProcIDs.end (),
           remoteGIDsView.begin (),
           ImportData_->remoteLIDs_.begin ());

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDs and remoteProcIDs_ are input;
    // exportGIDs and exportProcIDs_ are output arrays which are
    // allocated by createFromRecvs().
    Array<GO> exportGIDs;
    ImportData_->distributor_.createFromRecvs (remoteGIDsView ().getConst (),
                                               remoteProcIDs, exportGIDs,
                                               ImportData_->exportPIDs_);
    // if (debug_ && ! out_.is_null ()) {
    //   std::ostringstream os;
    //   const int myRank = source.getComm ()->getRank ();
    //   os << myRank << ": Import::setupExport: Getting LIDs" << endl;
    //   *out_ << os.str ();
    // }

    // Find the LIDs corresponding to the (outgoing) GIDs in
    // exportGIDs.  For sparse matrix-vector multiply, this tells the
    // calling process how to index into the source vector to get the
    // elements which it needs to send.
    //
    // NOTE (mfh 03 Mar 2014) This is now a candidate for a
    // thread-parallel kernel, but only if using the new thread-safe
    // Map implementation.
    const size_type numExportIDs = exportGIDs.size ();
    if (numExportIDs > 0) {
      ImportData_->exportLIDs_.resize (numExportIDs);
      ArrayView<const GO> expGIDs = exportGIDs ();
      ArrayView<LO> expLIDs = ImportData_->exportLIDs_ ();
      for (size_type k = 0; k < numExportIDs; ++k) {
        expLIDs[k] = source.getLocalElement (expGIDs[k]);
      }
    }

    if (debug_ && ! out_.is_null ()) {
      std::ostringstream os;
      const int myRank = source.getComm ()->getRank ();
      os << myRank << ": Import::setupExport: done" << endl;
      *out_ << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  findUnionTargetGIDs(Teuchos::Array<GlobalOrdinal>& unionTgtGIDs,
                      Teuchos::Array<std::pair<int,GlobalOrdinal>>& remotePGIDs,
                      typename Teuchos::Array<GlobalOrdinal>::size_type& numSameGIDs,
                      typename Teuchos::Array<GlobalOrdinal>::size_type& numPermuteGIDs,
                      typename Teuchos::Array<GlobalOrdinal>::size_type& numRemoteGIDs,
                      const Teuchos::ArrayView<const GlobalOrdinal>& sameGIDs1,
                      const Teuchos::ArrayView<const GlobalOrdinal>& sameGIDs2,
                      Teuchos::Array<GlobalOrdinal>& permuteGIDs1,
                      Teuchos::Array<GlobalOrdinal>& permuteGIDs2,
                      Teuchos::Array<GlobalOrdinal>& remoteGIDs1,
                      Teuchos::Array<GlobalOrdinal>& remoteGIDs2,
                      Teuchos::Array<int>& remotePIDs1,
                      Teuchos::Array<int>& remotePIDs2) const
  {

    typedef GlobalOrdinal GO;
    typedef typename Teuchos::Array<GO>::size_type size_type;

    const size_type numSameGIDs1 = sameGIDs1.size();
    const size_type numSameGIDs2 = sameGIDs2.size();

    // Sort the permute GIDs
    std::sort(permuteGIDs1.begin(), permuteGIDs1.end());
    std::sort(permuteGIDs2.begin(), permuteGIDs2.end());

    // Get the union of the two target maps
    // Reserve the maximum possible size to guard against reallocations from
    // push_back operations.
    unionTgtGIDs.reserve(numSameGIDs1 + numSameGIDs2 +
                         permuteGIDs1.size() + permuteGIDs2.size() +
                         remoteGIDs1.size() + remoteGIDs2.size());

    // Copy the same GIDs to unionTgtGIDs.  Cases for numSameGIDs1 !=
    // numSameGIDs2 must be treated separately.
    typename Teuchos::Array<GO>::iterator permuteGIDs1_end;
    typename Teuchos::Array<GO>::iterator permuteGIDs2_end;
    if (numSameGIDs2 > numSameGIDs1) {

      numSameGIDs = numSameGIDs2;
      permuteGIDs2_end = permuteGIDs2.end();

      // Copy the same GIDs from tgtGIDs to the union
      std::copy(sameGIDs2.begin(), sameGIDs2.end(), std::back_inserter(unionTgtGIDs));

      // Remove GIDs from permuteGIDs1 that have already been copied in to unionTgtGIDs
      // set_difference allows the last (output) argument to alias the first.
      permuteGIDs1_end = std::set_difference(permuteGIDs1.begin(), permuteGIDs1.end(),
          unionTgtGIDs.begin()+numSameGIDs1, unionTgtGIDs.end(),
          permuteGIDs1.begin());

    } else {

      numSameGIDs = numSameGIDs1;
      permuteGIDs1_end = permuteGIDs1.end();

      // Copy the same GIDs from tgtGIDs to the union
      std::copy(sameGIDs1.begin(), sameGIDs1.end(), std::back_inserter(unionTgtGIDs));

      // Remove GIDs from permuteGIDs2 that have already been copied in to unionTgtGIDs
      // set_difference allows the last (output) argument to alias the first.
      permuteGIDs2_end = std::set_difference(permuteGIDs2.begin(), permuteGIDs2.end(),
          unionTgtGIDs.begin()+numSameGIDs2, unionTgtGIDs.end(),
          permuteGIDs2.begin());

    }

    // Get the union of the permute GIDs and push it back on unionTgtGIDs
    std::set_union(permuteGIDs1.begin(), permuteGIDs1_end,
                   permuteGIDs2.begin(), permuteGIDs2_end,
                   std::back_inserter(unionTgtGIDs));

    // Sort the PID,GID pairs and find the unique set
    Teuchos::Array<std::pair<int,GO>> remotePGIDs1(remoteGIDs1.size());
    for (size_type k=0; k<remoteGIDs1.size(); k++)
      remotePGIDs1[k] = std::make_pair(remotePIDs1[k], remoteGIDs1[k]);
    std::sort(remotePGIDs1.begin(), remotePGIDs1.end());

    Teuchos::Array<std::pair<int,GO>> remotePGIDs2(remoteGIDs2.size());
    for (size_type k=0; k<remoteGIDs2.size(); k++)
      remotePGIDs2[k] = std::make_pair(remotePIDs2[k], remoteGIDs2[k]);
    std::sort(remotePGIDs2.begin(), remotePGIDs2.end());

    remotePGIDs.reserve(remotePGIDs1.size()+remotePGIDs2.size());
    std::merge(remotePGIDs1.begin(), remotePGIDs1.end(),
               remotePGIDs2.begin(), remotePGIDs2.end(),
               std::back_inserter(remotePGIDs));
    auto it = std::unique(remotePGIDs.begin(), remotePGIDs.end());
    remotePGIDs.resize(std::distance(remotePGIDs.begin(), it));

    // Finally, insert remote GIDs
    const size_type oldSize = unionTgtGIDs.size();
    unionTgtGIDs.resize(oldSize+remotePGIDs.size());
    for (size_type start=oldSize, k=0; k<remotePGIDs.size(); k++)
      unionTgtGIDs[start+k] = remotePGIDs[k].second;

    // Compute output only quantities
    numRemoteGIDs = remotePGIDs.size();
    numPermuteGIDs = unionTgtGIDs.size() - numSameGIDs - numRemoteGIDs;

    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setUnion (const Import<LocalOrdinal, GlobalOrdinal, Node>& rhs) const
  {
    typedef Tpetra::global_size_t GST;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Import<LO, GO, Node> import_type;
    typedef typename Array<GO>::size_type size_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    using Teuchos::TimeMonitor;
    std::string label = std::string("Tpetra::Import::setUnion");
    RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(label)));
    label = "Tpetra::Import::setUnion : Union GIDs";
    RCP<TimeMonitor> MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(label)));
#endif

    RCP<const map_type> srcMap = this->getSourceMap ();
    RCP<const map_type> tgtMap1 = this->getTargetMap ();
    RCP<const map_type> tgtMap2 = rhs.getTargetMap ();
    RCP<const Comm<int> > comm = srcMap->getComm ();

    const bool debug = Details::Behavior::debug("Tpetra::Import::setUnion");

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          ! srcMap->isSameAs (* (rhs.getSourceMap ())), std::invalid_argument,
          "Tpetra::Import::setUnion: The source Map of the input Import must be the "
          "same as (in the sense of Map::isSameAs) the source Map of this Import.");
      TEUCHOS_TEST_FOR_EXCEPTION(
          ! Details::congruent (* (tgtMap1->getComm ()), * (tgtMap2->getComm ())),
          std::invalid_argument, "Tpetra::Import::setUnion: "
          "The target Maps must have congruent communicators.");
    }

    // It's probably worth the one all-reduce to check whether the two
    // Maps are the same.  If so, we can just return a copy of *this.
    // isSameAs() bypasses the all-reduce if the pointers are equal.
    if (tgtMap1->isSameAs (*tgtMap2)) {
      return rcp (new import_type (*this));
    }

    // Alas, the two target Maps are not the same.  That means we have
    // to compute their union, and the union Import object.

    // Get the same GIDs (same GIDs are a subview of the first numSame target
    // GIDs)
    const size_type numSameGIDs1 = this->getNumSameIDs();
    ArrayView<const GO> sameGIDs1 = (tgtMap1->getNodeElementList())(0,numSameGIDs1);

    const size_type numSameGIDs2 = rhs.getNumSameIDs();
    ArrayView<const GO> sameGIDs2 = (tgtMap2->getNodeElementList())(0,numSameGIDs2);

    // Get permute GIDs
    ArrayView<const LO> permuteToLIDs1 = this->getPermuteToLIDs();
    Array<GO> permuteGIDs1(permuteToLIDs1.size());
    for (size_type k=0; k<permuteGIDs1.size(); k++)
      permuteGIDs1[k] = tgtMap1->getGlobalElement(permuteToLIDs1[k]);

    ArrayView<const LO> permuteToLIDs2 = rhs.getPermuteToLIDs();
    Array<GO> permuteGIDs2(permuteToLIDs2.size());
    for (size_type k=0; k<permuteGIDs2.size(); k++)
      permuteGIDs2[k] = tgtMap2->getGlobalElement(permuteToLIDs2[k]);

    // Get remote GIDs
    ArrayView<const LO> remoteLIDs1 = this->getRemoteLIDs();
    Array<GO> remoteGIDs1(remoteLIDs1.size());
    for (size_type k=0; k<remoteLIDs1.size(); k++)
      remoteGIDs1[k] = this->getTargetMap()->getGlobalElement(remoteLIDs1[k]);

    ArrayView<const LO> remoteLIDs2 = rhs.getRemoteLIDs();
    Array<GO> remoteGIDs2(remoteLIDs2.size());
    for (size_type k=0; k<remoteLIDs2.size(); k++)
      remoteGIDs2[k] = rhs.getTargetMap()->getGlobalElement(remoteLIDs2[k]);

    // Get remote PIDs
    Array<int> remotePIDs1;
    Tpetra::Import_Util::getRemotePIDs(*this, remotePIDs1);

    Array<int> remotePIDs2;
    Tpetra::Import_Util::getRemotePIDs(rhs, remotePIDs2);

    // Get the union of the target GIDs
    Array<GO> unionTgtGIDs;
    Array<std::pair<int,GO>> remotePGIDs;
    size_type numSameIDsUnion, numPermuteIDsUnion, numRemoteIDsUnion;

    findUnionTargetGIDs(unionTgtGIDs, remotePGIDs,
                        numSameIDsUnion, numPermuteIDsUnion, numRemoteIDsUnion,
                        sameGIDs1, sameGIDs2, permuteGIDs1, permuteGIDs2,
                        remoteGIDs1, remoteGIDs2, remotePIDs1, remotePIDs2);

    // Extract GIDs and compute LIDS, PIDs for the remotes in the union
    Array<LO> remoteLIDsUnion(numRemoteIDsUnion);
    Array<GO> remoteGIDsUnion(numRemoteIDsUnion);
    Array<int> remotePIDsUnion(numRemoteIDsUnion);
    const size_type unionRemoteIDsStart = numSameIDsUnion + numPermuteIDsUnion;
    for (size_type k = 0; k < numRemoteIDsUnion; ++k) {
      remoteLIDsUnion[k] = unionRemoteIDsStart + k;
      remotePIDsUnion[k] = remotePGIDs[k].first;
      remoteGIDsUnion[k] = remotePGIDs[k].second;
    }

    // Compute the permute-to LIDs (in the union target Map).
    // Convert the permute GIDs to permute-from LIDs in the source Map.
    Array<LO> permuteToLIDsUnion(numPermuteIDsUnion);
    Array<LO> permuteFromLIDsUnion(numPermuteIDsUnion);
    for (size_type k = 0; k < numPermuteIDsUnion; ++k) {
      size_type idx = numSameIDsUnion + k;
      permuteToLIDsUnion[k] = static_cast<LO>(idx);
      permuteFromLIDsUnion[k] = srcMap->getLocalElement(unionTgtGIDs[idx]);
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM2 = Teuchos::null;
    label = "Tpetra::Import::setUnion : Construct Tgt Map";
    MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(label)));
#endif

    // Create the union target Map.
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBaseUnion = std::min(tgtMap1->getIndexBase(), tgtMap2->getIndexBase());
    RCP<const map_type> unionTgtMap =
      rcp(new map_type(INVALID, unionTgtGIDs(), indexBaseUnion, comm, srcMap->getNode()));

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM2 = Teuchos::null;
    label = "Tpetra::Import::setUnion : Export GIDs";
    MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(label)));
#endif

    // Thus far, we have computed the following in the union Import:
    //   - numSameIDs
    //   - numPermuteIDs and permuteFromLIDs
    //   - numRemoteIDs, remoteGIDs, remoteLIDs, and remotePIDs
    //
    // Now it's time to compute the export IDs and initialize the
    // Distributor.

    Array<GO> exportGIDsUnion;
    Array<LO> exportLIDsUnion;
    Array<int> exportPIDsUnion;
    Distributor distributor (comm, this->out_);

#ifdef TPETRA_IMPORT_SETUNION_USE_CREATE_FROM_SENDS
    // Compute the export IDs without communication, by merging the
    // lists of (export LID, export PID) pairs from the two input
    // Import objects.  The export LIDs in both input Import objects
    // are LIDs in the source Map.  Then, use the export PIDs to
    // initialize the Distributor via createFromSends.

    // const size_type numExportIDs1 = this->getNumExportIDs ();
    ArrayView<const LO> exportLIDs1 = this->getExportLIDs ();
    ArrayView<const LO> exportPIDs1 = this->getExportPIDs ();

    // const size_type numExportIDs2 = rhs.getNumExportIDs ();
    ArrayView<const LO> exportLIDs2 = rhs.getExportLIDs ();
    ArrayView<const LO> exportPIDs2 = rhs.getExportPIDs ();

    // We have to keep the export LIDs in PID-sorted order, then merge
    // them.  So, first key-value merge (LID,PID) pairs, treating PIDs
    // as values, merging values by replacement.  Then, sort the
    // (LID,PID) pairs again by PID.

    // Sort (LID,PID) pairs by LID for the later merge, and make
    // each sequence unique by LID.
    Array<LO> exportLIDs1Copy (exportLIDs1.begin (), exportLIDs1.end ());
    Array<int> exportPIDs1Copy (exportLIDs1.begin (), exportLIDs1.end ());
    sort2 (exportLIDs1Copy.begin (), exportLIDs1Copy.end (),
           exportPIDs1Copy.begin ());
    typename ArrayView<LO>::iterator exportLIDs1_end = exportLIDs1Copy.end ();
    typename ArrayView<LO>::iterator exportPIDs1_end = exportPIDs1Copy.end ();
    merge2 (exportLIDs1_end, exportPIDs1_end,
            exportLIDs1Copy.begin (), exportLIDs1_end,
            exportPIDs1Copy.begin (), exportPIDs1_end,
            project1st<LO, LO> ());

    Array<LO> exportLIDs2Copy (exportLIDs2.begin (), exportLIDs2.end ());
    Array<int> exportPIDs2Copy (exportLIDs2.begin (), exportLIDs2.end ());
    sort2 (exportLIDs2Copy.begin (), exportLIDs2Copy.end (),
           exportPIDs2Copy.begin ());
    typename ArrayView<LO>::iterator exportLIDs2_end = exportLIDs2Copy.end ();
    typename ArrayView<LO>::iterator exportPIDs2_end = exportPIDs2Copy.end ();
    merge2 (exportLIDs2_end, exportPIDs2_end,
            exportLIDs2Copy.begin (), exportLIDs2_end,
            exportPIDs2Copy.begin (), exportPIDs2_end,
            project1st<LO, LO> ());

    // Merge export (LID,PID) pairs.  In this merge operation, the
    // LIDs are the "keys" and the PIDs their "values."  We combine
    // the "values" (PIDs) in the pairs by replacement, rather than
    // by adding them together.
    keyValueMerge (exportLIDs1Copy.begin (), exportLIDs1Copy.end (),
                   exportPIDs1Copy.begin (), exportPIDs1Copy.end (),
                   exportLIDs2Copy.begin (), exportLIDs2Copy.end (),
                   exportPIDs2Copy.begin (), exportPIDs2Copy.end (),
                   std::back_inserter (exportLIDsUnion),
                   std::back_inserter (exportPIDsUnion),
                   project1st<int, int> ());

    // Resort the merged (LID,PID) pairs by PID.
    sort2 (exportPIDsUnion.begin (), exportPIDsUnion.end (),
           exportLIDsUnion.begin ());

    // Initialize the Distributor.  Using createFromSends instead of
    // createFromRecvs avoids the initialization and use of a
    // temporary Distributor object.
    (void) distributor.createFromSends (exportPIDsUnion ().getConst ());
#else // NOT TPETRA_IMPORT_SETUNION_USE_CREATE_FROM_SENDS

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDsUnion and remotePIDsUnion are
    // input; exportGIDsUnion and exportPIDsUnion are output arrays
    // which are allocated by createFromRecvs().
    distributor.createFromRecvs (remoteGIDsUnion().getConst(),
                                 remotePIDsUnion().getConst(),
                                 exportGIDsUnion, exportPIDsUnion);

    // Find the (source Map) LIDs corresponding to the export GIDs.
    const size_type numExportIDsUnion = exportGIDsUnion.size ();
    exportLIDsUnion.resize (numExportIDsUnion);
    for (size_type k = 0; k < numExportIDsUnion; ++k) {
      exportLIDsUnion[k] = srcMap->getLocalElement (exportGIDsUnion[k]);
    }
#endif // TPETRA_IMPORT_SETUNION_USE_CREATE_FROM_SENDS

    // Create and return the union Import. This uses the "expert" constructor
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM2 = Teuchos::null;
    label = "Tpetra::Import::setUnion : Construct Import";
    MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(label)));
#endif
    RCP<const import_type> unionImport =
      rcp (new import_type (srcMap, unionTgtMap,
                            as<size_t> (numSameIDsUnion),
                            permuteToLIDsUnion, permuteFromLIDsUnion,
                            remoteLIDsUnion, exportLIDsUnion,
                            exportPIDsUnion, distributor, this->out_));

    return unionImport;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setUnion () const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > unionImport;
    RCP<const map_type> srcMap = this->getSourceMap ();
    RCP<const map_type> tgtMap = this->getTargetMap ();
    RCP<const Comm<int> > comm = srcMap->getComm ();

    ArrayView<const GO> srcGIDs = srcMap->getNodeElementList ();
    ArrayView<const GO> tgtGIDs = tgtMap->getNodeElementList ();

    // All elements in srcMap will be in the "new" target map, so...
    size_t numSameIDsNew    = srcMap->getNodeNumElements();
    size_t numRemoteIDsNew  = getNumRemoteIDs();
    Array<LO> permuteToLIDsNew, permuteFromLIDsNew; // empty on purpose

    // Grab some old data
    ArrayView<const LO> remoteLIDsOld = getRemoteLIDs();
    ArrayView<const LO> exportLIDsOld = getExportLIDs();

    // Build up the new map (same part)
    Array<GO> GIDs(numSameIDsNew + numRemoteIDsNew);
    for(size_t i=0; i<numSameIDsNew; i++)
      GIDs[i] = srcGIDs[i];

    // Build up the new map (remote part) and remotes list
    Array<LO> remoteLIDsNew(numRemoteIDsNew);
    for(size_t i=0; i<numRemoteIDsNew; i++) {
      GIDs[numSameIDsNew + i] = tgtGIDs[remoteLIDsOld[i]];
      remoteLIDsNew[i] = numSameIDsNew+i;
    }

    // Build the new target map
    GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    RCP<const map_type> targetMapNew = rcp(new map_type(GO_INVALID,GIDs,tgtMap->getIndexBase(),tgtMap->getComm(),tgtMap->getNode()));

    // Exports are trivial (since the sourcemap doesn't change)
    Array<int> exportPIDsnew(getExportPIDs());
    Array<LO> exportLIDsnew(getExportLIDs());

    // Copy the Distributor (due to how the Import constructor works)
    Distributor D(getDistributor());

    // Build the importer using the "expert" constructor
    unionImport = rcp(new Import<LocalOrdinal, GlobalOrdinal, Node>(srcMap,
                                                                    targetMapNew,
                                                                    numSameIDsNew,
                                                                    permuteToLIDsNew,
                                                                    permuteFromLIDsNew,
                                                                    remoteLIDsNew,
                                                                    exportLIDsnew,
                                                                    exportPIDsnew,D));

    return unionImport;
  }




  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  createRemoteOnlyImport(const Teuchos::RCP<const map_type>& remoteTarget) const
  {
    using Teuchos::rcp;
    typedef Import<LocalOrdinal,GlobalOrdinal,Node> import_type;

    const size_t NumRemotes = getNumRemoteIDs ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      NumRemotes != remoteTarget->getNodeNumElements (),
      std::runtime_error, "Tpetra::createRemoteOnlyImport: "
      "remoteTarget map ID count doesn't match.");

    // Compute the new Remote LIDs
    Teuchos::ArrayView<const LocalOrdinal> oldRemoteLIDs = getRemoteLIDs ();
    Teuchos::Array<LocalOrdinal> newRemoteLIDs (NumRemotes);
    for (size_t i = 0; i < NumRemotes; ++i) {
      newRemoteLIDs[i] = remoteTarget->getLocalElement (getTargetMap ()->getGlobalElement (oldRemoteLIDs[i]));
      // Now we make sure these guys are in sorted order (AztecOO-ML ordering)
      TEUCHOS_TEST_FOR_EXCEPTION(
        i > 0 && newRemoteLIDs[i] < newRemoteLIDs[i-1],
        std::runtime_error, "Tpetra::createRemoteOnlyImport: "
        "this and remoteTarget order don't match.");
    }

    // Copy ExportPIDs and such
    // NOTE: Be careful: The "Expert" Import constructor we use does a "swap"
    // for most of the LID/PID lists and the Distributor, meaning it
    // ruins the existing object if we pass things in directly.  Hence
    // we copy them first.
    Teuchos::Array<int> newExportPIDs (getExportPIDs ());
    Teuchos::Array<LocalOrdinal> newExportLIDs (getExportLIDs ());
    Teuchos::Array<LocalOrdinal> dummy;
    Distributor newDistor (getDistributor ());

    return rcp (new import_type (getSourceMap (), remoteTarget,
                                 static_cast<size_t> (0), dummy, dummy,
                                 newRemoteLIDs, newExportLIDs,
                                 newExportPIDs, newDistor));
  }

} // namespace Tpetra

#define TPETRA_IMPORT_CLASS_INSTANT(LO, GO, NODE) \
  \
  template class Import< LO , GO , NODE >;

// Explicit instantiation macro.
// Only invoke this when in the Tpetra namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_IMPORT_INSTANT(LO, GO, NODE) \
  TPETRA_IMPORT_CLASS_INSTANT(LO, GO, NODE)

#endif // TPETRA_IMPORT_DEF_HPP
