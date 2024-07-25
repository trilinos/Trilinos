// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_EXPORT_DEF_HPP
#define TPETRA_EXPORT_DEF_HPP


#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_ImportExportData.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Details_DualViewUtil.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include <memory>

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    base_type (source, target, out, plist, "Export")
  {
    using Teuchos::rcp;
    using std::endl;
    using ::Tpetra::Details::ProfilingRegion;
    ProfilingRegion regionExport ("Tpetra::Export::Export");

    if (this->verbose ()) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor" << endl;
      this->verboseOutputStream () << os.str ();
    }
    Teuchos::Array<GlobalOrdinal> exportGIDs;
    setupSamePermuteExport (exportGIDs);
    if (source->isDistributed ()) {
      setupRemote (exportGIDs);
    }

    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_host () );

    this->detectRemoteExportLIDsContiguous();

    if (this->verbose ()) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << myRank << ": Export ctor: done" << endl;
      this->verboseOutputStream () << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target) :
    Export (source, target, Teuchos::null, Teuchos::null)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out) :
    Export (source, target, out, Teuchos::null)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    Export (source, target, Teuchos::null, plist)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs) :
    base_type (rhs)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Import<LocalOrdinal,GlobalOrdinal,Node>& importer) :
    base_type (importer, typename base_type::reverse_tag ())
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    // Call the base class' method.  It does all the work.
    this->describeImpl (out, "Tpetra::Export", verbLevel);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Export<LocalOrdinal,GlobalOrdinal,Node>::
  print (std::ostream& os) const
  {
    auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));
    // "Print" traditionally meant "everything."
    this->describe (*out, Teuchos::VERB_EXTREME);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  setupSamePermuteExport (Teuchos::Array<GlobalOrdinal>& exportGIDs)
  {
    using ::Tpetra::Details::makeDualViewFromOwningHostView;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::view_alloc_no_init;
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::null;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using size_type = typename ArrayView<const GO>::size_type;
    const char tfecfFuncName[] = "setupSamePermuteExport: ";
    ProfilingRegion regionExport ("Tpetra::Export::setupSamePermuteExport");

    std::unique_ptr<std::string> prefix;
    if (this->verbose ()) {
      auto srcMap = this->getSourceMap ();
      auto comm = srcMap.is_null () ? Teuchos::null : srcMap->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();

      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::Export::setupSamePermuteExport: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));

      std::ostringstream os2;
      os2 << *prefix << "Start" << std::endl;
      this->verboseOutputStream () << os2.str ();
    }

    const map_type& source = * (this->getSourceMap ());
    const map_type& target = * (this->getTargetMap ());
    ArrayView<const GO> sourceGIDs = source.getLocalElementList ();
    ArrayView<const GO> targetGIDs = target.getLocalElementList ();

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
    for ( ; numSameGids < numGids &&
            rawSrcGids[numSameGids] == rawTgtGids[numSameGids];
          ++numSameGids)
      {} // third clause of 'for' does everything
    this->TransferData_->numSameIDs_ = numSameGids;

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "numIDs: " << numGids
         << ", numSameIDs: " << numSameGids << endl;
      this->verboseOutputStream () << os.str ();
    }

    // Compute permuteToLIDs_, permuteFromLIDs_, exportGIDs, and
    // exportLIDs_.  The first two arrays are IDs to be permuted, and
    // the latter two arrays are IDs to sent out ("exported"), called
    // "export" IDs.
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  IDs to send are in the
    // source Map, but not in the target Map.

    // Iterate over the source Map's LIDs, since we only need to do
    // GID -> LID lookups for the target Map.
    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO numSrcLids = static_cast<LO> (numSrcGids);
    LO numPermutes = 0;
    LO numExports = 0;

    for (LO srcLid = numSameGids; srcLid < numSrcLids; ++srcLid) {
      const GO curSrcGid = rawSrcGids[srcLid];
      // getLocalElement() returns LINVALID if the GID isn't in the
      // target Map.  This saves us a lookup (which
      // isNodeGlobalElement() would do).
      const LO tgtLid = target.getLocalElement (curSrcGid);
      if (tgtLid != LINVALID) { // if target.isNodeGlobalElement (curSrcGid)
        ++numPermutes;
      }
      else {
        ++numExports;
      }
    }
    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "numPermutes: " << numPermutes
         << ", numExports: " << numExports << endl;
      this->verboseOutputStream () << os.str ();
    }
    TEUCHOS_ASSERT( numPermutes + numExports ==
                    numSrcLids - numSameGids );

    typename decltype (this->TransferData_->permuteToLIDs_)::t_host
      permuteToLIDs (view_alloc_no_init ("permuteToLIDs"), numPermutes);
    typename decltype (this->TransferData_->permuteToLIDs_)::t_host
      permuteFromLIDs (view_alloc_no_init ("permuteFromLIDs"), numPermutes);
    typename decltype (this->TransferData_->permuteToLIDs_)::t_host
      exportLIDs (view_alloc_no_init ("exportLIDs"), numExports);

    // FIXME (mfh 03 Feb 2019) Replace with std::unique_ptr of array,
    // to avoid superfluous initialization on resize.
    exportGIDs.resize (numExports);

    {
      LO numPermutes2 = 0;
      LO numExports2 = 0;
      for (LO srcLid = numSameGids; srcLid < numSrcLids; ++srcLid) {
        const GO curSrcGid = rawSrcGids[srcLid];
        const LO tgtLid = target.getLocalElement (curSrcGid);
        if (tgtLid != LINVALID) {
          permuteToLIDs[numPermutes2] = tgtLid;
          permuteFromLIDs[numPermutes2] = srcLid;
          ++numPermutes2;
        }
        else {
          exportGIDs[numExports2] = curSrcGid;
          exportLIDs[numExports2] = srcLid;
          ++numExports2;
        }
      }
      TEUCHOS_ASSERT( numPermutes == numPermutes2 );
      TEUCHOS_ASSERT( numExports == numExports2 );
      TEUCHOS_ASSERT( size_t (numExports) == size_t (exportGIDs.size ()) );
    }

    // Defer making this->TransferData_->exportLIDs_ until after
    // getRemoteIndexList, since we might need to shrink it then.

    // exportLIDs is the list of this process' LIDs that it has to
    // send out.  Since this is an Export, and therefore the target
    // Map is nonoverlapping, we know that each export LID only needs
    // to be sent to one process.  However, the source Map may be
    // overlapping, so multiple processes might send to the same LID
    // on a receiving process.

    if (numExports != 0 && ! source.isDistributed ()) {
      // This Export has export LIDs, meaning that the source Map has
      // entries on this process that are not in the target Map on
      // this process.  However, the source Map is not distributed
      // globally.  This implies that this Import is not locally
      // complete on this process.
      this->TransferData_->isLocallyComplete_ = false;
      if (this->verbose ()) {
        std::ostringstream os;
        os << *prefix << "Export is not locally complete" << endl;
        this->verboseOutputStream () << os.str ();
      }
      // mfh 12 Sep 2016: I disagree that this is "abuse"; it may be
      // correct behavior, depending on the circumstances.
      TPETRA_ABUSE_WARNING
        (true, std::runtime_error, "::setupSamePermuteExport(): Source has "
         "export LIDs but Source is not distributed globally.  Exporting to "
         "a submap of the target map.");
    }

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
      if (this->verbose ()) {
        std::ostringstream os;
        os << *prefix << "Source Map is distributed; "
          "call targetMap.getRemoteiNdexList" << endl;
        this->verboseOutputStream () << os.str ();
      }
      this->TransferData_->exportPIDs_.resize(exportGIDs.size ());
      // This call will assign any GID in the target Map with no
      // corresponding process ID a fake process ID of -1.  We'll use
      // this below to remove exports for processses that don't exist.
      const LookupStatus lookup =
        target.getRemoteIndexList (exportGIDs(),
                                   this->TransferData_->exportPIDs_ ());
      // mfh 12 Sep 2016: I disagree that this is "abuse"; it may be
      // correct behavior, depending on the circumstances.
      TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
        "::setupSamePermuteExport(): The source Map has GIDs not found "
        "in the target Map.");

      // Get rid of process IDs not in the target Map.  This prevents
      // exporting to GIDs which don't belong to any process in the
      // target Map.
      if (lookup == IDNotPresent) {
        // There is at least one GID owned by the calling process in
        // the source Map, which is not owned by any process in the
        // target Map.
        this->TransferData_->isLocallyComplete_ = false;

        Teuchos::Array<int>& exportPIDs = this->TransferData_->exportPIDs_;

        const size_type totalNumExports = exportPIDs.size ();
        const size_type numInvalidExports =
          std::count_if (exportPIDs.begin (), exportPIDs.end (),
                         [] (const int procId) { return procId == -1; });
        if (this->verbose ()) {
          std::ostringstream os;
          os << *prefix << "totalNumExports: " << totalNumExports
             << ", numInvalidExports: " << numInvalidExports << endl;
          this->verboseOutputStream () << os.str ();
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numInvalidExports == 0, std::logic_error,
           "targetMap.getRemoteIndexList returned IDNotPresent, but no export "
           "PIDs are -1.  Please report this bug to the Tpetra developers.");

        // We know that at least one export ID is invalid, that is,
        // not in any process on the target Map.  If all export IDs
        // are invalid, we can delete all exports.  Otherwise, keep
        // the valid exports and discard the rest.  This is legit
        // Petra Object Model behavior, but it's a less common case.

        if (numInvalidExports == totalNumExports) {
          exportGIDs.resize (0);
          exportLIDs = decltype (exportLIDs) ();
          exportPIDs.resize (0);
        }
        else {
          size_type numValidExports = 0;
          for (size_type e = 0; e < totalNumExports; ++e) {
            if (this->TransferData_->exportPIDs_[e] != -1) {
              exportGIDs[numValidExports] = exportGIDs[e];
              exportLIDs[numValidExports] = exportLIDs[e];
              exportPIDs[numValidExports] = exportPIDs[e];
              ++numValidExports;
            }
          }
          exportGIDs.resize (numValidExports);
          Kokkos::resize (exportLIDs, numValidExports);
          exportPIDs.resize (numValidExports);
        }
      }
    }

    // FIXME (mfh 03 Feb 2019) These three DualViews could share a
    // single device allocation, in order to avoid high cudaMalloc
    // cost and device memory fragmentation.
    makeDualViewFromOwningHostView (this->TransferData_->permuteToLIDs_, permuteToLIDs);
    makeDualViewFromOwningHostView (this->TransferData_->permuteFromLIDs_, permuteFromLIDs);
    makeDualViewFromOwningHostView (this->TransferData_->exportLIDs_, exportLIDs);

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Done!" << std::endl;
      this->verboseOutputStream () << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  setupRemote (Teuchos::Array<GlobalOrdinal>& exportGIDs)
  {
    using ::Tpetra::Details::view_alloc_no_init;
    using ::Tpetra::Details::makeDualViewFromOwningHostView;
    using Teuchos::Array;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;

    std::unique_ptr<std::string> prefix;
    if (this->verbose ()) {
      auto srcMap = this->getSourceMap ();
      auto comm = srcMap.is_null () ? Teuchos::null : srcMap->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();

      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::Export::setupRemote: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));

      std::ostringstream os2;
      os2 << *prefix << "Start" << std::endl;
      this->verboseOutputStream () << os2.str ();
    }

    TEUCHOS_ASSERT( ! this->getTargetMap ().is_null () );
    const map_type& tgtMap = * (this->getTargetMap ());

    // Sort exportPIDs_ in ascending order, and apply the same
    // permutation to exportGIDs_ and exportLIDs_.  This ensures that
    // exportPIDs_[i], exportGIDs_[i], and exportLIDs_[i] all
    // refer to the same thing.
    {
      TEUCHOS_ASSERT( size_t (this->TransferData_->exportLIDs_.extent (0)) ==
                      size_t (this->TransferData_->exportPIDs_.size ()) );
      this->TransferData_->exportLIDs_.modify_host ();
      auto exportLIDs = this->TransferData_->exportLIDs_.view_host ();
      sort3 (this->TransferData_->exportPIDs_.begin (),
             this->TransferData_->exportPIDs_.end (),
             exportGIDs.getRawPtr (),
             exportLIDs.data ());
      this->TransferData_->exportLIDs_.sync_device ();
      // FIXME (mfh 03 Feb 2019) We actually end up sync'ing
      // exportLIDs_ to device twice, once in setupSamePermuteExport,
      // and once here.  We could avoid the first sync.
    }

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Call createFromSends" << endl;
      this->verboseOutputStream () << os.str ();
    }

    // Construct the list of entries that calling image needs to send
    // as a result of everyone asking for what it needs to receive.
    //
    // mfh 05 Jan 2012: I understand the above comment as follows:
    // Construct the communication plan from the list of image IDs to
    // which we need to send.
    Teuchos::Array<int>& exportPIDs = this->TransferData_->exportPIDs_;
    Distributor& distributor = this->TransferData_->distributor_;
    const size_t numRemoteIDs = distributor.createFromSends (exportPIDs ());

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "numRemoteIDs: " << numRemoteIDs
         << "; call doPostsAndWaits" << endl;
      this->verboseOutputStream () << os.str ();
    }

    // Use the communication plan with ExportGIDs to find out who is
    // sending to us and get the proper ordering of GIDs for incoming
    // remote entries (these will be converted to LIDs when done).

    Kokkos::View<const GO*, Kokkos::HostSpace> exportGIDsConst(exportGIDs.data(), exportGIDs.size());
    Kokkos::View<GO*, Kokkos::HostSpace> remoteGIDs("remoteGIDs", numRemoteIDs);
    distributor.doPostsAndWaits(exportGIDsConst, 1, remoteGIDs);

    // Remote (incoming) IDs come in as GIDs; convert to LIDs.  LIDs
    // tell this process where to store the incoming remote data.
    using host_remote_lids_type =
      typename decltype (this->TransferData_->remoteLIDs_)::t_host;
    host_remote_lids_type remoteLIDs
      (view_alloc_no_init ("remoteLIDs"), numRemoteIDs);

    for (LO j = 0; j < LO (numRemoteIDs); ++j) {
      remoteLIDs[j] = tgtMap.getLocalElement (remoteGIDs[j]);
    }
    makeDualViewFromOwningHostView (this->TransferData_->remoteLIDs_, remoteLIDs);

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl;
      this->verboseOutputStream () << os.str ();
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
  template class Export< LO , GO , NODE >;

#endif // TPETRA_EXPORT_DEF_HPP
