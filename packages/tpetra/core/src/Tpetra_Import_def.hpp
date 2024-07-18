// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_IMPORT_DEF_HPP
#define TPETRA_IMPORT_DEF_HPP

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_ImportExportData.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_DualViewUtil.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Teuchos_as.hpp"
#ifdef HAVE_TPETRA_MMM_TIMINGS
#include "Teuchos_TimeMonitor.hpp"
#endif
#include <array>
#include <memory>

namespace Teuchos {
  template<class T>
  std::string toString (const std::vector<T>& x)
  {
    std::ostringstream os;
    os << "[";
    const std::size_t N = x.size ();
    for (std::size_t k = 0; k < N; ++k) {
      os << x[k];
      if (k + std::size_t (1) < N) {
        os << ",";
      }
    }
    os << "]";
    return os.str ();
  }

  template<class ElementType, class DeviceType>
  std::string toString (const Kokkos::View<const ElementType*, DeviceType>& x)
  {
    std::ostringstream os;
    os << "[";
    const std::size_t N = std::size_t (x.extent (0));
    for (std::size_t k = 0; k < N; ++k) {
      os << x[k];
      if (k + std::size_t (1) < N) {
        os << ",";
      }
    }
    os << "]";
    return os.str ();
  }
} // namespace Teuchos

namespace Tpetra {

  // head:  init(source, target, true, remotePIDs, Teuchos::null);

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  init (const Teuchos::RCP<const map_type>& source,
        const Teuchos::RCP<const map_type>& /* target */,
        bool useRemotePIDs,
        Teuchos::Array<int> & remotePIDs,
        const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::Ptr;
    using Teuchos::rcp;
    using std::endl;
    ProfilingRegion regionImportInit ("Tpetra::Import::init");

    std::unique_ptr<std::string> verbPrefix;
    if (this->verbose ()) {
      std::ostringstream os;
      const int myRank = source->getComm ()->getRank ();
      os << "Proc " << myRank << ": Tpetra::Import::init: ";
      verbPrefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << endl;
      this->verboseOutputStream () << os.str ();
    }

    Array<GlobalOrdinal> remoteGIDs;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    using Teuchos::TimeMonitor;
    std::string label;
    if(!plist.is_null())
      label = plist->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(":iport_ctor:preIData: ");
#else
    (void)plist;
#endif
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      auto MM(*TimeMonitor::getNewTimer(prefix));
#endif
      if (this->verbose ()) {
        std::ostringstream os;
        os << *verbPrefix << "Call setupSamePermuteRemote" << endl;
        this->verboseOutputStream () << os.str ();
      }
      setupSamePermuteRemote (remoteGIDs);
    }
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      prefix = std::string("Tpetra ")+ label + std::string(":iport_ctor:preSetupExport: ");
      auto MM2(*TimeMonitor::getNewTimer(prefix));
#endif
      if (source->isDistributed ()) {
        if (this->verbose ()) {
          std::ostringstream os;
          os << *verbPrefix << "Call setupExport" << endl;
          this->verboseOutputStream () << os.str ();
        }
        setupExport (remoteGIDs,useRemotePIDs,remotePIDs);
      }
      else if (this->verbose ()) {
        std::ostringstream os;
        os << *verbPrefix << "Source Map not distributed; skip setupExport"
           << endl;
        this->verboseOutputStream () << os.str ();
      }
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
      os << *verbPrefix << "Done!" << endl;
      this->verboseOutputStream () << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type >& source,
          const Teuchos::RCP<const map_type >& target) :
    base_type (source, target, Teuchos::null, Teuchos::null, "Import")
  {
    Teuchos::Array<int> dummy;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    Teuchos::RCP<Teuchos::ParameterList> mypars = rcp(new Teuchos::ParameterList);
    mypars->set("Timer Label","Naive_tAFC");
    init(source, target, false, dummy, mypars);
#else
    init (source, target, false, dummy, Teuchos::null);
#endif
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out) :
    base_type (source, target, out, Teuchos::null, "Import")
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, Teuchos::null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    base_type (source, target, Teuchos::null, plist, "Import")
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    base_type (source, target, out, plist, "Import")
  {
    Teuchos::Array<int> dummy;
    init (source, target, false, dummy, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const map_type>& source,
          const Teuchos::RCP<const map_type>& target,
          Teuchos::Array<int>& remotePIDs,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    base_type (source, target, Teuchos::null, plist, "Import")
  {
    init (source, target, true, remotePIDs, plist);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Import<LocalOrdinal,GlobalOrdinal,Node>& rhs) :
    base_type (rhs)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter) :
    base_type (exporter, typename base_type::reverse_tag ())
  {}

  // cblcbl
  // This is the "createExpert" version of the constructor to be used with pid/gid pairs obtained from
  // reverse communication
   template <class LocalOrdinal, class GlobalOrdinal, class Node>
   Import<LocalOrdinal,GlobalOrdinal,Node>::
   Import (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& source,
           const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& target,
           const Teuchos::ArrayView<int> & userRemotePIDs,
           const Teuchos::ArrayView<const LocalOrdinal> & userExportLIDs,
           const Teuchos::ArrayView<const int>          & userExportPIDs,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<Teuchos::FancyOStream>& out) :
     base_type (source, target, out, plist, "Import")
  {
    using ::Tpetra::Details::makeDualViewFromArrayView;
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::rcp;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using size_type = Teuchos::Array<int>::size_type;

    std::unique_ptr<std::string> prefix;
    if (this->verbose ()) {
      auto comm = source.is_null () ? Teuchos::null : source->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::Import createExpert ctor: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << endl;
      this->verboseOutputStream () << os.str ();
    }

    ArrayView<const GO> sourceGIDs = source->getLocalElementList ();
    ArrayView<const GO> targetGIDs = target->getLocalElementList ();

    Array<GO> tRemoteGIDs;
    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Call setupSamePermuteRemote" << endl;
      this->verboseOutputStream () << os.str ();
    }
    setupSamePermuteRemote (tRemoteGIDs);

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Sort & filter IDs" << endl;
      this->verboseOutputStream () << os.str ();
    }

    auto tRemoteLIDs = this->TransferData_->remoteLIDs_.view_host ();
    this->TransferData_->remoteLIDs_.modify_host ();
    Teuchos::Array<int> tRemotePIDs (userRemotePIDs);

    if (this->verbose () && this->getNumRemoteIDs () > 0 && ! source->isDistributed ()) {
      std::ostringstream os;
      os << *prefix << "Target Map has remote LIDs but source Map is not "
        "distributed.  Importing to a submap of the target Map." << endl;
      this->verboseOutputStream () << os.str ();
    }
    // FIXME (mfh 03 Feb 2019) I don't see this as "abuse"; it's
    // perfectly valid Petra Object Model behavior.
    TPETRA_ABUSE_WARNING
      (getNumRemoteIDs () > 0 && ! source->isDistributed (),
       std::runtime_error,
       "::constructExpert: Target Map has remote LIDs but source Map "
       "is not distributed.  Importing to a submap of the target Map.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (tRemotePIDs.size () != tRemoteGIDs.size () ||
       size_t (tRemoteGIDs.size ()) != size_t (tRemoteLIDs.extent (0)),
       std::runtime_error, "Import::Import createExpert version: "
       "Size mismatch on userRemotePIDs, remoteGIDs, and remoteLIDs "
       "Array's to sort3.");

    sort3 (tRemotePIDs.begin (),
           tRemotePIDs.end (),
           tRemoteGIDs.begin (),
           tRemoteLIDs.data ());

    //Get rid of IDs that don't exist in SourceMap
    size_type cnt = 0;
    size_type indexIntoRemotePIDs = tRemotePIDs.size ();
    for (size_type i = 0; i < indexIntoRemotePIDs; ++i) {
      if (tRemotePIDs[i] == -1) {
        ++cnt;
      }
    }

    if (cnt == 0) { // done modifying remoteLIDs_
      this->TransferData_->remoteLIDs_.sync_device ();
    }
    else {
      if (indexIntoRemotePIDs - cnt > 0) {
        Array<GO>  newRemoteGIDs(indexIntoRemotePIDs-cnt);
        Array<LO>  newRemoteLIDs(indexIntoRemotePIDs-cnt);
        Array<int> newRemotePIDs(indexIntoRemotePIDs-cnt);
        cnt = 0;

        for (size_type j = 0; j < indexIntoRemotePIDs; ++j)
          if(tRemotePIDs[j] != -1) {
            newRemoteGIDs[cnt] = tRemoteGIDs[j];
            newRemotePIDs[cnt] = tRemotePIDs[j];
            newRemoteLIDs[cnt] = target->getLocalElement(tRemoteGIDs[j]);
            ++cnt;
          }
        indexIntoRemotePIDs = cnt;
        tRemoteGIDs = newRemoteGIDs;
        tRemotePIDs = newRemotePIDs;
        makeDualViewFromArrayView (this->TransferData_->remoteLIDs_,
                                   newRemoteLIDs ().getConst (),
                                   "remoteLIDs");
      }
      else { //valid RemoteIDs empty
        indexIntoRemotePIDs = 0;
        tRemoteGIDs.clear();
        tRemotePIDs.clear();
        this->TransferData_->remoteLIDs_ = decltype (this->TransferData_->remoteLIDs_) ();
      }
    }

    this->TransferData_->exportPIDs_ = Teuchos::Array<int> (userExportPIDs);
    makeDualViewFromArrayView (this->TransferData_->exportLIDs_,
                               userExportLIDs, "exportLIDs");

    bool locallyComplete = true;
    for (size_type i = 0; i < userExportPIDs.size () && locallyComplete; ++i) {
      if (userExportPIDs[i] == -1) {
        locallyComplete = false;
      }
    }
    this->TransferData_->isLocallyComplete_ = locallyComplete;

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "locallyComplete: "
         << (locallyComplete ? "true" : "false")
         << "; call createFromSendsAndRecvs" << endl;
      this->verboseOutputStream () << os.str ();
    }
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string mmm_prefix =
        std::string("Tpetra ") + std::string(":iport_ctor:cFSAR ");

      auto MM3(*Teuchos::TimeMonitor::getNewTimer(mmm_prefix));
#endif
      Distributor& distributor = this->TransferData_->distributor_;
      distributor.createFromSendsAndRecvs (this->TransferData_->exportPIDs_, tRemotePIDs);
    }

    this->detectRemoteExportLIDsContiguous();

    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_host () );
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
    base_type (source, target, out, plist, "Import")
  {
    using ::Tpetra::Details::makeDualViewFromArrayView;
    using std::endl;

    std::unique_ptr<std::string> prefix;
    if (this->verbose ()) {
      auto comm = source.is_null () ? Teuchos::null : source->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::Import export ctor: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << endl;
      this->verboseOutputStream () << os.str ();
    }

    bool locallyComplete = true;
    for (Teuchos::Array<int>::size_type i = 0; i < exportPIDs.size (); ++i) {
      if (exportPIDs[i] == -1) {
        locallyComplete = false;
      }
    }
    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "numSameIDs: " << numSameIDs << ", locallyComplete: "
         << (locallyComplete ? "true" : "false") << endl;
      this->verboseOutputStream () << os.str ();
    }

    this->TransferData_->isLocallyComplete_ = locallyComplete;
    this->TransferData_->numSameIDs_ = numSameIDs;

    makeDualViewFromArrayView (this->TransferData_->permuteToLIDs_,
                               permuteToLIDs ().getConst (),
                               "permuteToLIDs");
    TEUCHOS_ASSERT( size_t (this->TransferData_->permuteToLIDs_.extent (0)) ==
                    size_t (permuteToLIDs.size ()) );
    makeDualViewFromArrayView (this->TransferData_->permuteFromLIDs_,
                               permuteFromLIDs ().getConst (),
                               "permuteFromLIDs");
    TEUCHOS_ASSERT( size_t (this->TransferData_->permuteFromLIDs_.extent (0)) ==
                    size_t (permuteFromLIDs.size ()) );
    makeDualViewFromArrayView (this->TransferData_->remoteLIDs_,
                               remoteLIDs ().getConst (),
                               "remoteLIDs");
    TEUCHOS_ASSERT( size_t (this->TransferData_->remoteLIDs_.extent (0)) ==
                    size_t (remoteLIDs.size ()) );
    makeDualViewFromArrayView (this->TransferData_->exportLIDs_,
                               exportLIDs ().getConst (),
                               "exportLIDs");
    TEUCHOS_ASSERT( size_t (this->TransferData_->exportLIDs_.extent (0)) ==
                    size_t (exportLIDs.size ()) );
    this->TransferData_->exportPIDs_.swap (exportPIDs);
    this->TransferData_->distributor_.swap (distributor);

    this->detectRemoteExportLIDsContiguous();

    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteFromLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->permuteToLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->remoteLIDs_.need_sync_host () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_device () );
    TEUCHOS_ASSERT( ! this->TransferData_->exportLIDs_.need_sync_host () );
  }

  namespace { // (anonymous)

    template <class LO, class GO, class NT>
    struct ImportLocalSetupResult
    {
      Teuchos::RCP<const ::Tpetra::Map<LO, GO, NT> > targetMap;
      LO numSameIDs;
      // std::vector<LO> permuteToLIDs; // users aren't supposed to have permutes
      // std::vector<LO> permuteFromLIDs; // users aren't suppoosed to have permutes
      std::vector<GO> remoteGIDs;
      std::vector<LO> remoteLIDs;
      std::vector<int> remotePIDs;
      LO numPermutes; // users aren't supposed to have permutes
    };

    template<class T>
    void printArray (std::ostream& out, const T x[], const std::size_t N)
    {
      out << "[";
      for (std::size_t k = 0; k < N; ++k) {
        out << x[k];
        if (k + 1 < N) {
          out << ", ";
        }
      }
      out << "]";
    }

    template<class LO, class GO, class NT>
    ImportLocalSetupResult<LO, GO, NT>
    setupSamePermuteRemoteFromUserGlobalIndexList (const ::Tpetra::Map<LO, GO, NT>& sourceMap,
                                                   const GO targetMapRemoteOrPermuteGlobalIndices[],
                                                   const int targetMapRemoteOrPermuteProcessRanks[],
                                                   const LO numTargetMapRemoteOrPermuteGlobalIndices,
                                                   const bool mayReorderTargetMapIndicesLocally,
                                                   Teuchos::FancyOStream* out, // only valid if verbose
                                                   const std::string* verboseHeader, // only valid if verbose
                                                   const bool verbose,
                                                   const bool debug)
    {
      using std::endl;
      const int myRank = sourceMap.getComm ()->getRank ();
      ImportLocalSetupResult<LO, GO, NT> result;

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Import::setupSPR w/ remote GIDs & PIDs: " << endl
            << *verboseHeader << "  Input GIDs: ";
        printArray (os, targetMapRemoteOrPermuteGlobalIndices, numTargetMapRemoteOrPermuteGlobalIndices);
        os << endl << "  Input PIDs: ";
        printArray (os, targetMapRemoteOrPermuteProcessRanks, numTargetMapRemoteOrPermuteGlobalIndices);
        os << endl;
        *out << os.str ();
      }

      // In debug mode, check whether any of the input GIDs are
      // actually in the source Map on the calling process.  That's an
      // error, because it means duplicate GIDs on the calling
      // process.  Also check if any of the input PIDs are invalid.
      if (debug) {
        std::vector<GO> badGIDs;
        std::vector<int> badPIDs;
        const Teuchos::Comm<int>& comm = * (sourceMap.getComm ());
        const int numProcs = comm.getSize ();

        for (LO k = 0; k < numTargetMapRemoteOrPermuteGlobalIndices; ++k) {
          const GO tgtGID = targetMapRemoteOrPermuteGlobalIndices[k];
          if (sourceMap.isNodeGlobalElement (tgtGID)) {
            badGIDs.push_back (tgtGID);
          }
          const int tgtPID = targetMapRemoteOrPermuteProcessRanks[k];
          if (tgtPID < 0 || tgtPID >= numProcs) {
            badPIDs.push_back (tgtPID);
          }
        }

        std::array<int, 2> lclStatus {{
          badGIDs.size () == 0 ? 1 : 0,
          badPIDs.size () == 0 ? 1 : 0
        }};
        std::array<int, 2> gblStatus {{0, 0}}; // output argument
        Teuchos::reduceAll<int, int> (comm, Teuchos::REDUCE_MIN, 2,
                                      lclStatus.data (), gblStatus.data ());
        const bool good = gblStatus[0] == 1 && gblStatus[1] == 1;
        // Don't actually print all the "bad" GIDs and/or PIDs unless
        // in verbose mode, since there could be many of them.
        if (verbose && gblStatus[0] != 1) {
          std::ostringstream os;
          os << *verboseHeader << "- Some input GIDs are already in the source Map: ";
          printArray (os, badGIDs.data (), badGIDs.size ());
          os << endl;
          ::Tpetra::Details::gathervPrint (*out, os.str (), comm);
        }
        if (verbose && gblStatus[0] != 1) {
          std::ostringstream os;
          os << *verboseHeader << "- Some input PIDs are invalid: ";
          printArray (os, badPIDs.data (), badPIDs.size ());
          os << endl;
          ::Tpetra::Details::gathervPrint (*out, os.str (), comm);
        }

        if (! good) {
          std::ostringstream os;
          os << "Tpetra::Import constructor that takes remote GIDs and PIDs: ";
          if (gblStatus[0] != 1) {
            os << "Some input GIDs (global indices) are already in the source Map!  ";
          }
          if (gblStatus[1] != 1) {
            os << "Some input PIDs (process ranks) are invalid!  ";
          }
          os << "Rerun with the environment variable TPETRA_VERBOSE=Tpetra::Import "
            "to see what GIDs and/or PIDs are bad.";
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
        }
      }

      // Create list of GIDs to go into target Map.  We need to copy
      // the GIDs into this list anyway, so once we have them, we can
      // sort the "remotes" in place.
      const LO numLclSrcIDs = static_cast<LO> (sourceMap.getLocalNumElements ());
      const LO numLclTgtIDs = numLclSrcIDs + numTargetMapRemoteOrPermuteGlobalIndices;
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Copy source Map GIDs into target Map GID list: "
          "numLclSrcIDs=" << numLclSrcIDs
           << ", numTargetMapRemoteOrPermuteGlobalIndices="
           << numTargetMapRemoteOrPermuteGlobalIndices << endl;
        *out << os.str ();
      }
      std::vector<GO> tgtGIDs (numLclTgtIDs); // will go into target Map ctor
      if (sourceMap.isContiguous ()) {
        GO curTgtGID = sourceMap.getMinGlobalIndex ();
        for (LO k = 0; k < numLclSrcIDs; ++k, ++curTgtGID) {
          tgtGIDs[k] = curTgtGID;
        }
      }
      else { // avoid calling getLocalElementList on a contiguous Map
        auto srcGIDs = sourceMap.getLocalElementList (); // Teuchos::ArrayView has a different
        for (LO k = 0; k < numLclSrcIDs; ++k) {         // iterator type, so can't std::copy
          tgtGIDs[k] = srcGIDs[k];
        }
      }
      std::copy (targetMapRemoteOrPermuteGlobalIndices,
                 targetMapRemoteOrPermuteGlobalIndices + numTargetMapRemoteOrPermuteGlobalIndices,
                 tgtGIDs.begin () + numLclSrcIDs);

      // Optionally, sort input by process rank, so that remotes
      // coming from the same process are grouped together.  Only sort
      // remote GIDs.  While doing so, detect permutes (input "remote"
      // GIDs whose rank is the same as that of the calling process).
      //
      // Permutes are actually an error.  We normally detect them in
      // debug mode, but if we sort, we have a nearly free opportunity
      // to do so.  We may also safely ignore permutes as duplicates.
      //
      // NOTE: tgtPIDs only includes remotes, not source Map entries.
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Sort by PID? "
           << (mayReorderTargetMapIndicesLocally ? "true" : "false")  << endl;
        *out << os.str ();
      }
      std::vector<int> tgtPIDs (targetMapRemoteOrPermuteProcessRanks,
                                targetMapRemoteOrPermuteProcessRanks + numTargetMapRemoteOrPermuteGlobalIndices);
      result.numPermutes = 0;
      if (mayReorderTargetMapIndicesLocally) {
        Tpetra::sort2 (tgtPIDs.begin (), tgtPIDs.end (), tgtGIDs.begin () + numLclSrcIDs);
        auto range = std::equal_range (tgtPIDs.begin (), tgtPIDs.end (), myRank); // binary search
        if (range.second > range.first) {
          result.numPermutes = static_cast<LO> (range.second - range.first);
        }
      }
      else { // don't sort; linear search to count permutes
        result.numPermutes = static_cast<LO> (std::count (tgtPIDs.begin (), tgtPIDs.end (), myRank));
      }
      // The _actual_ number of remotes.
      const LO numRemotes = numTargetMapRemoteOrPermuteGlobalIndices - result.numPermutes;
      result.numSameIDs = static_cast<LO> (sourceMap.getLocalNumElements ());

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- numSame=" << result.numSameIDs
           << ", numPermutes=" << result.numPermutes
           << ", numRemotes=" << numRemotes << endl;
        *out << os.str ();
      }

      if (result.numPermutes == 0) {
        if (verbose) {
          std::ostringstream os;
          os << *verboseHeader << "- No permutes" << endl;
          *out << os.str ();
        }
        result.remoteGIDs = std::vector<GO> (tgtGIDs.begin () + numLclSrcIDs, tgtGIDs.end ());
        result.remotePIDs.swap (tgtPIDs);
        result.remoteLIDs.resize (numRemotes);
        for (LO k = 0; k < numRemotes; ++k) {
          const LO tgtLid = result.numSameIDs + k;
          result.remoteLIDs[k] = tgtLid;
        }
        if (verbose) {
          std::ostringstream os;
          os << *verboseHeader << "- Remote GIDs: "
             << Teuchos::toString (result.remoteGIDs) << endl;
          os << *verboseHeader << "- Remote PIDs: "
             << Teuchos::toString (result.remotePIDs) << endl;
          os << *verboseHeader << "- Remote LIDs: "
             << Teuchos::toString (result.remoteLIDs) << endl;
          *out << os.str ();
        }
      }
      else { // separate permutes from remotes
        // This case doesn't need to be optimal; it just needs to be
        // correct.  Users really shouldn't give permutes to this
        // Import constructor.
        result.remoteGIDs.reserve (numRemotes);
        result.remoteLIDs.reserve (numRemotes);
        result.remotePIDs.reserve (numRemotes);
        for (LO k = 0; k < numTargetMapRemoteOrPermuteGlobalIndices; ++k) {
          const LO tgtLid = result.numSameIDs + k;
          const GO tgtGid = tgtGIDs[numLclSrcIDs + k];
          const int tgtPid = tgtPIDs[k];

          if (tgtPid != myRank) { // it's a remote
            result.remoteGIDs.push_back (tgtGid);
            result.remoteLIDs.push_back (tgtLid);
            result.remotePIDs.push_back (tgtPid);
          }
        }
        if (verbose) {
          std::ostringstream os;
          os << *verboseHeader << "- Some permutes" << endl;
          *out << os.str ();
        }
      }

      if (sourceMap.isDistributed ()) {
        if (verbose) {
          std::ostringstream os;
          os << *verboseHeader << "- Sort remotes by PID, as Import always does"
             << endl
             << *verboseHeader << "-- remotePIDs before: "
             << Teuchos::toString (result.remotePIDs) << endl
             << *verboseHeader << "-- remoteGIDs before: "
             << Teuchos::toString (result.remoteGIDs) << endl
             << *verboseHeader << "-- remoteLIDs before: "
             << Teuchos::toString (result.remoteLIDs) << endl;
          *out << os.str ();
        }
        // Import always sorts these, regardless of what the user wanted.
        sort3 (result.remotePIDs.begin (),
               result.remotePIDs.end (),
               result.remoteGIDs.begin (),
               result.remoteLIDs.begin ());
        if (verbose) {
          std::ostringstream os;
          os << *verboseHeader << "-- remotePIDs after: "
             << Teuchos::toString (result.remotePIDs) << endl
             << *verboseHeader << "-- remoteGIDs after: "
             << Teuchos::toString (result.remoteGIDs) << endl
             << *verboseHeader << "-- remoteLIDs after: "
             << Teuchos::toString (result.remoteLIDs) << endl;
          std::cerr << os.str ();
        }
      }

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Make target Map" << endl;
        *out << os.str ();
      }
      using ::Teuchos::rcp;
      typedef ::Tpetra::Map<LO, GO, NT> map_type;
      typedef ::Tpetra::global_size_t GST;
      const GST MAP_COMPUTES_GLOBAL_COUNT = ::Teuchos::OrdinalTraits<GST>::invalid ();
      result.targetMap = rcp (new map_type (MAP_COMPUTES_GLOBAL_COUNT,
                                            tgtGIDs.data (),
                                            numLclTgtIDs,
                                            sourceMap.getIndexBase (),
                                            sourceMap.getComm ()));
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Done with sameSPR..." << endl;
        *out << os.str ();
      }
      return result;
    }
  } // namespace (anonymous)

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& sourceMap,
          const GlobalOrdinal targetMapRemoteOrPermuteGlobalIndices[],
          const int targetMapRemoteOrPermuteProcessRanks[],
          const LocalOrdinal numTargetMapRemoteOrPermuteGlobalIndices,
          const bool mayReorderTargetMapIndicesLocally,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<Teuchos::FancyOStream>& debugOutput) :
    // Special case: target Map is null on base_type construction.
    // It's worthwhile for invariants like out_ not being null.
    // We'll set TransferData_ again below.
    base_type (sourceMap, Teuchos::null, debugOutput, plist, "Import")
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::makeDualViewFromOwningHostView;
    using ::Tpetra::Details::makeDualViewFromVector;
    using ::Tpetra::Details::printDualView;
    using ::Tpetra::Details::view_alloc_no_init;
    using Teuchos::ArrayView;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Node NT;

    const bool debug = Behavior::debug ("Import") ||
      Behavior::debug ("Tpetra::Import");

    std::unique_ptr<std::string> verbPfx;
    if (this->verbose ()) {
      std::ostringstream os;
      const int myRank = sourceMap->getComm ()->getRank ();
      os << "Proc " << myRank << ": Tpetra::Import ctor from remotes: ";
      verbPfx = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "mayReorder=" << (mayReorderTargetMapIndicesLocally ? "true" : "false")
         << endl;
      this->verboseOutputStream () << os.str ();
    }

    TEUCHOS_ASSERT( ! this->TransferData_.is_null () );
    ImportLocalSetupResult<LO, GO, NT> localSetupResult =
      setupSamePermuteRemoteFromUserGlobalIndexList<LO, GO, NT>
      (*sourceMap,
       targetMapRemoteOrPermuteGlobalIndices,
       targetMapRemoteOrPermuteProcessRanks,
       numTargetMapRemoteOrPermuteGlobalIndices,
       mayReorderTargetMapIndicesLocally,
       this->TransferData_->out_.getRawPtr (),
       verbPfx.get (),
       this->verbose (),
       debug);

    // Since we invoked the base_type constructor above, we know that
    // out_ is nonnull, so we don't have to waste time creating it
    // again.
    using data_type = ImportExportData<LO, GO, NT>;
    TEUCHOS_ASSERT( ! this->TransferData_.is_null () );
    this->TransferData_ = rcp (new data_type (sourceMap,
                                              localSetupResult.targetMap,
                                              this->TransferData_->out_,
                                              plist));
    this->TransferData_->numSameIDs_ = localSetupResult.numSameIDs;
    // Skip permutes; they are user error, because they duplicate
    // non-remote indices.
    makeDualViewFromVector (this->TransferData_->remoteLIDs_,
                            localSetupResult.remoteLIDs,
                            "remoteLIDs");
    // "Is locally complete" for an Import means that all target Map
    // indices on the calling process exist on at least one process
    // (not necessarily this one) in the source Map.  For this
    // constructor, this is true if and only if all input target PIDs
    // are valid PIDs in the communicator.
    //
    // FIXME (mfh 20 Feb 2018) For now, assume this is always true.
    this->TransferData_->isLocallyComplete_ = true;

    Teuchos::Array<GO> exportGIDs;
    if (sourceMap->isDistributed ()) {
      if (this->verbose ()) {
        std::ostringstream os;
        os << *verbPfx << "Make Distributor (createFromRecvs)" << endl;
        this->verboseOutputStream () << os.str ();
      }
      ArrayView<const GO> remoteGIDs (localSetupResult.remoteGIDs.data (),
                                      localSetupResult.remoteGIDs.size ());
      ArrayView<const int> remotePIDs (localSetupResult.remotePIDs.data (),
                                       localSetupResult.remotePIDs.size ());
      // Call Distributor::createFromRecvs to turn the remote GIDs and
      // their owning PIDs into a send-and-receive communication plan.
      // remoteGIDs and remotePIDs are input; exportGIDs and
      // exportPIDs are output arrays that createFromRecvs allocates.
      Distributor& distributor = this->TransferData_->distributor_;
      distributor.createFromRecvs (remoteGIDs,
                                   remotePIDs,
                                   exportGIDs,
                                   this->TransferData_->exportPIDs_);
      // Find the LIDs corresponding to the (outgoing) GIDs in
      // exportGIDs.  For sparse matrix-vector multiply, this tells
      // the calling process how to index into the source vector to
      // get the elements which it needs to send.
      //
      // NOTE (mfh 03 Mar 2014) This is now a candidate for a
      // thread-parallel kernel, but only if using the new thread-safe
      // Map implementation.
      if (this->verbose ()) {
        std::ostringstream os;
        os << *verbPfx << "Compute exportLIDs" << endl;
        this->verboseOutputStream () << os.str ();
      }
      using size_type = typename Teuchos::Array<GO>::size_type;
      const size_type numExportIDs = exportGIDs.size ();

      typename decltype (this->TransferData_->exportLIDs_)::t_host
        exportLIDs (view_alloc_no_init ("exportLIDs"), numExportIDs);

      for (size_type k = 0; k < numExportIDs; ++k) {
        exportLIDs[k] = sourceMap->getLocalElement (exportGIDs[k]);
      }
      makeDualViewFromOwningHostView (this->TransferData_->exportLIDs_, exportLIDs);
    }

    if (this->verbose ()) {
      std::ostringstream os;
      os << *verbPfx;
      printDualView (os, this->TransferData_->remoteLIDs_,
                     "ImportExportData::remoteLIDs_");
      os << endl;
      this->verboseOutputStream () << os.str ();
    }
    if (this->verbose ()) {
      std::ostringstream os;
      os << *verbPfx << "Done!" << endl;
      this->verboseOutputStream () << os.str ();
    }
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
    using ::Tpetra::Details::makeDualViewFromOwningHostView;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::view_alloc_no_init;
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::null;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;
    ProfilingRegion regionExport ("Tpetra::Import::setupSamePermuteRemote");

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
    // numSameIDs_ is for the common case of an Import where all the
    // overlapping GIDs are at the end of the target Map, but
    // otherwise the source and target Maps are the same.  This allows
    // a fast contiguous copy for the initial "same IDs."
    size_type numSameGids = 0;
    for ( ; numSameGids < numGids && rawSrcGids[numSameGids] == rawTgtGids[numSameGids]; ++numSameGids)
      {} // third clause of 'for' does everything
    this->TransferData_->numSameIDs_ = numSameGids;

    // Compute permuteToLIDs_, permuteFromLIDs_, remoteGIDs, and
    // remoteLIDs_.  The first two arrays are IDs to be permuted, and
    // the latter two arrays are IDs to be received ("imported"),
    // called "remote" IDs.
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  IDs to receive are in the
    // target Map, but not the source Map.

    // Iterate over the target Map's LIDs, since we only need to do
    // GID -> LID lookups for the source Map.
    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO numTgtLids = as<LO> (numTgtGids);
    LO numPermutes = 0;

    for (LO tgtLid = numSameGids; tgtLid < numTgtLids; ++tgtLid) {
      const GO curTargetGid = rawTgtGids[tgtLid];
      // getLocalElement() returns LINVALID if the GID isn't in the
      // source Map.  This saves us a lookup (which
      // isNodeGlobalElement() would do).
      const LO srcLid = source.getLocalElement (curTargetGid);
      if (srcLid != LINVALID) { // if source.isNodeGlobalElement (curTargetGid)
        ++numPermutes;
      }
    }
    const LO numRemotes = (numTgtLids - numSameGids) - numPermutes;

    using host_perm_type =
      typename decltype (this->TransferData_->permuteToLIDs_)::t_host;
    host_perm_type permuteToLIDs
      (view_alloc_no_init ("permuteToLIDs"), numPermutes);
    host_perm_type permuteFromLIDs
      (view_alloc_no_init ("permuteFromLIDs"), numPermutes);
    typename decltype (this->TransferData_->remoteLIDs_)::t_host remoteLIDs
      (view_alloc_no_init ("permuteFromLIDs"), numRemotes);

    {
      LO numPermutes2 = 0;
      LO numRemotes2 = 0;
      for (LO tgtLid = numSameGids; tgtLid < numTgtLids; ++tgtLid) {
        const GO curTargetGid = rawTgtGids[tgtLid];
        const LO srcLid = source.getLocalElement (curTargetGid);
        if (srcLid != LINVALID) {
          permuteToLIDs[numPermutes2] = tgtLid;
          permuteFromLIDs[numPermutes2] = srcLid;
          ++numPermutes2;
        }
        else {
          remoteGIDs.push_back (curTargetGid);
          remoteLIDs[numRemotes2] = tgtLid;
          ++numRemotes2;
        }
      }
      TEUCHOS_ASSERT( numPermutes == numPermutes2 );
      TEUCHOS_ASSERT( numRemotes == numRemotes2 );
      TEUCHOS_ASSERT( size_t (numPermutes) + remoteGIDs.size () == size_t (numTgtLids - numSameGids) );
    }

    makeDualViewFromOwningHostView (this->TransferData_->permuteToLIDs_, permuteToLIDs);
    makeDualViewFromOwningHostView (this->TransferData_->permuteFromLIDs_, permuteFromLIDs);
    makeDualViewFromOwningHostView (this->TransferData_->remoteLIDs_, remoteLIDs);
    if (remoteLIDs.extent (0) != 0 && ! source.isDistributed ()) {
      // This Import has remote LIDs, meaning that the target Map has
      // entries on this process that are not in the source Map on
      // this process.  However, the source Map is not distributed
      // globally.  This implies that this Import is not locally
      // complete on this process.
      this->TransferData_->isLocallyComplete_ = false;
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
               Teuchos::Array<int>& userRemotePIDs,
               const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using ::Tpetra::Details::makeDualViewFromOwningHostView;
    using ::Tpetra::Details::view_alloc_no_init;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using std::endl;
    using GO = GlobalOrdinal;
    typedef typename Array<int>::difference_type size_type;
    const char tfecfFuncName[] = "setupExport: ";
    const char suffix[] = "  Please report this bug to the Tpetra developers.";

    std::unique_ptr<std::string> prefix;
    if (this->verbose ()) {
      auto srcMap = this->getSourceMap ();
      auto comm = srcMap.is_null () ? Teuchos::null : srcMap->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::Import::setupExport: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << std::endl;
      this->verboseOutputStream () << os.str ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->getSourceMap ().is_null (), std::logic_error,
       "Source Map is null.  " << suffix);
    const map_type& source = * (this->getSourceMap ());

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
      if (this->verbose ()) {
        std::ostringstream os;
        os << *prefix << "Call sourceMap.getRemoteIndexList" << endl;
        this->verboseOutputStream () << os.str ();
      }
      lookup = source.getRemoteIndexList (remoteGIDsView, newRemotePIDs ());
    }
    Array<int>& remoteProcIDs = useRemotePIDs ? userRemotePIDs : newRemotePIDs;

    if (lookup == IDNotPresent) {
      // There is at least one GID owned by the calling process in the
      // target Map, which is not owned by any process in the source
      // Map.
      this->TransferData_->isLocallyComplete_ = false;

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
                       std::bind (std::equal_to<int> (), -1, std::placeholders::_1));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numInvalidRemote == 0, std::logic_error, "Calling getRemoteIndexList "
         "on the source Map returned IDNotPresent, but none of the returned "
         "\"remote\" process ranks are -1.  Please report this bug to the "
         "Tpetra developers.");

#ifdef HAVE_TPETRA_MMM_TIMINGS
      using Teuchos::TimeMonitor;
      std::string label;
      if(!plist.is_null())
        label = plist->get("Timer Label",label);
      std::string prefix = std::string("Tpetra ")+ label + std::string(":iport_ctor:setupExport:1 ");
      auto MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix)));
#else
      (void)plist;
#endif

      // If all of them are invalid, we can delete the whole array.
      const size_type totalNumRemote = this->getNumRemoteIDs ();
      if (numInvalidRemote == totalNumRemote) {
        // all remotes are invalid; we have no remotes; we can delete the remotes
        remoteProcIDs.clear ();
        remoteGIDs.clear (); // invalidates remoteGIDsView
        this->TransferData_->remoteLIDs_ =
          decltype (this->TransferData_->remoteLIDs_) ();
      }
      else {
        // Some remotes are valid; we need to keep the valid ones.
        // Pack and resize remoteProcIDs, remoteGIDs, and remoteLIDs_.
        size_type numValidRemote = 0;
#ifdef HAVE_TPETRA_DEBUG
        ArrayView<GO> remoteGIDsPtr = remoteGIDsView;
#else
        GO* const remoteGIDsPtr = remoteGIDsView.getRawPtr ();
#endif // HAVE_TPETRA_DEBUG

        // Don't mark the DualView modified, since we'll reallocate it.
        auto remoteLIDs = this->TransferData_->remoteLIDs_.view_host ();

        for (size_type r = 0; r < totalNumRemote; ++r) {
          // Pack in all the valid remote PIDs and GIDs.
          if (remoteProcIDs[r] != -1) {
            remoteProcIDs[numValidRemote] = remoteProcIDs[r];
            remoteGIDsPtr[numValidRemote] = remoteGIDsPtr[r];
            remoteLIDs[numValidRemote] = remoteLIDs[r];
            ++numValidRemote;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numValidRemote != totalNumRemote - numInvalidRemote,
           std::logic_error, "After removing invalid remote GIDs and packing "
           "the valid remote GIDs, numValidRemote = " << numValidRemote
           << " != totalNumRemote - numInvalidRemote = "
           << totalNumRemote - numInvalidRemote
           << ".  Please report this bug to the Tpetra developers.");

        remoteProcIDs.resize (numValidRemote);
        remoteGIDs.resize (numValidRemote);

        Kokkos::resize (remoteLIDs, numValidRemote);
        this->TransferData_->remoteLIDs_ = decltype (this->TransferData_->remoteLIDs_) ();
        makeDualViewFromOwningHostView (this->TransferData_->remoteLIDs_, remoteLIDs);
      }
      // Revalidate the view after clear or resize.
      remoteGIDsView = remoteGIDs ();
    }

    // Sort remoteProcIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs and remoteLIDs_.  This ensures that
    // remoteProcIDs[i], remoteGIDs[i], and remoteLIDs_[i] all refer
    // to the same thing.
    {
      this->TransferData_->remoteLIDs_.modify_host ();
      auto remoteLIDs = this->TransferData_->remoteLIDs_.view_host ();
      sort3 (remoteProcIDs.begin (),
             remoteProcIDs.end (),
             remoteGIDsView.getRawPtr (),
             remoteLIDs.data ());
      this->TransferData_->remoteLIDs_.sync_device ();
    }

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDs and remoteProcIDs_ are input;
    // exportGIDs and exportProcIDs_ are output arrays which are
    // allocated by createFromRecvs().
    Array<GO> exportGIDs;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    using Teuchos::TimeMonitor;
    std::string label;
    if(!plist.is_null())
      label = plist->get("Timer Label",label);
    std::string prefix2 = std::string("Tpetra ")+ label + std::string(":iport_ctor:setupExport:3 ");
    auto MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2)));
#endif

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Call createFromRecvs" << endl;
      this->verboseOutputStream () << endl;
    }
    this->TransferData_->distributor_.createFromRecvs (remoteGIDsView ().getConst (),
                                               remoteProcIDs, exportGIDs,
                                               this->TransferData_->exportPIDs_);

    // Find the LIDs corresponding to the (outgoing) GIDs in
    // exportGIDs.  For sparse matrix-vector multiply, this tells the
    // calling process how to index into the source vector to get the
    // elements which it needs to send.
#ifdef HAVE_TPETRA_MMM_TIMINGS
    prefix2 = std::string("Tpetra ")+ label + std::string(":iport_ctor:setupExport:4 ");
    MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2)));
#endif

    // NOTE (mfh 03 Mar 2014) This is now a candidate for a
    // thread-parallel kernel, but only if using the new thread-safe
    // Map implementation.
    const size_type numExportIDs = exportGIDs.size ();
    if (numExportIDs > 0) {
      typename decltype (this->TransferData_->exportLIDs_)::t_host
        exportLIDs (view_alloc_no_init ("exportLIDs"), numExportIDs);
      ArrayView<const GO> expGIDs = exportGIDs ();

      for (size_type k = 0; k < numExportIDs; ++k) {
        exportLIDs[k] = source.getLocalElement (expGIDs[k]);
      }
      makeDualViewFromOwningHostView (this->TransferData_->exportLIDs_, exportLIDs);
    }

    if (this->verbose ()) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl;
      this->verboseOutputStream () << os.str ();
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
    using ::Tpetra::Details::Behavior;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using GST = Tpetra::global_size_t;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using import_type = Import<LO, GO, Node>;
    using size_type = typename Array<GO>::size_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    using Teuchos::TimeMonitor;
    std::string label = std::string("Tpetra::Import::setUnion");
    TimeMonitor MM(*TimeMonitor::getNewTimer(label));
#endif

    RCP<const map_type> srcMap = this->getSourceMap ();
    RCP<const map_type> tgtMap1 = this->getTargetMap ();
    RCP<const map_type> tgtMap2 = rhs.getTargetMap ();
    RCP<const Comm<int> > comm = srcMap->getComm ();

    const bool debug = Behavior::debug ("Import::setUnion") ||
      Behavior::debug ("Tpetra::Import::setUnion");

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! srcMap->isSameAs (* (rhs.getSourceMap ())), std::invalid_argument,
         "Tpetra::Import::setUnion: The source Map of the input Import must be the "
         "same as (in the sense of Map::isSameAs) the source Map of this Import.");
      const Comm<int>& comm1 = * (tgtMap1->getComm ());
      const Comm<int>& comm2 = * (tgtMap2->getComm ());
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ::Tpetra::Details::congruent (comm1, comm2),
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
    ArrayView<const GO> sameGIDs1 = (tgtMap1->getLocalElementList())(0,numSameGIDs1);

    const size_type numSameGIDs2 = rhs.getNumSameIDs();
    ArrayView<const GO> sameGIDs2 = (tgtMap2->getLocalElementList())(0,numSameGIDs2);

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
    MM.disableTimer(label);
    label = "Tpetra::Import::setUnion : Construct Target Map";
    TimeMonitor MM2(*TimeMonitor::getNewTimer(label));
#endif

    // Create the union target Map.
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBaseUnion = std::min(tgtMap1->getIndexBase(), tgtMap2->getIndexBase());
    RCP<const map_type> unionTgtMap =
      rcp(new map_type(INVALID, unionTgtGIDs(), indexBaseUnion, comm));

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM2.disableTimer(label);
    label = "Tpetra::Import::setUnion : Export GIDs";
    TimeMonitor MM3(*TimeMonitor::getNewTimer(label));
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
    Distributor distributor (comm, this->TransferData_->out_);

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
    MM3.disableTimer(label);
    label = "Tpetra::Import::setUnion : Construct Import";
    TimeMonitor MM4(*TimeMonitor::getNewTimer(label));
#endif
    RCP<const import_type> unionImport =
      rcp (new import_type (srcMap, unionTgtMap,
                            as<size_t> (numSameIDsUnion),
                            permuteToLIDsUnion, permuteFromLIDsUnion,
                            remoteLIDsUnion, exportLIDsUnion,
                            exportPIDsUnion, distributor,
                            this->TransferData_->out_));
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

    ArrayView<const GO> srcGIDs = srcMap->getLocalElementList ();
    ArrayView<const GO> tgtGIDs = tgtMap->getLocalElementList ();

    // All elements in srcMap will be in the "new" target map, so...
    size_t numSameIDsNew    = srcMap->getLocalNumElements ();
    size_t numRemoteIDsNew  = this->getNumRemoteIDs ();
    Array<LO> permuteToLIDsNew, permuteFromLIDsNew; // empty on purpose

    // Grab some old data
    ArrayView<const LO> remoteLIDsOld = this->getRemoteLIDs ();
    ArrayView<const LO> exportLIDsOld = this->getExportLIDs ();

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
    RCP<const map_type> targetMapNew =
      rcp (new map_type (GO_INVALID, GIDs, tgtMap->getIndexBase (),
                         tgtMap->getComm ()));

    // Exports are trivial (since the sourcemap doesn't change)
    Array<int> exportPIDsnew (this->getExportPIDs ());
    Array<LO> exportLIDsnew (this->getExportLIDs ());

    // Copy the Distributor (due to how the Import constructor works)
    Distributor D (this->getDistributor ());

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
  createRemoteOnlyImport (const Teuchos::RCP<const map_type>& remoteTarget) const
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::gathervPrint;
    using Teuchos::outArg;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using import_type = Import<LocalOrdinal,GlobalOrdinal,Node>;

    const char funcPrefix[] = "Tpetra::createRemoteOnlyImport: ";
    int lclSuccess = 1;
    int gblSuccess = 1;
    const bool debug = Behavior::debug ();

    const size_t NumRemotes = this->getNumRemoteIDs ();
    std::unique_ptr<std::string> procPrefix;
    Teuchos::RCP<const Teuchos::Comm<int>> comm;
    if (debug) {
      comm = remoteTarget.is_null () ? Teuchos::null :
        remoteTarget->getComm ();
      std::ostringstream os;
      os << "Proc ";
      if (comm.is_null ()) {
        os << "?";
      }
      else {
        os << comm->getRank ();
      }
      os << ": ";
      procPrefix = std::unique_ptr<std::string> (new std::string (os.str ()));
    }

    if (debug) {
      std::ostringstream lclErr;
      if (remoteTarget.is_null ()) {
        lclSuccess = -1;
      }
      else if (NumRemotes != remoteTarget->getLocalNumElements ()) {
        lclSuccess = 0;
        lclErr << *procPrefix << "getNumRemoteIDs() = " << NumRemotes
               << " != remoteTarget->getLocalNumElements() = "
               << remoteTarget->getLocalNumElements () << "." << endl;
      }

      if (comm.is_null ()) {
        lclSuccess = gblSuccess;
      }
      else {
        reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (gblSuccess == -1, std::invalid_argument, funcPrefix
         << "Input target Map is null on at least one process.");

      if (gblSuccess != 1) {
        if (comm.is_null ()) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, lclErr.str ());
        }
        else {
          std::ostringstream gblErr;
          gblErr << funcPrefix << endl;
          gathervPrint (gblErr, lclErr.str (), *comm);
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, gblErr.str ());
        }
      }
    }

    // Compute the new Remote LIDs
    Teuchos::ArrayView<const LO> oldRemoteLIDs = this->getRemoteLIDs ();
    Teuchos::Array<LO> newRemoteLIDs (NumRemotes);
    const map_type& tgtMap = * (this->getTargetMap ());
    size_t badCount = 0;

    std::unique_ptr<std::vector<size_t>> badIndices;
    if (debug) {
      badIndices = std::unique_ptr<std::vector<size_t>> (new std::vector<size_t>);
    }

    for (size_t i = 0; i < NumRemotes; ++i) {
      const LO oldLclInd = oldRemoteLIDs[i];
      if (oldLclInd == Teuchos::OrdinalTraits<LO>::invalid ()) {
        ++badCount;
        if (debug) { badIndices->push_back (i); }
        continue;
      }
      const GO gblInd = tgtMap.getGlobalElement (oldLclInd);
      if (gblInd == Teuchos::OrdinalTraits<GO>::invalid ()) {
        ++badCount;
        if (debug) { badIndices->push_back (i); }
        continue;
      }
      const LO newLclInd = remoteTarget->getLocalElement (gblInd);
      if (newLclInd == Teuchos::OrdinalTraits<LO>::invalid ()) {
        ++badCount;
        if (debug) { badIndices->push_back (i); }
        continue;
      }
      newRemoteLIDs[i] = newLclInd;
      // Now we make sure these guys are in sorted order (AztecOO-ML ordering)
      if (i > 0 && newRemoteLIDs[i] < newRemoteLIDs[i-1]) {
        ++badCount;
        if (debug) { badIndices->push_back (i); }
      }
    }

    if (badCount != 0) {
      lclSuccess = 0;
    }

    if (debug) {
      if (comm.is_null ()) {
        lclSuccess = gblSuccess;
      }
      else {
        reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      }
      std::ostringstream lclErr;
      if (lclSuccess != 1) {
        lclErr << *procPrefix << "Count of bad indices: " << badCount
               << ", bad indices: [";
        // TODO (mfh 04 Sep 2019) Limit the maximum number of bad
        // indices to print.
        for (size_t k = 0; k < badCount; ++k) {
          const size_t badIndex = (*badIndices)[k];
          lclErr << "(" << badIndex << ","
                 << oldRemoteLIDs[badIndex] << ")";
          if (k + size_t (1) < badCount) {
            lclErr << ", ";
          }
        }
        lclErr << "]" << endl;
      }

      if (gblSuccess != 1) {
        std::ostringstream gblErr;
        gblErr << funcPrefix << "this->getRemoteLIDs() has \"bad\" "
          "indices on one or more processes.  \"Bad\" means that the "
          "indices are invalid, they don't exist in the target Map, "
          "they don't exist in remoteTarget, or they are not in "
          "sorted order.  In what follows, I will show the \"bad\" "
          "indices as (k, LID) pairs, where k is the zero-based "
          "index of the LID in this->getRemoteLIDs()." << endl;
        if (comm.is_null ()) {
          gblErr << lclErr.str ();
        }
        else {
          gathervPrint (gblErr, lclErr.str (), *comm);
        }
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::runtime_error, gblErr.str ());
      }
    }
    else { // not debug
      TEUCHOS_TEST_FOR_EXCEPTION
        (lclSuccess == 0, std::runtime_error, funcPrefix
         << "this->getRemoteLIDs() has " << badCount
         << "ind" << (badCount == 1 ? "ex" : "ices")
         << " \"bad\" indices on this process." << endl);
    }

    // Copy ExportPIDs and such
    // NOTE: Be careful: The "Expert" Import constructor we use does a "swap"
    // for most of the LID/PID lists and the Distributor, meaning it
    // ruins the existing object if we pass things in directly.  Hence
    // we copy them first.
    Teuchos::Array<int> newExportPIDs (this->getExportPIDs ());
    Teuchos::Array<LO> newExportLIDs (this->getExportLIDs ());
    Teuchos::Array<LO> dummy;
    Distributor newDistor (this->getDistributor ());

    return rcp (new import_type (this->getSourceMap (), remoteTarget,
                                 static_cast<size_t> (0), dummy, dummy,
                                 newRemoteLIDs, newExportLIDs,
                                 newExportPIDs, newDistor));
  }

} // namespace Tpetra

#define TPETRA_IMPORT_CLASS_INSTANT(LO, GO, NODE) \
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
