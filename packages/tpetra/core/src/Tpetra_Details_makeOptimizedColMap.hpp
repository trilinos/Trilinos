// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_makeOptimizedColMap.hpp
/// \brief Optimize your graph's or matrix's column Map for better
///   communication performance.
///
/// \warning Any function (indeed, anything) in the Tpetra::Details
///   namespace comes with <i>no guarantee</i> of backwards
///   compatibility.  It may disappear or change at any time.
///
/// The functions in this file take a domain Map and a column Map of a
/// distributed graph (Tpetra::CrsGraph) or matrix (e.g.,
/// Tpetra::CrsMatrix).  They then create a new column Map, which
/// optimizes the performance of an Import operation from the domain
/// Map to the new column Map.  The function makeOptimizedColMapAndImport
/// also creates and returns that Import.  Creating the new column Map
/// and its Import at the same time saves some communication, since
/// making the Import requires some of the same information that
/// optimizing the column Map does.

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <memory>

namespace Tpetra {
namespace Details {

  /// \class OptColMap
  /// \brief Implementation detail of makeOptimizedColMap, and
  ///   makeOptimizedColMapAndImport.
  ///
  /// \warning Please use the free functions (not part of a class)
  ///   makeOptimizedColMap or makeOptimizedColMapAndImport.
  template<class MapType>
  class OptColMap {
  public:
    using local_ordinal_type = typename MapType::local_ordinal_type;
    using global_ordinal_type = typename MapType::global_ordinal_type;
    using node_type = typename MapType::node_type;
    using map_type = ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;
    using import_type = ::Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>;

    static Teuchos::RCP<const map_type>
    makeOptColMap (std::ostream& errStream,
                   bool& lclErr,
                   const map_type& domMap,
                   const map_type& colMap,
                   const import_type* /* oldImport */)
    {
      using ::Tpetra::Details::Behavior;
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::FancyOStream;
      using Teuchos::getFancyOStream;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      using std::endl;
      using LO = local_ordinal_type;
      using GO = global_ordinal_type;
      const char prefix[] = "Tpetra::Details::makeOptimizedColMap: ";

      RCP<const Teuchos::Comm<int> > comm = colMap.getComm ();
      std::ostream& err = errStream;

      const bool verbose = Behavior::verbose ("Tpetra::Details::makeOptimizedColMap");

      RCP<FancyOStream> outPtr = getFancyOStream (rcpFromRef (std::cerr));
      TEUCHOS_TEST_FOR_EXCEPTION
        (outPtr.is_null (), std::logic_error,
         "outPtr is null; this should never happen!");
      FancyOStream& out = *outPtr;
      Teuchos::OSTab tab1 (out);

      std::unique_ptr<std::string> verboseHeader;
      if (verbose) {
        std::ostringstream os;
        const int myRank = comm->getRank ();
        os << "Proc " << myRank << ": ";
        verboseHeader = std::unique_ptr<std::string> (new std::string (os.str ()));
      }
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "Tpetra::Details::makeOptimizedColMap" << endl;
        out << os.str ();
      }

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "Domain Map GIDs: [";
        const LO domMapLclNumInds = static_cast<LO> (domMap.getLocalNumElements ());
        for (LO lid = 0; lid < domMapLclNumInds; ++lid) {
          const GO gid = domMap.getGlobalElement (lid);
          os << gid;
          if (lid + LO (1) < domMapLclNumInds) {
            os << ", ";
          }
        }
        os << "]" << endl;
        out << os.str ();
      }

      const LO colMapLclNumInds = static_cast<LO> (colMap.getLocalNumElements ());

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "Column Map GIDs: [";
        for (LO lid = 0; lid < colMapLclNumInds; ++lid) {
          const GO gid = colMap.getGlobalElement (lid);
          os << gid;
          if (lid + LO (1) < colMapLclNumInds) {
            os << ", ";
          }
        }
        os << "]" << endl;
        out << os.str ();
      }

      // Count remote GIDs.
      LO numOwnedGids = 0;
      LO numRemoteGids = 0;
      for (LO colMapLid = 0; colMapLid < colMapLclNumInds; ++colMapLid) {
        const GO colMapGid = colMap.getGlobalElement (colMapLid);
        if (domMap.isNodeGlobalElement (colMapGid)) {
          ++numOwnedGids;
        }
        else {
          ++numRemoteGids;
        }
      }

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- numOwnedGids: " << numOwnedGids << endl
           << *verboseHeader << "- numRemoteGids: " << numRemoteGids << endl;
        out << os.str ();
      }

      // Put all colMap GIDs on the calling process in a single array.
      // Owned GIDs go in front, and remote GIDs at the end.
      Array<GO> allGids (numOwnedGids + numRemoteGids);
      ArrayView<GO> ownedGids = allGids.view (0, numOwnedGids);
      ArrayView<GO> remoteGids = allGids.view (numOwnedGids, numRemoteGids);

      // Fill ownedGids and remoteGids (and therefore allGids).  We use
      // two loops, one to count (above) and one to fill (here), in
      // order to avoid dynamic memory allocation during the loop (in
      // this case, lots of calls to push_back()).  That will simplify
      // use of Kokkos to parallelize these loops later.
      LO ownedPos = 0;
      LO remotePos = 0;
      for (LO colMapLid = 0; colMapLid < colMapLclNumInds; ++colMapLid) {
        const GO colMapGid = colMap.getGlobalElement (colMapLid);
        if (domMap.isNodeGlobalElement (colMapGid)) {
          ownedGids[ownedPos++] = colMapGid;
        }
        else {
          remoteGids[remotePos++] = colMapGid;
        }
      }

      // If, for some reason, the running count doesn't match the
      // orignal count, fill in any remaining GID spots with an
      // obviously invalid value.  We don't want to stop yet, because
      // other processes might not have noticed this error; Map
      // construction is a collective, so we can't stop now.
      if (ownedPos != numOwnedGids) {
        lclErr = true;
        err << prefix << "On Process " << comm->getRank () << ", ownedPos = "
            << ownedPos << " != numOwnedGids = " << numOwnedGids << endl;
        for (LO colMapLid = ownedPos; colMapLid < numOwnedGids; ++colMapLid) {
          ownedGids[colMapLid] = Teuchos::OrdinalTraits<GO>::invalid ();
        }
      }
      if (remotePos != numRemoteGids) {
        lclErr = true;
        err << prefix << "On Process " << comm->getRank () << ", remotePos = "
            << remotePos << " != numRemoteGids = " << numRemoteGids << endl;
        for (LO colMapLid = remotePos; colMapLid < numRemoteGids; ++colMapLid) {
          remoteGids[colMapLid] = Teuchos::OrdinalTraits<GO>::invalid ();
        }
      }

      // Figure out what processes own what GIDs in the domain Map.
      // Initialize the output array of remote PIDs with the "invalid
      // process rank" -1, to help us test whether getRemoteIndexList
      // did its job.
      Array<int> remotePids (numRemoteGids, -1);
      const LookupStatus lookupStatus =
        domMap.getRemoteIndexList (remoteGids, remotePids ());

      // If any process returns IDNotPresent, then at least one of the
      // remote indices was not present in the domain Map.  This means
      // that the Import object cannot be constructed, because of
      // incongruity between the column Map and domain Map.  This
      // means that either the column Map or domain Map, or both, is
      // incorrect.
      const bool getRemoteIndexListFailed = (lookupStatus == IDNotPresent);
      if (getRemoteIndexListFailed) {
        lclErr = true;
        err << prefix << "On Process " << comm->getRank () << ", some indices "
          "in the input colMap (the original column Map) are not in domMap (the "
          "domain Map).  Either these indices or the domain Map is invalid.  "
          "Likely cause: For a nonsquare matrix, you must give the domain and "
          "range Maps as input to fillComplete." << endl;
      }

      // Check that getRemoteIndexList actually worked, by making sure
      // that none of the remote PIDs are -1.
      for (LO k = 0; k < numRemoteGids; ++k) {
        bool foundInvalidPid = false;
        if (remotePids[k] == -1) {
          foundInvalidPid = true;
          break;
        }
        if (foundInvalidPid) {
          lclErr = true;
          err << prefix << "On Process " << comm->getRank () << ", "
            "getRemoteIndexList returned -1 for the process ranks of "
            "one or more GIDs on this process." << endl;
        }
      }

      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- Before sort2:" << endl
           << *verboseHeader << "-- ownedGids: " << Teuchos::toString (ownedGids) << endl
           << *verboseHeader << "-- remoteGids: " << Teuchos::toString (remoteGids) << endl
           << *verboseHeader << "-- allGids: " << Teuchos::toString (allGids ()) << endl;
        out << os.str ();
      }
      using Tpetra::sort2;
      sort2 (remotePids.begin (), remotePids.end (), remoteGids.begin (), true);
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "- After sort2:" << endl
           << *verboseHeader << "-- ownedGids: " << Teuchos::toString (ownedGids) << endl
           << *verboseHeader << "-- remoteGids: " << Teuchos::toString (remoteGids) << endl
           << *verboseHeader << "-- allGids: " << Teuchos::toString (allGids ()) << endl;
        out << os.str ();
      }

      auto optColMap = rcp (new map_type (colMap.getGlobalNumElements (),
                                          allGids (),
                                          colMap.getIndexBase (),
                                          comm));
      if (verbose) {
        std::ostringstream os;
        os << *verboseHeader << "Tpetra::Details::makeOptimizedColMap: Done" << endl;
        out << os.str ();
      }
      return optColMap;
    }

    /// \brief Return an optimized reordering of the given column Map.
    ///   Optionally, recompute an Import from the input domain Map to
    ///   the new column Map.
    /// \tparam MapType A specialization of Map.
    ///
    /// See the documentation of the free function
    /// makeOptimizedColMapAndImport().
    ///
    /// \param errStream [out] Output stream for human-readable error
    ///   reporting.  This is local to the calling process and may
    ///   differ on different processes.
    /// \param lclErr [out] On output: true if anything went wrong on
    ///   the calling process.  This value is local to the calling
    ///   process and may differ on different processes.
    /// \param domMap [in] Domain Map of a CrsGraph or CrsMatrix.
    /// \param colMap [in] <i>Original</i> column Map of the same
    ///   CrsGraph or CrsMatrix as \c domMap.
    /// \param oldImport [in] Optional pointer to the "original
    ///   Import: an Import from \c domMap to \c colMap.  This is not
    ///   required, but if you supply this, this function may use it
    ///   to avoid some communication and/or work when setting up the
    ///   new Import object. 
    ///
    /// \return The possibly reordered column Map \c newColMap, and the
    ///   corresponding Import from \c domMap to \c newColMap. 
    ///
    /// \pre \c domMap and \c colMap must have the same or congruent
    ///   communicators.
    /// \pre On all calling processes, the indices in \c colMap must be
    ///   a subset of the indices in \c domMap.
    static std::pair<Teuchos::RCP<const map_type>,
                     Teuchos::RCP<import_type> >
    makeOptColMapAndImport (std::ostream& errStream,
                            bool& lclErr,
                            const map_type& domMap,
                            const map_type& colMap,
                            const import_type* oldImport)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;

      // mfh 15 May 2018: For now, just call makeOptColMap, and use
      // the conventional two-Map (source and target) Import
      // constructor.
      RCP<const map_type> newColMap =
        makeOptColMap (errStream, lclErr, domMap, colMap, oldImport);
      RCP<import_type> imp (new import_type (rcp (new map_type (domMap)), newColMap));

      // FIXME (mfh 06 Jul 2014) This constructor throws a runtime
      // error, so I'm not using it for now.
      //
      // imp = rcp (new import_type (domMap, newColMap, remoteGids,
      //                             remotePids (), remoteLids (),
      //                             Teuchos::null, Teuchos::null));
      return std::make_pair (newColMap, imp);
    }
  };

  /// \brief Return an optimized reordering of the given column Map.
  /// \tparam MapType A specialization of Map.
  ///
  /// \param err [out] Output stream for human-readable error
  ///   reporting.  This is local to the calling process and may
  ///   differ on different processes.
  /// \param lclErr [out] On output: true if anything went wrong on
  ///   the calling process.  This value is local to the calling
  ///   process and may differ on different processes.
  /// \param domMap [in] Domain Map of a CrsGraph or CrsMatrix.
  /// \param colMap [in] <i>Original</i> column Map of the same
  ///   CrsGraph or CrsMatrix as \c domMap.
  /// \param oldImport [in] Optional pointer to the "original
  ///   Import: an Import from \c domMap to \c colMap.  This is not
  ///   required, but if you supply this, this function may use it
  ///   to avoid some communication and/or work when setting up the
  ///   new Import object. 
  ///
  /// \return The possibly reordered column Map \c newColMap.
  ///
  /// This is a convenience wrapper for makeOptimizedColMapAndImport().
  /// (Please refer to that function's documentation in this file.)
  /// It does everything that that function does, except that it does
  /// not compute a new Import.
  template<class MapType>
  Teuchos::RCP<const MapType>
  makeOptimizedColMap (std::ostream& errStream,
                       bool& lclErr,
                       const MapType& domMap,
                       const MapType& colMap,
                       const Tpetra::Import<
                         typename MapType::local_ordinal_type,
                         typename MapType::global_ordinal_type,
                         typename MapType::node_type>* oldImport = nullptr)
  {
    using map_type = ::Tpetra::Map<
      typename MapType::local_ordinal_type,
      typename MapType::global_ordinal_type,
      typename MapType::node_type>;
    using impl_type = OptColMap<map_type>;
    auto mapPtr = impl_type::makeOptColMap (errStream, lclErr,
                                            domMap, colMap, oldImport);
    return mapPtr;
  }

  /// \brief Return an optimized reordering of the given column Map.
  ///   Optionally, recompute an Import from the input domain Map to
  ///   the new column Map.
  /// \tparam MapType A specialization of Map.
  ///
  /// This function takes a domain Map and a column Map of a
  /// distributed graph (Tpetra::CrsGraph) or matrix (e.g.,
  /// Tpetra::CrsMatrix).  It then creates a new column Map, which
  /// optimizes the performance of an Import operation from the domain
  /// Map to the new column Map.  This function also optionally
  /// creates that Import.  Creating the new column Map and its Import
  /// at the same time saves some communication, since making the
  /// Import requires some of the same information that optimizing the
  /// column Map does.
  ///
  /// \param err [out] Output stream for human-readable error
  ///   reporting.  This is local to the calling process and may
  ///   differ on different processes.
  /// \param lclErr [out] On output: true if anything went wrong on
  ///   the calling process.  This value is local to the calling
  ///   process and may differ on different processes.
  /// \param domMap [in] Domain Map of a CrsGraph or CrsMatrix.
  /// \param colMap [in] <i>Original</i> column Map of the same
  ///   CrsGraph or CrsMatrix as \c domMap.
  /// \param oldImport [in] Optional pointer to the "original
  ///   Import: an Import from \c domMap to \c colMap.  This is not
  ///   required, but if you supply this, this function may use it
  ///   to avoid some communication and/or work when setting up the
  ///   new Import object. 
  ///
  /// \return The possibly reordered column Map \c newColMap, and the
  ///   corresponding Import from \c domMap to \c newColMap. 
  ///
  /// \pre \c domMap and \c colMap must have the same or congruent
  ///   communicators.
  /// \pre On all calling processes, the indices in \c colMap must be
  ///   a subset of the indices in \c domMap.
  ///
  /// The returned column Map's global indices (GIDs) will have the
  /// following order on all calling processes:
  ///
  /// <ol>
  /// <li> GIDs that occur in both \c colMap and \c domMap (on the
  ///      calling process) go first.
  /// </li>
  /// <li> GIDs in \c colMap on the calling process, but not in the
  ///      domain Map on the calling process, follow.  They are
  ///      ordered first contiguously by their owning process rank
  ///      (in the domain Map), then in increasing order within that.
  /// </li>
  /// </ol>
  ///
  /// This imitates the ordering used by AztecOO and Epetra.  Storing
  /// indices contiguously that are owned by the same process (in the
  /// domain Map) permits the use of contiguous send and receive
  /// buffers in Distributor, which is used in an Import operation.
  template<class MapType>
  std::pair<Teuchos::RCP<const MapType>,
            Teuchos::RCP<typename OptColMap<MapType>::import_type> >
  makeOptimizedColMapAndImport (std::ostream& errStream,
                                bool& lclErr,
                                const MapType& domMap,
                                const MapType& colMap,
                                const typename OptColMap<MapType>::import_type* oldImport = nullptr)
  {
    using local_ordinal_type = typename MapType::local_ordinal_type;
    using global_ordinal_type = typename MapType::global_ordinal_type;
    using node_type = typename MapType::node_type;
    using map_type = ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;
    using impl_type = OptColMap<map_type>;

    auto mapAndImp = impl_type::makeOptColMapAndImport (errStream, lclErr, domMap, colMap, oldImport);
    return std::make_pair (mapAndImp.first, mapAndImp.second);
  }

} // namespace Details
} // namespace Tpetra
