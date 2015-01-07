/*
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
*/

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

#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Util.hpp>

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
    typedef MapType map_type;
    typedef typename MapType::local_ordinal_type local_ordinal_type;
    typedef typename MapType::global_ordinal_type global_ordinal_type;
    typedef typename MapType::node_type node_type;
    typedef Import<local_ordinal_type,
                   global_ordinal_type,
                   node_type> import_type;

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
    ///   new Import object.  This function will <i>only</i> look at
    ///   this pointer if \c makeImport is true.
    /// \param makeImport [in] Whether to make and return an Import from
    ///   the input domain Map to the new column Map.
    ///
    /// \return The possibly reordered column Map \c newColMap, and the
    ///   corresponding Import from \c domMap to \c newColMap.  The
    ///   latter is nonnull if and only if \c makeImport is true.
    ///
    /// \pre \c domMap and \c colMap must have the same or congruent
    ///   communicators.
    /// \pre On all calling processes, the indices in \c colMap must be
    ///   a subset of the indices in \c domMap.
    static std::pair<map_type, Teuchos::RCP<import_type> >
    make (std::ostream& errStream,
          bool& lclErr,
          const map_type& domMap,
          const map_type& colMap,
          const import_type* oldImport,
          const bool makeImport)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using std::endl;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      const char prefix[] = "Tpetra::makeOptimizedColMapAndImport: ";
      std::ostream& err = errStream;

      (void) oldImport; // We don't currently use this argument.

      RCP<const Teuchos::Comm<int> > comm = colMap.getComm ();
      const LO colMapMinLid = colMap.getMinLocalIndex ();
      const LO colMapMaxLid = colMap.getMaxLocalIndex ();

      // Count the numbers of GIDs in colMap that are in and not in
      // domMap on the calling process.  Check for zero indices on the
      // calling process first, because if it's true, then we shouldn't
      // trust [getMinLocalIndex(), getMaxLocalIndex()] to return a
      // correct range.
      LO numOwnedGids = 0;
      LO numRemoteGids = 0;
      if (colMap.getNodeNumElements () != 0) {
        for (LO colMapLid = colMapMinLid; colMapLid <= colMapMaxLid; ++colMapLid) {
          const GO colMapGid = colMap.getGlobalElement (colMapLid);
          if (domMap.isNodeLocalElement (colMapGid)) {
            ++numOwnedGids;
          } else {
            ++numRemoteGids;
          }
        }
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
      if (colMap.getNodeNumElements () != 0) {
        for (LO colMapLid = colMapMinLid; colMapLid <= colMapMaxLid; ++colMapLid) {
          const GO colMapGid = colMap.getGlobalElement (colMapLid);
          if (domMap.isNodeLocalElement (colMapGid)) {
            ownedGids[ownedPos++] = colMapGid;
          } else {
            remoteGids[remotePos++] = colMapGid;
          }
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
      Array<LO> remoteLids;
      if (makeImport) {
        remoteLids.resize (numRemoteGids);
        std::fill (remoteLids.begin (), remoteLids.end (),
                   Teuchos::OrdinalTraits<LO>::invalid ());
      }
      LookupStatus lookupStatus;
      if (makeImport) {
        lookupStatus = domMap.getRemoteIndexList (remoteGids, remotePids (),
                                                  remoteLids ());
      } else {
        lookupStatus = domMap.getRemoteIndexList (remoteGids, remotePids ());
      }

      // If any process returns IDNotPresent, then at least one of the
      // remote indices was not present in the domain Map.  This means
      // that the Import object cannot be constructed, because of
      // incongruity between the column Map and domain Map.  This means
      // that either the column Map or domain Map, or both, is
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

      // Sort incoming remote column Map indices so that all columns
      // coming from a given remote process are contiguous.  This means
      // the Import's Distributor doesn't need to reorder data.
      if (makeImport) {
        sort2 (remotePids.begin (), remotePids.end (), remoteGids.begin ());
      }
      else {
        sort3 (remotePids.begin (), remotePids.end (),
               remoteGids.begin (),
               remoteLids.begin ());
      }
      // Make the new column Map.
      MapType newColMap (colMap.getGlobalNumElements (), allGids (),
                         colMap.getIndexBase (), comm, colMap.getNode ());
      // Optionally, make the new Import object.
      RCP<import_type> imp;
      if (makeImport) {
        imp = rcp (new import_type (rcp (new map_type (domMap)),
                                    rcp (new map_type (newColMap))));
        // FIXME (mfh 06 Jul 2014) This constructor throws a runtime
        // error, so I'm not using it for now.
        //
        // imp = rcp (new import_type (domMap, newColMap, remoteGids,
        //                             remotePids (), remoteLids (),
        //                             Teuchos::null, Teuchos::null));
      }
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
  ///
  /// \return The possibly reordered column Map \c newColMap.
  ///
  /// This is a convenience wrapper for makeOptimizedColMapAndImport().
  /// (Please refer to that function's documentation in this file.)
  /// It does everything that that function does, except that it does
  /// not compute a new Import.
  template<class MapType>
  MapType
  makeOptimizedColMap (std::ostream& errStream,
                       bool& lclErr,
                       const MapType& domMap,
                       const MapType& colMap)
  {
    typedef typename MapType::local_ordinal_type LO;
    typedef typename MapType::global_ordinal_type GO;
    typedef typename MapType::node_type NT;
    typedef ::Tpetra::Import<LO, GO, NT> import_type;

    const bool makeImport = false;
    std::pair<MapType, Teuchos::RCP<import_type> > ret =
      OptColMap<MapType>::make (errStream, lclErr, domMap, colMap,
                                NULL, makeImport);
    return ret.first;
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
  /// \param makeImport [in] Whether to make and return an Import from
  ///   the input domain Map to the new column Map.
  ///
  /// \return The possibly reordered column Map \c newColMap, and the
  ///   corresponding Import from \c domMap to \c newColMap.  The
  ///   latter is nonnull if and only if \c makeImport is true.
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
  std::pair<MapType, Teuchos::RCP<typename OptColMap<MapType>::import_type> >
  makeOptimizedColMapAndImport (std::ostream& errStream,
                                bool& lclErr,
                                const MapType& domMap,
                                const MapType& colMap,
                                const typename OptColMap<MapType>::import_type* oldImport,
                                const bool makeImport)
  {
    return OptColMap<MapType>::make (errStream, lclErr, domMap, colMap,
                                     oldImport, makeImport);
  }

} // namespace Details
} // namespace Tpetra
