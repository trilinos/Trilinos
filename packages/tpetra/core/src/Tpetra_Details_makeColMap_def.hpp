#ifndef TPETRA_DETAILS_MAKECOLMAP_DEF_HPP
#define TPETRA_DETAILS_MAKECOLMAP_DEF_HPP

// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_makeColMap_def.hpp
/// \brief Definition of Tpetra::Details::makeColMap, a function for
///   creating the column Map of a Tpetra::CrsGraph
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra.  Users may not rely on this file or its contents.
///   They may change or disappear at any time.
///
/// This file defines the Tpetra::Details::makeColMap function, which
/// creates the column Map of a Tpetra::CrsGraph.

#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Util.hpp"
#include "Teuchos_Array.hpp"
#include "Kokkos_Bitset.hpp"
#include <set>
#include <vector>

namespace Tpetra {
namespace Details {

template <class LO, class GO, class NT>
int
makeColMapImpl(Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& colMap,
            Teuchos::Array<int>& remotePIDs,
            const Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& domMap,
            size_t numLocalColGIDs,
            size_t numRemoteColGIDs,
            std::set<GO>& RemoteGIDSet,
            std::vector<GO>& RemoteGIDUnorderedVector,
            std::vector<bool>& GIDisLocal, 
            const bool sortEachProcsGids,
            std::ostream* errStrm)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::rcp;
  using std::endl;
  int errCode = 0;
  const char prefix[] = "Tpetra::Details::makeColMapImpl: ";
  typedef ::Tpetra::Map<LO, GO, NT> map_type;
  // Possible short-circuit for serial scenario:
  //
  // If all domain GIDs are present as column indices, then set
  // ColMap=DomainMap.  By construction, LocalGIDs is a subset of
  // DomainGIDs.
  //
  // If we have
  //   * Number of remote GIDs is 0, so that ColGIDs == LocalGIDs,
  // and
  //   * Number of local GIDs is number of domain GIDs
  // then
  //   * LocalGIDs \subset DomainGIDs && size(LocalGIDs) ==
  //     size(DomainGIDs) => DomainGIDs == LocalGIDs == ColGIDs
  // on the calling process.
  //
  // We will concern ourselves only with the special case of a
  // serial DomainMap, obviating the need for communication.
  //
  // If
  //   * DomainMap has a serial communicator
  // then we can set the column Map as the domain Map
  // return. Benefit: this graph won't need an Import object
  // later.
  //
  // Note, for a serial domain map, there can be no RemoteGIDs,
  // because there are no remote processes.  Likely explanations
  // for this are:
  //  * user submitted erroneous column indices
  //  * user submitted erroneous domain Map
  if (domMap->getComm ()->getSize () == 1) {
    if (numRemoteColGIDs != 0) {
      errCode = -2;
      if (errStrm != NULL) {
        *errStrm << prefix << "The domain Map only has one process, but "
                 << numRemoteColGIDs << " column "
                 << (numRemoteColGIDs != 1 ? "indices are" : "index is")
                 << " not in the domain Map.  Either these indices are "
          "invalid or the domain Map is invalid.  Remember that nonsquare "
          "matrices, or matrices where the row and range Maps differ, "
          "require calling the version of fillComplete that takes the "
          "domain and range Maps as input." << endl;
      }
    }
    if (numLocalColGIDs == domMap->getLocalNumElements ()) {
      colMap = domMap; // shallow copy
      return errCode;
    }
  }

  // Populate myColumns with a list of all column GIDs.  Put
  // locally owned (in the domain Map) GIDs at the front: they
  // correspond to "same" and "permuted" entries between the
  // column Map and the domain Map.  Put remote GIDs at the back.
  Array<GO> myColumns(numLocalColGIDs + numRemoteColGIDs);
  // get pointers into myColumns for each part
  ArrayView<GO> LocalColGIDs  = myColumns (0, numLocalColGIDs);
  ArrayView<GO> remoteColGIDs = myColumns (numLocalColGIDs, numRemoteColGIDs);

  // Copy the remote GIDs into myColumns
  if (sortEachProcsGids) {
    // The std::set puts GIDs in increasing order.
    std::copy (RemoteGIDSet.begin(), RemoteGIDSet.end(),
               remoteColGIDs.begin());
  }
  else {
    // Respect the originally encountered order.
    std::copy (RemoteGIDUnorderedVector.begin(),
               RemoteGIDUnorderedVector.end(), remoteColGIDs.begin());
  }

  // Make a list of process ranks corresponding to the remote GIDs.
  // remotePIDs is an output argument of getRemoteIndexList below;
  // its initial contents don't matter.
  if (static_cast<size_t> (remotePIDs.size ()) != numRemoteColGIDs) {
    remotePIDs.resize (numRemoteColGIDs);
  }
  // Look up the remote process' ranks in the domain Map.
  {
    const LookupStatus stat =
      domMap->getRemoteIndexList (remoteColGIDs, remotePIDs ());

    // If any process returns IDNotPresent, then at least one of
    // the remote indices was not present in the domain Map.  This
    // means that the Import object cannot be constructed, because
    // of incongruity between the column Map and domain Map.
    // This has two likely causes:
    //   - The user has made a mistake in the column indices
    //   - The user has made a mistake with respect to the domain Map
    if (stat == IDNotPresent) {
      if (errStrm != NULL) {
        *errStrm << prefix << "Some column indices are not in the domain Map."
          "Either these column indices are invalid or the domain Map is "
          "invalid.  Likely cause: For a nonsquare matrix, you must give the "
          "domain and range Maps as input to fillComplete." << endl;
      }
      // Don't return yet, because not all processes may have
      // encountered this error state.  This function ends with an
      // all-reduce, so we have to make sure that everybody gets to
      // that point.  The resulting Map may be wrong, but at least
      // nothing should crash.
      errCode = -3;
    }
  }

  // Sort incoming remote column indices by their owning process
  // rank, so that all columns coming from a given remote process
  // are contiguous.  This means the Import's Distributor doesn't
  // need to reorder data.
  //
  // NOTE (mfh 02 Sep 2014) This needs to be a stable sort, so that
  // it respects either of the possible orderings of GIDs (sorted,
  // or original order) specified above.
  sort2 (remotePIDs.begin (), remotePIDs.end (), remoteColGIDs.begin (), true);

  // Copy the local GIDs into myColumns. Two cases:
  // 1. If the number of Local column GIDs is the same as the number
  //    of Local domain GIDs, we can simply read the domain GIDs
  //    into the front part of ColIndices (see logic above from the
  //    serial short circuit case)
  // 2. We step through the GIDs of the DomainMap, checking to see
  //    if each domain GID is a column GID.  We want to do this to
  //    maintain a consistent ordering of GIDs between the columns
  //    and the domain.

  const size_t numDomainElts = domMap->getLocalNumElements ();
  if (numLocalColGIDs == numDomainElts) {
    // If the number of locally owned GIDs are the same as the
    // number of local domain Map elements, then the local domain
    // Map elements are the same as the locally owned GIDs.
    if (domMap->isContiguous ()) {
      // NOTE (mfh 03 Mar 2013, 02 Sep 2014) In the common case that
      // the domain Map is contiguous, it's more efficient to avoid
      // calling getLocalElementList(), since that permanently
      // constructs and caches the GID list in the contiguous Map.
      GO curColMapGid = domMap->getMinGlobalIndex ();
      for (size_t k = 0; k < numLocalColGIDs; ++k, ++curColMapGid) {
        LocalColGIDs[k] = curColMapGid;
      }
    }
    else {
      ArrayView<const GO> domainElts = domMap->getLocalElementList ();
      std::copy (domainElts.begin(), domainElts.end(), LocalColGIDs.begin());
    }
  }
  else {
    // Count the number of locally owned GIDs, both to keep track of
    // the current array index, and as a sanity check.
    size_t numLocalCount = 0;
    if (domMap->isContiguous ()) {
      // NOTE (mfh 03 Mar 2013, 02 Sep 2014) In the common case that
      // the domain Map is contiguous, it's more efficient to avoid
      // calling getLocalElementList(), since that permanently
      // constructs and caches the GID list in the contiguous Map.
      GO curColMapGid = domMap->getMinGlobalIndex ();
      for (size_t i = 0; i < numDomainElts; ++i, ++curColMapGid) {
        if (GIDisLocal[i]) {
          LocalColGIDs[numLocalCount++] = curColMapGid;
        }
      }
    }
    else {
      ArrayView<const GO> domainElts = domMap->getLocalElementList ();
      for (size_t i = 0; i < numDomainElts; ++i) {
        if (GIDisLocal[i]) {
          LocalColGIDs[numLocalCount++] = domainElts[i];
        }
      }
    }
    if (numLocalCount != numLocalColGIDs) {
      if (errStrm != NULL) {
        *errStrm << prefix << "numLocalCount = " << numLocalCount
                 << " != numLocalColGIDs = " << numLocalColGIDs
                 << ".  This should never happen.  "
          "Please report this bug to the Tpetra developers." << endl;
      }
      // Don't return yet, because not all processes may have
      // encountered this error state.  This function ends with an
      // all-reduce, so we have to make sure that everybody gets to
      // that point.
      errCode = -4;
    }
  }

  // FIXME (mfh 03 Apr 2013) Now would be a good time to use the
  // information we collected above to construct the Import.  In
  // particular, building an Import requires:
  //
  // 1. numSameIDs (length of initial contiguous sequence of GIDs
  //    on this process that are the same in both Maps; this
  //    equals the number of domain Map elements on this process)
  //
  // 2. permuteToLIDs and permuteFromLIDs (both empty in this
  //    case, since there's no permutation going on; the column
  //    Map starts with the domain Map's GIDs, and immediately
  //    after them come the remote GIDs)
  //
  // 3. remoteGIDs (exactly those GIDs that we found out above
  //    were not in the domain Map) and remoteLIDs (which we could
  //    have gotten above by using the three-argument version of
  //    getRemoteIndexList() that computes local indices as well
  //    as process ranks, instead of the two-argument version that
  //    was used above)
  //
  // 4. remotePIDs (which we have from the getRemoteIndexList() call
  //    above) -- NOTE (mfh 14 Sep 2017) Fix for Trilinos GitHub
  //    Issue #628 (https://github.com/trilinos/Trilinos/issues/628)
  //    addresses this.  remotePIDs is now an output argument of
  //    both this function and CrsGraph::makeColMap, and
  //    CrsGraph::makeImportExport now has an option to use that
  //    information, if available (that is, if makeColMap was
  //    actually called, which it would not be if the graph already
  //    had a column Map).
  //
  // 5. Apply the permutation from sorting remotePIDs to both
  //    remoteGIDs and remoteLIDs (by calling sort3 above instead of
  //    sort2), instead of just to remoteLIDs alone.
  //
  // 6. Everything after the sort3 call in Import::setupExport():
  //    a. Create the Distributor via createFromRecvs(), which
  //       computes exportGIDs and exportPIDs
  //    b. Compute exportLIDs from exportGIDs (by asking the
  //       source Map, in this case the domain Map, to convert
  //       global to local)
  //
  // Steps 1-5 come for free, since we must do that work anyway in
  // order to compute the column Map.  In particular, Step 3 is
  // even more expensive than Step 6a, since it involves both
  // creating and using a new Distributor object.
  const global_size_t INV =
    Tpetra::Details::OrdinalTraits<global_size_t>::invalid ();
  // FIXME (mfh 05 Mar 2014) Doesn't the index base of a Map have to
  // be the same as the Map's min GID? If the first column is empty
  // (contains no entries), then the column Map's min GID won't
  // necessarily be the same as the domain Map's index base.
  const GO indexBase = domMap->getIndexBase ();
  colMap = rcp (new map_type (INV, myColumns, indexBase, domMap->getComm ()));
  return errCode;
}

template <class LO, class GO, class NT>
int
makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& colMap,
            Teuchos::Array<int>& remotePIDs,
            const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& domMap,
            const RowGraph<LO, GO, NT>& graph,
            const bool sortEachProcsGids,
            std::ostream* errStrm)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::rcp;
  using std::endl;
  typedef ::Tpetra::Map<LO, GO, NT> map_type;
  const char prefix[] = "Tpetra::Details::makeColMap: ";
  int errCode = 0;

  // If the input domain Map or its communicator is null on the
  // calling process, then the calling process does not participate in
  // the returned column Map.  Thus, we can set the returned column
  // Map to null on those processes, and return immediately.  This is
  // not an error condition, as long as when domMap and its
  // communicator are NOT null, the graph's other Maps and
  // communicator are not also null.
  if (domMap.is_null () || domMap->getComm ().is_null ()) {
    colMap = Teuchos::null;
    return errCode;
  }

  Array<GO> myColumns;
  if (graph.isLocallyIndexed ()) {
    colMap = graph.getColMap ();
    // If the graph is locally indexed, it had better have a column Map.
    // The extra check for ! graph.hasColMap() is conservative.
    if (colMap.is_null () || ! graph.hasColMap ()) {
      errCode = -1;
      if (errStrm != NULL) {
        *errStrm << prefix << "The graph is locally indexed on the calling "
          "process, but has no column Map (either getColMap() returns null, "
          "or hasColMap() returns false)." << endl;
      }
      // Under this error condition, this process will not fill
      // myColumns.  The resulting new column Map will be incorrect,
      // but at least we won't hang, and this process will report the
      // error.
    }
    else {
      // The graph already has a column Map, and is locally indexed on
      // the calling process.  However, it may be globally indexed (or
      // neither locally nor globally indexed) on other processes.
      // Assume that we want to recreate the column Map.
      if (colMap->isContiguous ()) {
        // The number of indices on each process must fit in LO.
        const LO numCurGids = static_cast<LO> (colMap->getLocalNumElements ());
        myColumns.resize (numCurGids);
        const GO myFirstGblInd = colMap->getMinGlobalIndex ();
        for (LO k = 0; k < numCurGids; ++k) {
          myColumns[k] = myFirstGblInd + static_cast<GO> (k);
        }
      }
      else { // the column Map is NOT contiguous
        ArrayView<const GO> curGids = graph.getColMap ()->getLocalElementList ();
        // The number of indices on each process must fit in LO.
        const LO numCurGids = static_cast<LO> (curGids.size ());
        myColumns.resize (numCurGids);
        for (LO k = 0; k < numCurGids; ++k) {
          myColumns[k] = curGids[k];
        }
      } // whether the graph's current column Map is contiguous
    } // does the graph currently have a column Map?
  }
  else if (graph.isGloballyIndexed ()) {
    // Go through all the rows, finding the populated column indices.
    //
    // Our final list of indices for the column Map constructor will
    // have the following properties (all of which are with respect to
    // the calling process):
    //
    // 1. Indices in the domain Map go first.
    // 2. Indices not in the domain Map follow, ordered first
    //    contiguously by their owning process rank (in the domain
    //    Map), then in increasing order within that.
    // 3. No duplicate indices.
    //
    // This imitates the ordering used by Aztec(OO) and Epetra.
    // Storing indices owned by the same process (in the domain Map)
    // contiguously permits the use of contiguous send and receive
    // buffers.
    //
    // We begin by partitioning the column indices into "local" GIDs
    // (owned by the domain Map) and "remote" GIDs (not owned by the
    // domain Map).  We use the same order for local GIDs as the
    // domain Map, so we track them in place in their array.  We use
    // an std::set (RemoteGIDSet) to keep track of remote GIDs, so
    // that we don't have to merge duplicates later.
    const LO LINV = Tpetra::Details::OrdinalTraits<LO>::invalid ();
    size_t numLocalColGIDs = 0;
    size_t numRemoteColGIDs = 0;

    // GIDisLocal[lid] is false if and only if local index lid in the
    // domain Map is remote (not local).
    std::vector<bool> GIDisLocal (domMap->getLocalNumElements (), false);
    std::set<GO> RemoteGIDSet;
    // This preserves the not-sorted Epetra order of GIDs.
    // We only use this if sortEachProcsGids is false.
    std::vector<GO> RemoteGIDUnorderedVector;

    if (! graph.getRowMap ().is_null ()) {
      const Tpetra::Map<LO, GO, NT>& rowMap = * (graph.getRowMap ());
      const LO lclNumRows = rowMap.getLocalNumElements ();
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const GO gblRow = rowMap.getGlobalElement (lclRow);
        typename RowGraph<LO,GO,NT>::global_inds_host_view_type rowGids;
        graph.getGlobalRowView (gblRow, rowGids);

        const LO numEnt = static_cast<LO> (rowGids.size ());
        if (numEnt != 0) {
          for (LO k = 0; k < numEnt; ++k) {
            const GO gid = rowGids[k];
            const LO lid = domMap->getLocalElement (gid);
            if (lid != LINV) {
              const bool alreadyFound = GIDisLocal[lid];
              if (! alreadyFound) {
                GIDisLocal[lid] = true;
                ++numLocalColGIDs;
              }
            }
            else {
              const bool notAlreadyFound = RemoteGIDSet.insert (gid).second;
              if (notAlreadyFound) { // gid did not exist in the set before
                if (! sortEachProcsGids) {
                  // The user doesn't want to sort remote GIDs (for each
                  // remote process); they want us to keep remote GIDs
                  // in their original order.  We do this by stuffing
                  // each remote GID into an array as we encounter it
                  // for the first time.  The std::set helpfully tracks
                  // first encounters.
                  RemoteGIDUnorderedVector.push_back (gid);
                }
                ++numRemoteColGIDs;
              }
            }
          } // for each entry k in row r
        } // if row r contains > 0 entries
      } // for each locally owned row r
    } // if the graph has a nonnull row Map

    return makeColMapImpl<LO, GO, NT>(
        colMap, remotePIDs,
        domMap,
        numLocalColGIDs, numRemoteColGIDs,
        RemoteGIDSet, RemoteGIDUnorderedVector, GIDisLocal,
        sortEachProcsGids, errStrm);

  } // if the graph is globally indexed
  else {
    // If we reach this point, the graph is neither locally nor
    // globally indexed.  Thus, the graph is empty on this process
    // (per the usual legacy Petra convention), so myColumns will be
    // left empty.
    ; // do nothing
  }

  const global_size_t INV =
    Tpetra::Details::OrdinalTraits<global_size_t>::invalid ();
  // FIXME (mfh 05 Mar 2014) Doesn't the index base of a Map have to
  // be the same as the Map's min GID? If the first column is empty
  // (contains no entries), then the column Map's min GID won't
  // necessarily be the same as the domain Map's index base.
  const GO indexBase = domMap->getIndexBase ();
  colMap = rcp (new map_type (INV, myColumns, indexBase, domMap->getComm ()));
  return errCode;
}

template<typename GOView, typename bitset_t>
struct GatherPresentEntries
{
  using GO = typename GOView::non_const_value_type;

  GatherPresentEntries(GO minGID_, const GOView& gids_, const bitset_t& present_)
    : minGID(minGID_), gids(gids_), present(present_)
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const GO i) const
  {
    present.set(gids(i) - minGID);
  }

  GO minGID;
  GOView gids;
  bitset_t present;
};

template<typename LO, typename GO, typename device_t, typename LocalMapType, typename const_bitset_t, bool doingRemotes>
struct ListGIDs
{
  using mem_space = typename device_t::memory_space;
  using GOView = Kokkos::View<GO*, mem_space>;
  using SingleView = Kokkos::View<GO, mem_space>;

  ListGIDs(GO minGID_, GOView& gidList_, SingleView& numElems_, const_bitset_t& present_, const LocalMapType& localDomainMap_)
    : minGID(minGID_), gidList(gidList_), numElems(numElems_), present(present_), localDomainMap(localDomainMap_)
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const GO i, GO& lcount, const bool finalPass) const
  {
    bool isRemote = localDomainMap.getLocalElement(i + minGID) == ::Tpetra::Details::OrdinalTraits<LO>::invalid();
    if(present.test(i) && doingRemotes == isRemote)
    {
      if(finalPass)
      {
        //lcount is the index where this GID should be inserted in gidList.
        gidList(lcount) = minGID + i;
      }
      lcount++;
    }
    if((i == static_cast<GO>(present.size() - 1)) && finalPass)
    {
      //Set the number of inserted indices in a single-element view
      numElems() = lcount;
    }
  }

  GO minGID;
  GOView gidList;
  SingleView numElems;
  const_bitset_t present;
  const LocalMapType localDomainMap;
};

template<typename GO, typename mem_space>
struct MinMaxReduceFunctor
{
  using MinMaxValue = typename Kokkos::MinMax<GO>::value_type;
  using GOView = Kokkos::View<GO*, mem_space>;

  MinMaxReduceFunctor(const GOView& gids_)
    : gids(gids_)
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const GO i, MinMaxValue& lminmax) const
  {
    GO gid = gids(i);
    if(gid < lminmax.min_val)
      lminmax.min_val = gid;
    if(gid > lminmax.max_val)
      lminmax.max_val = gid;
  };

  const GOView gids;
};

template <class LO, class GO, class NT>
int
makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& colMap,
            const Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& domMap,
            Kokkos::View<GO*, typename NT::memory_space> gids,
            std::ostream* errStrm)
{
  using Teuchos::RCP;
  using Teuchos::Array;
  using Kokkos::RangePolicy;
  using device_t = typename NT::device_type;
  using exec_space = typename device_t::execution_space;
  using memory_space = typename device_t::memory_space;
  // Note BMK 5-2021: this is deliberately not just device_t.
  // Bitset cannot use HIPHostPinnedSpace currently, so this needs to
  // use the default memory space for HIP (HIPSpace). Using the default mem
  // space is fine for all other backends too. This bitset type is only used
  // in this function so it won't cause type mismatches.
  using bitset_t = Kokkos::Bitset<typename exec_space::memory_space>;
  using const_bitset_t = Kokkos::ConstBitset<typename exec_space::memory_space>;
  using GOView = Kokkos::View<GO*, memory_space>;
  using SingleView = Kokkos::View<GO, memory_space>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using LocalMap = typename map_type::local_map_type;
  GO nentries = gids.extent(0);
  GO minGID = Teuchos::OrdinalTraits<GO>::max();
  GO maxGID = 0;
  using MinMaxValue = typename Kokkos::MinMax<GO>::value_type;
  MinMaxValue minMaxGID = {minGID, maxGID};
  Kokkos::parallel_reduce(RangePolicy<exec_space>(0, nentries),
      MinMaxReduceFunctor<GO, memory_space>(gids),
      Kokkos::MinMax<GO>(minMaxGID));
  minGID = minMaxGID.min_val; maxGID = minMaxGID.max_val;
  //Now, know the full range of input GIDs.
  //Determine the set of GIDs in the column map using a dense bitset, which corresponds to the range [minGID, maxGID]
  bitset_t presentGIDs(maxGID - minGID + 1);
  Kokkos::parallel_for(RangePolicy<exec_space>(0, nentries), GatherPresentEntries<GOView, bitset_t>(minGID, gids, presentGIDs));
  const_bitset_t constPresentGIDs(presentGIDs);
  //Get the set of local and remote GIDs on device
  SingleView numLocals("Num local GIDs");
  SingleView numRemotes("Num remote GIDs");
  GOView localGIDView(Kokkos::ViewAllocateWithoutInitializing("Local GIDs"), constPresentGIDs.count());
  GOView remoteGIDView(Kokkos::ViewAllocateWithoutInitializing("Remote GIDs"), constPresentGIDs.count());
  LocalMap localDomMap = domMap->getLocalMap();
  //This lists the locally owned GIDs in localGIDView
  Kokkos::parallel_scan(RangePolicy<exec_space>(0, constPresentGIDs.size()),
      ListGIDs<LO, GO, device_t, LocalMap, const_bitset_t, false>
      (minGID, localGIDView, numLocals, constPresentGIDs, localDomMap));
  //And this lists the remote GIDs in remoteGIDView
  Kokkos::parallel_scan(RangePolicy<exec_space>(0, constPresentGIDs.size()),
      ListGIDs<LO, GO, device_t, LocalMap, const_bitset_t, true>
      (minGID, remoteGIDView, numRemotes, constPresentGIDs, localDomMap));
  //Pull down the sizes
  GO numLocalColGIDs = 0;
  GO numRemoteColGIDs = 0;
  // DEEP_COPY REVIEW - DEVICE-TO-VALUE
  Kokkos::deep_copy(exec_space(), numLocalColGIDs, numLocals);
  // DEEP_COPY REVIEW - DEVICE-TO-numLocalColGIDs
  Kokkos::deep_copy(exec_space(), numRemoteColGIDs, numRemotes);
  //Pull down the remote lists
  auto localsHost = Kokkos::create_mirror_view(localGIDView);
  auto remotesHost = Kokkos::create_mirror_view(remoteGIDView);
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  Kokkos::deep_copy(exec_space(), localsHost, localGIDView);
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  Kokkos::deep_copy(exec_space(), remotesHost, remoteGIDView);
  // CAG: This fence was found to be required on Cuda with UVM=on.
  Kokkos::fence("Tpetra::makeColMap");
  //Finally, populate the STL structures which hold the index lists
  std::set<GO> RemoteGIDSet;
  std::vector<GO> RemoteGIDUnorderedVector;
  std::vector<bool> GIDisLocal (domMap->getLocalNumElements (), false);
  for(GO i = 0; i < numLocalColGIDs; i++)
  {
    GO gid = localsHost(i);
    //Already know that gid is locally owned, so this index will never be invalid().
    //makeColMapImpl uses this and the domain map to get the the local GID list.
    GIDisLocal[domMap->getLocalElement(gid)] = true;
  }
  RemoteGIDUnorderedVector.reserve(numRemoteColGIDs);
  for(GO i = 0; i < numRemoteColGIDs; i++)
  {
    GO gid = remotesHost(i);
    RemoteGIDSet.insert(gid);
    RemoteGIDUnorderedVector.push_back(gid);
  }
  //remotePIDs will be discarded in this version of makeColMap
  Array<int> remotePIDs;
  //Find the min and max GID
  return makeColMapImpl<LO, GO, NT>(
            colMap,
            remotePIDs,
            domMap,
            static_cast<size_t>(numLocalColGIDs),
            static_cast<size_t>(numRemoteColGIDs),
            RemoteGIDSet,
            RemoteGIDUnorderedVector,
            GIDisLocal, 
            true, //always sort remotes
            errStrm);
}

} // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_DETAILS_MAKECOLMAP_INSTANT(LO,GO,NT) \
  namespace Details { \
    template int \
    makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >&, \
                Teuchos::Array<int>&, \
                const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >&, \
                const RowGraph<LO, GO, NT>&, \
                const bool, \
                std::ostream*); \
    template int \
    makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >&, \
                const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >&, \
                Kokkos::View<GO*, typename NT::memory_space>, \
                std::ostream*); \
  }

#endif // TPETRA_DETAILS_MAKECOLMAP_DEF_HPP
