// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Details_makeOptimizedColMap.hpp>

// Test cases:
//
// 1. Contiguous uniform domain Map; old column Map is the same.  No
//    remote indices.  New column Map must be the same as the old
//    column Map.
//
// 2. Domain Map is in reverse order (locally, on each process) of the
//    contiguous uniform old column Map.  No remote indices.  New
//    column Map must be the same as the old column Map.
//
// 3. Contiguous uniform domain Map; column Map has its indices in
//    their original order, plus at the end, one index from the "next"
//    (wrap around) process to the right, and then one index from the
//    "next" (wrap around) process to the left.  New column Map must
//    reverse those remote indices, on all processes except Process 0.
//
// 4. Like #3, except that the domain Map indices are in reverse order
//    of their counterparts in the old column Map.
//
// 5. Same as #3, but column Map starts with the remote indices.  This
//    tests the function for the case when there are no "same" indices
//    on any process.

namespace {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;
  typedef Tpetra::global_size_t GST;
  typedef Array<int>::size_type size_type;

  // Whether the two Import objects are the same on the calling process ("locally").
  template<class ImportType>
  bool importsLocallySame (const ImportType& imp1, const ImportType& imp2)
  {
    typedef typename ImportType::map_type map_type;
    typedef typename map_type::local_ordinal_type LO;

    // 'congruent' is declared in Tpetra_Util.hpp.  It calls
    // MPI_Comm_compare if the two Comms are Teuchos::MpiComm
    // instances.  This is NOT a collective and does NOT require
    // communication.  See the MPI 3.0 standard, Section 6.4:
    // "Operations that access communicators are local and their
    // execution does not require interprocess communication."  This
    // includes MPI_Comm_compare (Section 6.4.1).
    if (! Tpetra::Details::congruent (* (imp1.getSourceMap ()->getComm ()),
                                      * (imp2.getSourceMap ()->getComm ()))) {
      return false;
    }
    else if (! Tpetra::Details::congruent (* (imp1.getTargetMap ()->getComm ()),
                                           * (imp2.getTargetMap ()->getComm ()))) {
      return false;
    }
    else if (imp1.getNumSameIDs () != imp2.getNumSameIDs ()) {
      return false;
    }
    else if (imp1.getNumPermuteIDs () != imp2.getNumPermuteIDs ()) {
      return false;
    }
    else if (imp1.getNumRemoteIDs () != imp2.getNumRemoteIDs ()) {
      return false;
    }
    else if (imp1.getNumExportIDs () != imp2.getNumExportIDs ()) {
      return false;
    }
    else {
      // Check lists of LIDs in the source Maps that are permuted.
      ArrayView<const LO> imp1permFromLids = imp1.getPermuteFromLIDs ();
      ArrayView<const LO> imp2permFromLids = imp2.getPermuteFromLIDs ();
      if (imp1permFromLids.size () != imp2permFromLids.size ()) {
        return false;
      }
      else if (! std::equal (imp1permFromLids.begin (), imp1permFromLids.end (),
                             imp2permFromLids.begin ())) {
        return false;
      }

      auto imp1_permuteFromLIDs = imp1.getPermuteFromLIDs_dv ();
      auto imp2_permuteFromLIDs = imp2.getPermuteFromLIDs_dv ();
      if (imp1_permuteFromLIDs.extent (0) != imp2_permuteFromLIDs.extent (0)) {
	return false;
      }
      else {
	auto imp1_ptr = imp1_permuteFromLIDs.view_host ().data ();
	const auto size = imp1_permuteFromLIDs.view_host ().extent (0);
	auto imp2_ptr = imp2_permuteFromLIDs.view_host ().data ();
      
	if (! std::equal (imp1_ptr, imp1_ptr + size, imp2_ptr)) {
	  return false;
	}
      }

      // Check lists of LIDs in the target Maps that are permuted.
      ArrayView<const LO> imp1permToLids = imp1.getPermuteToLIDs ();
      ArrayView<const LO> imp2permToLids = imp2.getPermuteToLIDs ();
      if (imp1permToLids.size () != imp2permToLids.size ()) {
        return false;
      }
      else if (! std::equal (imp1permToLids.begin (), imp1permToLids.end (),
                             imp2permToLids.begin ())) {
        return false;
      }

      auto imp1_permuteToLIDs = imp1.getPermuteToLIDs_dv ();
      auto imp2_permuteToLIDs = imp2.getPermuteToLIDs_dv ();
      if (imp1_permuteToLIDs.extent (0) != imp2_permuteToLIDs.extent (0)) {
	return false;
      }
      else {
	auto imp1_ptr = imp1_permuteToLIDs.view_host ().data ();
	const auto size = imp1_permuteToLIDs.view_host ().extent (0);
	auto imp2_ptr = imp2_permuteToLIDs.view_host ().data ();
      
	if (! std::equal (imp1_ptr, imp1_ptr + size, imp2_ptr)) {
	  return false;
	}
      }

      // Check LIDs in the target Maps to receive from other processes.
      ArrayView<const LO> imp1remoteLids = imp1.getRemoteLIDs ();
      ArrayView<const LO> imp2remoteLids = imp2.getRemoteLIDs ();
      if (imp1remoteLids.size () != imp2remoteLids.size ()) {
        return false;
      }
      else if (! std::equal (imp1remoteLids.begin (), imp1remoteLids.end (),
                             imp2remoteLids.begin ())) {
        return false;
      }

      auto imp1_remoteLIDs = imp1.getRemoteLIDs_dv ();
      auto imp2_remoteLIDs = imp2.getRemoteLIDs_dv ();
      if (imp1_remoteLIDs.extent (0) != imp2_remoteLIDs.extent (0)) {
	return false;
      }
      else {
	auto imp1_ptr = imp1_remoteLIDs.view_host ().data ();
	const auto size = imp1_remoteLIDs.view_host ().extent (0);
	auto imp2_ptr = imp2_remoteLIDs.view_host ().data ();
      
	if (! std::equal (imp1_ptr, imp1_ptr + size, imp2_ptr)) {
	  return false;
	}
      }
      
      // Check LIDs in the source Map that will be sent to other processes.
      ArrayView<const LO> imp1exportLids = imp1.getExportLIDs ();
      ArrayView<const LO> imp2exportLids = imp2.getExportLIDs ();
      if (imp1exportLids.size () != imp2exportLids.size ()) {
        return false;
      }
      else if (! std::equal (imp1exportLids.begin (), imp1exportLids.end (),
                             imp2exportLids.begin ())) {
        return false;
      }

      auto imp1_exportLIDs = imp1.getExportLIDs_dv ();
      auto imp2_exportLIDs = imp2.getExportLIDs_dv ();
      if (imp1_exportLIDs.extent (0) != imp2_exportLIDs.extent (0)) {
	return false;
      }
      else {
	auto imp1_ptr = imp1_exportLIDs.view_host ().data ();
	const auto size = imp1_exportLIDs.view_host ().extent (0);
	auto imp2_ptr = imp2_exportLIDs.view_host ().data ();
      
	if (! std::equal (imp1_ptr, imp1_ptr + size, imp2_ptr)) {
	  return false;
	}
      }

      // Check list of process ranks to which to send entries.
      ArrayView<const int> imp1exportPids = imp1.getExportPIDs ();
      ArrayView<const int> imp2exportPids = imp2.getExportPIDs ();
      if (imp1exportPids.size () != imp2exportPids.size ()) {
        return false;
      }
      else if (! std::equal (imp1exportPids.begin (), imp1exportPids.end (),
                             imp2exportPids.begin ())) {
        return false;
      }

      // The Imports have passed through the gauntlet of sameness tests.
      return true;
    }
  }

  // Whether the two Import objects are the same on all processes.
  template<class ImportType>
  bool importsGloballySame (const ImportType& imp1, const ImportType& imp2)
  {
    const bool lclSame = importsLocallySame<ImportType> (imp1, imp2);
    int lcl = lclSame ? 1 : 0;
    int gbl = 1;
    reduceAll<int, int> (* (imp1.getSourceMap ()->getComm ()),
                         REDUCE_MIN, lcl, outArg (gbl));
    return (gbl == 1);
  }

  // Report errors from all processes.
  void
  reportErrors (const Teuchos::Comm<int>& comm,
                const std::ostringstream& errStream,
                const bool lclErr)
  {
    const int myRank = comm.getRank ();
    const int numProcs = comm.getSize ();

    for (int p = 0; p < numProcs; ++p) {
      if (p == myRank && lclErr) {
        // FIXME (mfh 05 Jun 2014) This assumes that all processes
        // can print to cerr.  This is likely to be true on test
        // machines.  A fix would have to send error strings to Proc
        // 0 first for printing.
        cerr << "Error message from Process " << myRank << ": "
             << errStream.str () << endl;
      }
      comm.barrier (); // wait for output to finish
      comm.barrier ();
      comm.barrier ();
    }
  }

  // Mapping from a given Map specialization to its corresponding
  // Import specialization.  For example:
  //
  // typedef ::Tpetra::Map<LO, GO, NT> map_type;
  // // Get the corresponding Import type.
  // typedef typename GetImportType<map_type>::import_type import_type;
  template<class MapType>
  class GetImportType {
  public:
    using local_ordinal_type = typename MapType::local_ordinal_type;
    using global_ordinal_type = typename MapType::global_ordinal_type;
    using node_type = typename MapType::node_type;
    using import_type = ::Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>;
  };

  // For the given domain Map and (original) column Map, test
  // makeOptimizedColMap and makeOptimizedColMapAndImport.
  //
  // We provide 'out' and 'success' so that we can use macros like
  // TEST_ASSERT and TEST_EQUALITY.  Those macros use them implicitly,
  // so they must exist in the scope where the macros are used.
  template<class MapType>
  std::pair<Teuchos::RCP<const MapType>,
            Teuchos::RCP<typename GetImportType<MapType>::import_type> >
  testMakeOptColMap (Teuchos::FancyOStream& out,
                     bool& success,
                     const MapType& domMap,
                     const MapType& oldColMap)
  {
    using ::Tpetra::Details::makeOptimizedColMap;
    using ::Tpetra::Details::makeOptimizedColMapAndImport;
    typedef MapType map_type;
    // typedef typename MapType::local_ordinal_type LO;
    // typedef typename MapType::global_ordinal_type GO;
    typedef typename GetImportType<MapType>::import_type import_type;

    const Comm<int>& comm = * (domMap.getComm ());
    int lclSuccess = 1;
    int gblSuccess = 1;

    // Make sure that the domain Map and the (original) column Map
    // have congruent communicators.  If not, none of the following
    // tests will make sense.
    {
      out << "Check that the original column Map has a nonnull communicator."
          << endl;
      RCP<const Comm<int> > colMapComm = oldColMap.getComm ();
      lclSuccess = oldColMap.getComm ().is_null () ? 0 : 1;
      gblSuccess = 1;
      reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (gblSuccess == 1) {
        out << "Check that the domain Map and the original column Map "
          "have congruent communicators." << endl;
        const bool congruent = Tpetra::Details::congruent (comm, *colMapComm);
        lclSuccess = congruent ? 1 : 0;
        gblSuccess = 1;
        reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEST_EQUALITY( gblSuccess, 1 );
      }
      // Don't bother continuing the test if either of the above checks fail.
      TEUCHOS_TEST_FOR_EXCEPTION(
        gblSuccess != 1, std::logic_error, "The original column Map's "
        "communicator is either null, or is not congruent with the domain Map."
        "  In that case, it doesn't make sense for tests to continue.");
    }

    // Call makeOptimizedColMap to get a new column Map.
    out << "Calling makeOptimizedColMap" << endl;
    std::ostringstream errStrm1;
    bool lclErr1 = false;
    Teuchos::RCP<const map_type> newColMap =
      makeOptimizedColMap<map_type> (errStrm1, lclErr1, domMap, oldColMap);
    // Make sure that all processes succeeded.
    lclSuccess = lclErr1 ? 0 : 1;
    gblSuccess = 1;
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      reportErrors (comm, errStrm1, lclErr1);
    }

    // The resulting new column Map must be only a _local_ permutation
    // of the original column Map.  "Only a local permutation" means
    // the same thing as "compatible."
    out << "Test that the new column Map is only a local permutation "
      "of the original column Map." << endl;
    {
      const bool colMapsCompat = oldColMap.isCompatible (*newColMap);
      TEST_ASSERT( colMapsCompat );
    }

    // Call makeOptimizedColMapAndImport.  This
    // will make both the optimized column Map, and its corresponding
    // Import (from the domain Map 'domMap' to the new column Map).
    // Note that the function only promises to make an Import if
    // necessary, that is, if the domain Map and the new column Map
    // are not the same.
    out << "Calling makeOptimizedColMapAndImport"
        << endl;
    std::ostringstream errStrm3;
    bool lclErr3 = false;
    std::pair<RCP<const map_type>, RCP<import_type> > result3 =
      makeOptimizedColMapAndImport<map_type> (errStrm3, lclErr3, domMap,
                                              oldColMap, NULL);
    // Make sure that all processes succeeded.
    lclSuccess = lclErr3 ? 0 : 1;
    gblSuccess = 1;
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      reportErrors (comm, errStrm3, lclErr3);
    }
    // Check that calling makeOptimizedColMapAndImport 
    // produces the same column Map as makeOptimizedColMap.
    {
      const bool sameMaps = result3.first->isSameAs (*newColMap);
      TEST_ASSERT( sameMaps );
    }

    out << "Check that either the domain and new column Maps are the same, "
      "or that the returned Import is nontrivial (nonnull)." << endl;
    // Make sure that on all processes, either the domain and new
    // column Maps are the same, or the returned Import is nontrivial
    // (nonnull).
    lclSuccess = (result3.first->isSameAs (domMap) || ! result3.second.is_null ());
    gblSuccess = 1;
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );

    if (gblSuccess && ! result3.second.is_null ()) {
      // The tests below only make sense if the returned Import is
      // nontrivial (nonnull).  Checking gblSuccess as well as the
      // Import ensures that the above condition is true on all
      // processes.  That way, any calls below to collectives (like
      // Map::isSameAs or importsGloballySame) won't deadlock.

      // If the returned Import is nontrivial, make sure that it has
      // nonnull source and target Maps.  We need to test this on all
      // processes, since the isSameAs and importsGloballySame calls
      // below are collective.
      out << "The returned Import is nontrivial.  "
        "Check that its source and target Maps are nonnull." << endl;
      lclSuccess = ! result3.second->getSourceMap ().is_null () &&
        ! result3.second->getTargetMap ().is_null ();
      gblSuccess = 1;
      reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      const bool returnedSrcAndTgtMapsNonnull = (gblSuccess == 1);

      // If the returned Import is nontrivial and it has nonnull source
      // and target Maps, make sure that its source Map and domMap are
      // the same, and that its target Map and newColMap are the same.
      if (returnedSrcAndTgtMapsNonnull) {
        out << "The returned Import's source and target Maps are nonnull.  "
          "Check that its source Map is the same as domMap, and that its "
          "target Map is the same as the result of makeOptimizedColMap."
            << endl;

        // These are collectives.  We may safely call them without
        // risk of deadlock, because of the global checks above.
        TEST_ASSERT( result3.second->getSourceMap ()->isSameAs (domMap) );
        TEST_ASSERT( result3.second->getTargetMap ()->isSameAs (*newColMap) );
      }

      // Create an Import from domMap to newColMap in the usual way.
      // Make sure that it matches the Import returned by
      // makeOptimizedColMapAndImport.
      out << "Compare the returned Import against an Import created in the "
        "conventional way." << endl;
      import_type newImport (Teuchos::rcp (new map_type (domMap)),
                             result3.first);
      // It should always be the case that an Import's source and target
      // Maps are nonnull, especially if the Import was created in the
      // usual way.  It's worth checking, though.
      lclSuccess = ! newImport.getSourceMap ().is_null () &&
        ! newImport.getTargetMap ().is_null ();
      gblSuccess = 1;
      reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      const bool newSrcAndTgtMapsNonnull = (gblSuccess == 1);

      if (returnedSrcAndTgtMapsNonnull && newSrcAndTgtMapsNonnull) {
        // The three "sameness" tests below are collectives.  We may
        // safely call them without risk of deadlock, because of the
        // global checks above.

        // Check that both Imports' source Maps are the same.
        const map_type& srcMap = * (result3.second->getSourceMap ());
        const bool srcMapsSame = srcMap.isSameAs (* (newImport.getSourceMap ()));
        TEST_ASSERT( srcMapsSame );

        // Check that both Imports' target Maps are the same.
        const map_type& tgtMap = * (result3.second->getTargetMap ());
        const bool tgtMapsSame = tgtMap.isSameAs (* (newImport.getTargetMap ()));
        TEST_ASSERT( tgtMapsSame );

        // Check that both Imports are the same.
        const bool impSame =
          importsGloballySame<import_type> (* (result3.second), newImport);
        TEST_ASSERT( impSame );
      }
    }

    // The calling test may want to check the returned Map and Import
    // in other ways, so we return them.
    return result3;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MakeOptColMap, Test1, LO, GO )
  {
    typedef Tpetra::Map<LO, GO> map_type;
    typedef typename GetImportType<map_type>::import_type import_type;

    Teuchos::OSTab tab0 (out);
    out << "\"Tpetra::makeOptimizedColMap\": Test 1" << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const int myRank = comm->getRank ();
    //const int numProcs = comm->getSize ();

    const LO lclNumGids = static_cast<LO> (10);
    // const GO gblNumGids = static_cast<GO> (numProcs * lclNumGids);
    Array<GO> domMapGids (lclNumGids);
    for (LO lid = 0; lid < lclNumGids; ++lid) {
      const GO gid = static_cast<GO> (lid) + myRank * lclNumGids;
      domMapGids[lid] = gid;
    }
    ArrayView<const GO> oldColMapGids = domMapGids ();

    out << "Making the Maps" << endl;
    map_type domMap (INVALID, domMapGids (), indexBase, comm);
    map_type oldColMap (INVALID, oldColMapGids, indexBase, comm);

    // 'out' and 'success' are declared in all Teuchos unit tests.
    std::pair<RCP<const map_type>, RCP<import_type> > result =
      testMakeOptColMap<map_type> (out, success, domMap, oldColMap);

    // Specific requirement of this test.
    TEST_ASSERT( oldColMap.isSameAs (*result.first) );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MakeOptColMap, Test2, LO, GO )
  {
    typedef Tpetra::Map<LO, GO> map_type;
    typedef typename GetImportType<map_type>::import_type import_type;

    Teuchos::OSTab tab0 (out);
    out << "\"Tpetra::makeOptimizedColMap\": Test 2" << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    // const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    const LO lclNumGids = static_cast<LO> (10);
    const GO gblNumGids = static_cast<GO> (numProcs * lclNumGids);
    Array<GO> domMapGids (lclNumGids);

    for (LO k = 0; k < lclNumGids; ++k) {
      // GIDs occur in locally reverse order in this domain Map.
      const LO lid = (lclNumGids - static_cast<LO> (1)) - k;
      const GO gid = static_cast<GO> (lid) + myRank + lclNumGids;
      domMapGids[lid] = gid;
    }

    out << "Making the Maps" << endl;
    map_type domMap (static_cast<GST> (gblNumGids), domMapGids (),
                     indexBase, comm);
    map_type oldColMap (static_cast<GST> (gblNumGids),
                        static_cast<size_t> (lclNumGids),
                        indexBase, comm);

    // 'out' and 'success' are declared in all Teuchos unit tests.
    std::pair<RCP<const map_type>, RCP<import_type> > result =
      testMakeOptColMap<map_type> (out, success, domMap, oldColMap);

    // Specific requirement of this test.
    TEST_ASSERT( oldColMap.isSameAs (*result.first) );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MakeOptColMap, Test3, LO, GO )
  {
    using map_type = Tpetra::Map<LO, GO>;
    using import_type = typename GetImportType<map_type>::import_type;

    Teuchos::OSTab tab0 (out);
    out << "\"Tpetra::makeOptimizedColMap\": Test 3" << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    const LO lclNumGids = static_cast<LO> (10);
    const GO gblNumGids = static_cast<GO> (numProcs * lclNumGids);
    const GO gidOffset = static_cast<GO> (myRank) + static_cast<GO> (lclNumGids);
    Array<GO> domMapGids (lclNumGids);
    for (LO lid = 0; lid < lclNumGids; ++lid) {
      const GO gid = static_cast<GO> (lid) + myRank + lclNumGids;
      domMapGids[lid] = gid;
    }

    Array<GO> oldColMapGids (lclNumGids + 2);
    for (LO lid = 0; lid < lclNumGids; ++lid) {
      oldColMapGids[lid] = domMapGids[lid];
    }
    // GID from the process to the "right" (wrap around).
    oldColMapGids[lclNumGids] = (gidOffset + static_cast<GO> (lclNumGids)) % gblNumGids;
    if (oldColMapGids[lclNumGids] < indexBase) {
      oldColMapGids[lclNumGids] += gblNumGids;
    }
    // GID from the process to the "left" (wrap around).
    oldColMapGids[lclNumGids+1] = (gidOffset - static_cast<GO> (1)) % gblNumGids;
    if (oldColMapGids[lclNumGids+1] < indexBase) {
      oldColMapGids[lclNumGids+1] += gblNumGids;
    }

    out << "Making the Maps" << endl;
    map_type domMap (INVALID, domMapGids (), indexBase, comm);
    map_type oldColMap (INVALID, oldColMapGids, indexBase, comm);

    // 'out' and 'success' are declared in all Teuchos unit tests.
    std::pair<RCP<const map_type>, RCP<import_type> > result =
      testMakeOptColMap<map_type> (out, success, domMap, oldColMap);

    //
    // Specific requirements of this test.
    //

    TEST_ASSERT( oldColMap.isSameAs (*result.first) );
    TEST_ASSERT( ! result.second.is_null () );
    TEST_ASSERT( ! result.second.is_null () && ! result.second->getTargetMap ().is_null () );
    if (! result.second.is_null () && ! result.second->getTargetMap ().is_null ()) {
      const import_type& newImport = * (result.second);
      const map_type& newColMap = * (result.second->getTargetMap ());
      const size_t expectedNumRemoteIDs = static_cast<size_t> (2);

      TEST_ASSERT( newImport.getNumRemoteIDs () == expectedNumRemoteIDs );
      // Remote LIDs are local indices with respect to the target Map.
      // We will convert them to global indices below.
      ArrayView<const LO> remoteLids = newImport.getRemoteLIDs ();

      TEST_ASSERT( newImport.getNumRemoteIDs () == static_cast<size_t> (remoteLids.size ()) );
      if (newImport.getNumRemoteIDs () == static_cast<size_t> (remoteLids.size ()) &&
          remoteLids.size () == static_cast<size_type> (expectedNumRemoteIDs)) {
        Array<GO> remoteGids (remoteLids.size ());
        for (size_type k = 0; k < remoteLids.size (); ++k) {
          // Remote LIDs are local indices with respect to the target
          // Map.  This means we will need to convert them to global
          // indices, using the target Map, in order to determine
          // whether they are correct.
          remoteGids[k] = newColMap.getGlobalElement (remoteLids[k]);
        }

        if (myRank == 0 || myRank == numProcs - 1) {
          TEST_EQUALITY( remoteGids[0], oldColMapGids[lclNumGids+1] );
          TEST_EQUALITY( remoteGids[1], oldColMapGids[lclNumGids] );
        } else {
          TEST_EQUALITY( remoteGids[0], oldColMapGids[lclNumGids+1] );
          TEST_EQUALITY( remoteGids[1], oldColMapGids[lclNumGids] );
        }
      }
    }
  }

  //
  // INSTANTIATIONS (template tests must be instantiated in the same
  // anonymous namespace as where the tests were defined)
  //

#define UNIT_TEST_GROUP(LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MakeOptColMap, Test1, LO, GO )

  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MakeOptColMap, Test2, LO, GO )
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MakeOptColMap, Test3, LO, GO )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

} // namespace (anonymous)



