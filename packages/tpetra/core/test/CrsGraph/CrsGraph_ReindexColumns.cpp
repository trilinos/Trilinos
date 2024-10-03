// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsGraph.hpp"

namespace {
  //using Teuchos::broadcast;
  //using std::endl;

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;
  typedef Tpetra::global_size_t GST;
  typedef Teuchos::Array<size_t>::size_type size_type;

  //
  // UNIT TESTS
  //

  //
  // Test CrsGraph::reindexColumns, with only the column Map (no
  // Import) and with sorting on.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphReindexColumns, ColMapOnlySortingOn, LO, GO, Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    using lids_type = typename graph_type::nonconst_local_inds_host_view_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    int gblSuccess = 0;
    int lclSuccess = 1;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Collect all the debug output on each process here, so we can
    // print it coherently across processes at the end.
    std::ostringstream os;

    // Create the graph's row Map.
    // const size_t numLocalIndices = 5;
    const GST numGlobalIndices = static_cast<GST> (5 * numProcs);
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (numGlobalIndices, indexBase, comm,
                         Tpetra::GloballyDistributed));
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! rowMap->isContiguous (), std::logic_error, "The row Map is supposed "
      "to be contiguous, but is not.");

    const size_t maxNumEntPerRow = 3;
    graph_type graph (rowMap, maxNumEntPerRow);

    // Make the usual tridiagonal graph.  Let the graph create its own
    // column Map.  We'll use that column Map to create a new column
    // Map, and give that column Map to reindexColumns().

    if (rowMap->getLocalNumElements () > 0) {
      const GO myMinGblInd = rowMap->getMinGlobalIndex ();
      const GO myMaxGblInd = rowMap->getMaxGlobalIndex ();
      const GO gblMinGblInd = rowMap->getMinAllGlobalIndex ();
      const GO gblMaxGblInd = rowMap->getMaxAllGlobalIndex ();

      Array<GO> gblInds (maxNumEntPerRow);
      for (GO gblRowInd = myMinGblInd; gblRowInd <= myMaxGblInd; ++gblRowInd) {
        size_t numInds = 0;
        if (gblRowInd == gblMinGblInd) {
          if (gblRowInd < gblMaxGblInd) {
            numInds = 2;
            gblInds[0] = gblRowInd;
            gblInds[1] = gblRowInd + 1;
          } else { // special case of 1 x 1 graph
            numInds = 1;
            gblInds[0] = gblRowInd;
          }
        } else if (gblRowInd == gblMaxGblInd) {
          if (gblRowInd > gblMinGblInd) {
            numInds = 2;
            gblInds[0] = gblRowInd - 1;
            gblInds[1] = gblRowInd;
          } else { // special case of 1 x 1 graph
            numInds = 1;
            gblInds[0] = gblRowInd;
          }
        } else {
          numInds = 3;
          gblInds[0] = gblRowInd - 1;
          gblInds[1] = gblRowInd;
          gblInds[2] = gblRowInd + 1;
        }
        ArrayView<const GO> gblIndsView = gblInds (0, numInds);
        graph.insertGlobalIndices (gblRowInd, gblIndsView);
      }
    }
    graph.fillComplete ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isFillComplete (), std::logic_error, "The graph claims that it "
      "is not fill complete, after fillComplete was called.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.hasColMap () || graph.getColMap ().is_null (), std::logic_error,
      "The graph is fill complete, but doesn't have a column Map.");

    // Make a deep copy of the graph, to check later that the
    // conversion was correct.
    RCP<graph_type> graph2 = rcp(new graph_type(graph));
    TEST_ASSERT( ! graph2.is_null () );

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Create a new column Map, which has all the global indices on
    // each process (_locally_) in reverse order of the graph's
    // current column Map.

    // NOTE (mfh 21 Aug 2014) Don't get this as a reference!  Get it
    // as an RCP!  Remember that the reference comes from the result
    // of getColMap(), and the call to reindexColumns() will
    // invalidate the graph's current column Map.
    RCP<const map_type> curColMap = graph.getColMap ();
    Array<GO> newGblInds (curColMap->getLocalNumElements ());
    if (curColMap->isContiguous ()) {
      // const GO myMinGblInd = curColMap->getMinGlobalIndex ();
      const GO myMaxGblInd = curColMap->getMaxGlobalIndex ();

      const size_type myNumInds =
        static_cast<size_type> (curColMap->getLocalNumElements ());
      if (myNumInds > 0) {
        GO curGblInd = myMaxGblInd;
        for (size_type k = 0; k < myNumInds; ++k, --curGblInd) {
          newGblInds[k] = curGblInd;
        }
      }
    } else { // original column Map is not contiguous
      ArrayView<const GO> curGblInds = curColMap->getLocalElementList ();
      for (size_type k = 0; k < curGblInds.size (); ++k) {
        const size_type k_opposite = (newGblInds.size () - 1) - k;
        newGblInds[k] = curGblInds[k_opposite];
      }
    }

    RCP<const map_type> newColMap =
      rcp (new map_type (INVALID, newGblInds (), indexBase, comm));

    // Print both old and new column Maps.
    {
      comm->barrier ();
      RCP<Teuchos::FancyOStream> errStream =
        Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));

      if (myRank == 0) {
        cerr << "Original column Map:" << endl;
      }
      curColMap->describe (*errStream, Teuchos::VERB_EXTREME);

      if (myRank == 0) {
        cerr << endl << "New column Map:" << endl;
      }
      newColMap->describe (*errStream, Teuchos::VERB_EXTREME);
      comm->barrier ();
    }

    comm->barrier ();
    if (myRank == 2) {
      os << "Proc 2: checking global indices [10, 11, 12] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {10, 11, 12};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    else if (myRank == 3) {
      os << "Proc 3: checking global indices [15, 16, 14] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {15, 16, 14};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    comm->barrier ();

    TEST_ASSERT( graph.isFillComplete () );
    TEST_ASSERT( graph.isLocallyIndexed () );
    TEST_ASSERT( graph.hasColMap () );

    // reindexColumns() changes the graph, so we have to resume fill.
    graph.resumeFill ();

    TEST_ASSERT( ! graph.isFillComplete () );
    TEST_ASSERT( graph.isLocallyIndexed () );
    TEST_ASSERT( graph.hasColMap () );

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    const bool sortGraph = true;

    // Call the reindexColumns() method: the moment of truth!
    try {
      graph.reindexColumns (newColMap, Teuchos::null, sortGraph);
    } catch (std::exception& e) {
      success = false;
      os << "Proc " << myRank << ": reindexColumns() threw an exception: "
         << e.what () << endl;
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Now call fillComplete to compute the new Import, if necessary.
    graph.fillComplete ();

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Does the graph now have the right column Map?
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    // FIXME (mfh 18 Aug 2014) Some of these tests may hang if the
    // graph's column Map is null on some, but not all processes.
    if (! graph.getColMap ().is_null ()) {
      TEST_ASSERT( graph.getColMap ()->isSameAs (*newColMap) );
    }

    RCP<const import_type> theImport = graph.getImporter ();

    // If there's only one process in the communicator, the graph
    // won't have an Import object.  But if there's more than one
    // process, this particular graph should have one.
    TEST_ASSERT( comm->getSize () == 1 || ! theImport.is_null () );
    // FIXME (mfh 18 Aug 2014) Some of these tests may hang if the
    // graph's Import object is null on some, but not all processes.
    if (! theImport.is_null ()) {
      TEST_ASSERT( ! theImport->getSourceMap ().is_null () );
      TEST_ASSERT( ! theImport->getTargetMap ().is_null () );
      if (! theImport->getSourceMap ().is_null ()) {
        if (! graph.getDomainMap ().is_null ()) {
          TEST_ASSERT( theImport->getSourceMap ()->isSameAs (* (graph.getDomainMap ())) );
        }
      }
      if (! theImport->getTargetMap ().is_null ()) {
        TEST_ASSERT( theImport->getTargetMap ()->isSameAs (* newColMap) );
      }
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Check that all the graph's indices are correct.  We know that
    // this is a local-only column Map change, so we don't have to
    // call getRemoteIndexList to do this; just convert the local
    // indices to global in the new column Map, and then back to local
    // in the old column Map, and compare with those in the original
    // graph.
    const LO myNumRows = static_cast<LO> (rowMap->getLocalNumElements ());
    if (myNumRows > 0) {
      for (LO lclRowInd = 0; lclRowInd < myNumRows; ++lclRowInd) {
        os << "Proc " << myRank << ": Row: " << lclRowInd;

        if (graph.getNumEntriesInLocalRow (lclRowInd) !=
            graph2->getNumEntriesInLocalRow (lclRowInd)) {
          success = false;
          os << ": # entries differ: "
             << graph.getNumEntriesInLocalRow (lclRowInd) << " != "
             << graph2->getNumEntriesInLocalRow (lclRowInd) << endl;
          continue; // don't even bother with the rest
        }
        const size_type numEnt =
          static_cast<size_type> (graph.getNumEntriesInLocalRow (lclRowInd));

        // Get the "new" local column indices that resulted from the
        // call to reindexColumns.  Get by copy, not by view, so we
        // can sort it.
        lids_type newLclColInds ("newLclColIds",numEnt);
        {
          size_t actualNumEnt = 0;
          graph.getLocalRowCopy (lclRowInd, newLclColInds, actualNumEnt);
          if (static_cast<size_t> (numEnt) != actualNumEnt) {
            os << ", graph.getLocalRowCopy(...) reported different # entries"
               << endl;
            success = false;
            continue; // don't even bother with the rest
          }
        }
        //        os << ", newLclInds: " << Teuchos::toString (newLclColInds);

        // Use the new column Map to convert them to global indices.
        Array<GO> gblColInds (numEnt);
        for (size_type k = 0; k < numEnt; ++k) {
          if (newLclColInds[k] == Teuchos::OrdinalTraits<LO>::invalid ()) {
            success = false;
          }
          gblColInds[k] = newColMap->getGlobalElement (newLclColInds[k]);
          if (gblColInds[k] == Teuchos::OrdinalTraits<GO>::invalid ()) {
            success = false;
          }
        }
        os << ", gblInds: " << Teuchos::toString (gblColInds ());

        // Convert those global indices to the original column Map's
        // local indices.  Those should match the original local
        // indices in the (cloned) original graph.
        Array<LO> oldLclColInds (numEnt);
        for (size_type k = 0; k < numEnt; ++k) {
          const GO gblColInd = gblColInds[k];
          if (! curColMap->isNodeGlobalElement (gblColInd)) {
            os << ", " << gblColInd << " NOT in curColMap!";
            success = false;
          }
          if (! newColMap->isNodeGlobalElement (gblColInd)) {
            os << ", " << gblColInd << " NOT in newColMap!";
            success = false;
          }
          oldLclColInds[k] = curColMap->getLocalElement (gblColInd);
          if (oldLclColInds[k] == Teuchos::OrdinalTraits<LO>::invalid ()) {
            success = false;
          }
        }
        os << ", oldLclInds: " << Teuchos::toString (oldLclColInds);

        // Get the original local indices from the original graph.
        lids_type origLclColInds ("origLclColIds",numEnt);
        {
          size_t actualNumEnt = 0;
          graph2->getLocalRowCopy (lclRowInd, origLclColInds, actualNumEnt);
          if (static_cast<size_t> (numEnt) != actualNumEnt) {
            os << ", graph2.getLocalRowCopy(...) reported different # entries"
               << endl;
            success = false;
            continue; // don't even bother with the rest
          }
        }
        //        os << ", origLclInds: " << Teuchos::toString (origLclColInds);

        // The indices in both graphs don't need to be in the same
        // order; they just need to be the same indices.
        Tpetra::sort (origLclColInds, origLclColInds.extent(0));
        std::sort (oldLclColInds.begin (), oldLclColInds.end ());

        // Compare the two sets of indices.
        bool arraysSame = true;
        if ((size_t)oldLclColInds.size() != (size_t)origLclColInds.extent(0)) {
          arraysSame = false;
        } else {
          for (size_type k = 0; k < oldLclColInds.size (); ++k) {
            if (oldLclColInds[k] != origLclColInds[k]) {
              arraysSame = false;
            }
          }
        }
        if (! arraysSame) {
          success = false;
          os << ", WRONG!";
        }
        os << endl;
      }
    }

    comm->barrier ();
    os << endl;
    if (myRank == 2) {
      os << "Proc 2: checking global indices [10, 11, 12] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {10, 11, 12};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    else if (myRank == 3) {
      os << "Proc 3: checking global indices [15, 16, 14] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {15, 16, 14};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    comm->barrier ();

    if (false) {
      comm->barrier ();
      RCP<Teuchos::FancyOStream> errStream =
        Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));

      if (myRank == 0) {
        cerr << "Original column Map:" << endl;
      }
      curColMap->describe (*errStream, Teuchos::VERB_EXTREME);

      if (myRank == 0) {
        cerr << endl << "New column Map:" << endl;
      }
      newColMap->describe (*errStream, Teuchos::VERB_EXTREME);
      comm->barrier ();
    }

    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << os.str ();
      }
      comm->barrier (); // let output complete
      comm->barrier ();
      comm->barrier ();
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  //
  // Test CrsGraph::reindexColumns, with both the column Map and the
  // Import provided, and with sorting off.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphReindexColumns, ColMapAndImportSortingOff, LO, GO, Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    using lids_type = typename graph_type::nonconst_local_inds_host_view_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    int gblSuccess = 0;
    int lclSuccess = 1;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Collect all the debug output on each process here, so we can
    // print it coherently across processes at the end.
    std::ostringstream os;

    // Create the graph's row Map.
    // const size_t numLocalIndices = 5;
    const GST numGlobalIndices = static_cast<GST> (5 * numProcs);
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (numGlobalIndices, indexBase, comm,
                         Tpetra::GloballyDistributed));
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! rowMap->isContiguous (), std::logic_error, "The row Map is supposed "
      "to be contiguous, but is not.");

    const size_t maxNumEntPerRow = 3;
    graph_type graph (rowMap, maxNumEntPerRow);

    // Make the usual tridiagonal graph.  Let the graph create its own
    // column Map.  We'll use that column Map to create a new column
    // Map, and give that column Map to reindexColumns().

    if (rowMap->getLocalNumElements () > 0) {
      const GO myMinGblInd = rowMap->getMinGlobalIndex ();
      const GO myMaxGblInd = rowMap->getMaxGlobalIndex ();
      const GO gblMinGblInd = rowMap->getMinAllGlobalIndex ();
      const GO gblMaxGblInd = rowMap->getMaxAllGlobalIndex ();

      Array<GO> gblInds (maxNumEntPerRow);
      for (GO gblRowInd = myMinGblInd; gblRowInd <= myMaxGblInd; ++gblRowInd) {
        size_t numInds = 0;
        if (gblRowInd == gblMinGblInd) {
          if (gblRowInd < gblMaxGblInd) {
            numInds = 2;
            gblInds[0] = gblRowInd;
            gblInds[1] = gblRowInd + 1;
          } else { // special case of 1 x 1 graph
            numInds = 1;
            gblInds[0] = gblRowInd;
          }
        } else if (gblRowInd == gblMaxGblInd) {
          if (gblRowInd > gblMinGblInd) {
            numInds = 2;
            gblInds[0] = gblRowInd - 1;
            gblInds[1] = gblRowInd;
          } else { // special case of 1 x 1 graph
            numInds = 1;
            gblInds[0] = gblRowInd;
          }
        } else {
          numInds = 3;
          gblInds[0] = gblRowInd - 1;
          gblInds[1] = gblRowInd;
          gblInds[2] = gblRowInd + 1;
        }
        ArrayView<const GO> gblIndsView = gblInds (0, numInds);
        graph.insertGlobalIndices (gblRowInd, gblIndsView);
      }
    }
    graph.fillComplete ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isFillComplete (), std::logic_error, "The graph claims that it "
      "is not fill complete, after fillComplete was called.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.hasColMap () || graph.getColMap ().is_null (), std::logic_error,
      "The graph is fill complete, but doesn't have a column Map.");

    // Make a deep copy of the graph, to check later that the
    // conversion was correct.
    RCP<graph_type> graph2 = rcp(new graph_type(graph));
    TEST_ASSERT( ! graph2.is_null () );

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Create a new column Map, which has all the global indices on
    // each process (_locally_) in reverse order of the graph's
    // current column Map.

    // NOTE (mfh 21 Aug 2014) Don't get this as a reference!  Get it
    // as an RCP!  Remember that the reference comes from the result
    // of getColMap(), and the call to reindexColumns() will
    // invalidate the graph's current column Map.
    RCP<const map_type> curColMap = graph.getColMap ();
    Array<GO> newGblInds (curColMap->getLocalNumElements ());
    if (curColMap->isContiguous ()) {
      // const GO myMinGblInd = curColMap->getMinGlobalIndex ();
      const GO myMaxGblInd = curColMap->getMaxGlobalIndex ();

      const size_type myNumInds =
        static_cast<size_type> (curColMap->getLocalNumElements ());
      if (myNumInds > 0) {
        GO curGblInd = myMaxGblInd;
        for (size_type k = 0; k < myNumInds; ++k, --curGblInd) {
          newGblInds[k] = curGblInd;
        }
      }
    } else { // original column Map is not contiguous
      ArrayView<const GO> curGblInds = curColMap->getLocalElementList ();
      for (size_type k = 0; k < curGblInds.size (); ++k) {
        const size_type k_opposite = (newGblInds.size () - 1) - k;
        newGblInds[k] = curGblInds[k_opposite];
      }
    }

    RCP<const map_type> newColMap =
      rcp (new map_type (INVALID, newGblInds (), indexBase, comm));

    // Print both old and new column Maps.
    {
      comm->barrier ();
      RCP<Teuchos::FancyOStream> errStream =
        Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));

      if (myRank == 0) {
        cerr << "Original column Map:" << endl;
      }
      curColMap->describe (*errStream, Teuchos::VERB_EXTREME);

      if (myRank == 0) {
        cerr << endl << "New column Map:" << endl;
      }
      newColMap->describe (*errStream, Teuchos::VERB_EXTREME);
      comm->barrier ();
    }

    comm->barrier ();
    if (myRank == 2) {
      os << "Proc 2: checking global indices [10, 11, 12] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {10, 11, 12};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    else if (myRank == 3) {
      os << "Proc 3: checking global indices [15, 16, 14] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {15, 16, 14};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    comm->barrier ();

    TEST_ASSERT( graph.isFillComplete () );
    TEST_ASSERT( graph.isLocallyIndexed () );
    TEST_ASSERT( graph.hasColMap () );

    // reindexColumns() changes the graph, so we have to resume fill.
    graph.resumeFill ();

    TEST_ASSERT( ! graph.isFillComplete () );
    TEST_ASSERT( graph.isLocallyIndexed () );
    TEST_ASSERT( graph.hasColMap () );

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Create the new Import object (whose target Map is the new column Map).
    RCP<import_type> newImport =
      rcp (new import_type (graph.getDomainMap (), newColMap));

    const bool sortGraph = false;

    // Call the reindexColumns() method: the moment of truth!
    try {
      graph.reindexColumns (newColMap, newImport, sortGraph);
    } catch (std::exception& e) {
      success = false;
      os << "Proc " << myRank << ": reindexColumns() threw an exception: "
         << e.what () << endl;
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Now call fillComplete to compute the new Import, if necessary.
    graph.fillComplete ();

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Does the graph now have the right column Map?
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    // FIXME (mfh 18 Aug 2014) Some of these tests may hang if the
    // graph's column Map is null on some, but not all processes.
    if (! graph.getColMap ().is_null ()) {
      TEST_ASSERT( graph.getColMap ()->isSameAs (*newColMap) );
    }

    RCP<const import_type> theImport = graph.getImporter ();

    // If there's only one process in the communicator, the graph
    // won't have an Import object.  But if there's more than one
    // process, this particular graph should have one.
    TEST_ASSERT( comm->getSize () == 1 || ! theImport.is_null () );
    // FIXME (mfh 18 Aug 2014) Some of these tests may hang if the
    // graph's Import object is null on some, but not all processes.
    if (! theImport.is_null ()) {
      TEST_ASSERT( ! theImport->getSourceMap ().is_null () );
      TEST_ASSERT( ! theImport->getTargetMap ().is_null () );
      if (! theImport->getSourceMap ().is_null ()) {
        if (! graph.getDomainMap ().is_null ()) {
          TEST_ASSERT( theImport->getSourceMap ()->isSameAs (* (graph.getDomainMap ())) );
        }
      }
      if (! theImport->getTargetMap ().is_null ()) {
        TEST_ASSERT( theImport->getTargetMap ()->isSameAs (* newColMap) );
      }
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Check that all the graph's indices are correct.  We know that
    // this is a local-only column Map change, so we don't have to
    // call getRemoteIndexList to do this; just convert the local
    // indices to global in the new column Map, and then back to local
    // in the old column Map, and compare with those in the original
    // graph.
    const LO myNumRows = static_cast<LO> (rowMap->getLocalNumElements ());
    if (myNumRows > 0) {
      for (LO lclRowInd = 0; lclRowInd < myNumRows; ++lclRowInd) {
        os << "Proc " << myRank << ": Row: " << lclRowInd;

        if (graph.getNumEntriesInLocalRow (lclRowInd) !=
            graph2->getNumEntriesInLocalRow (lclRowInd)) {
          success = false;
          os << ": # entries differ: "
             << graph.getNumEntriesInLocalRow (lclRowInd) << " != "
             << graph2->getNumEntriesInLocalRow (lclRowInd) << endl;
          continue; // don't even bother with the rest
        }
        const size_type numEnt =
          static_cast<size_type> (graph.getNumEntriesInLocalRow (lclRowInd));

        // Get the "new" local column indices that resulted from the
        // call to reindexColumns.  Get by copy, not by view, so we
        // can sort it.
        lids_type newLclColInds ("newLclColInds",numEnt);
        {
          size_t actualNumEnt = 0;
          graph.getLocalRowCopy (lclRowInd, newLclColInds, actualNumEnt);
          if (static_cast<size_t> (numEnt) != actualNumEnt) {
            os << ", graph.getLocalRowCopy(...) reported different # entries"
               << endl;
            success = false;
            continue; // don't even bother with the rest
          }
        }
        //        os << ", newLclInds: " << Teuchos::toString (newLclColInds);

        // Use the new column Map to convert them to global indices.
        Array<GO> gblColInds (numEnt);
        for (size_type k = 0; k < numEnt; ++k) {
          if (newLclColInds[k] == Teuchos::OrdinalTraits<LO>::invalid ()) {
            success = false;
          }
          gblColInds[k] = newColMap->getGlobalElement (newLclColInds[k]);
          if (gblColInds[k] == Teuchos::OrdinalTraits<GO>::invalid ()) {
            success = false;
          }
        }
        os << ", gblInds: " << Teuchos::toString (gblColInds ());

        // Convert those global indices to the original column Map's
        // local indices.  Those should match the original local
        // indices in the (cloned) original graph.
        Array<LO> oldLclColInds (numEnt);
        for (size_type k = 0; k < numEnt; ++k) {
          const GO gblColInd = gblColInds[k];
          if (! curColMap->isNodeGlobalElement (gblColInd)) {
            os << ", " << gblColInd << " NOT in curColMap!";
            success = false;
          }
          if (! newColMap->isNodeGlobalElement (gblColInd)) {
            os << ", " << gblColInd << " NOT in newColMap!";
            success = false;
          }
          oldLclColInds[k] = curColMap->getLocalElement (gblColInd);
          if (oldLclColInds[k] == Teuchos::OrdinalTraits<LO>::invalid ()) {
            success = false;
          }
        }
        os << ", oldLclInds: " << Teuchos::toString (oldLclColInds);

        // Get the original local indices from the original graph.
        lids_type origLclColInds("origLclColIds",numEnt);
        {
          size_t actualNumEnt = 0;
          graph2->getLocalRowCopy (lclRowInd, origLclColInds, actualNumEnt);
          if (static_cast<size_t> (numEnt) != actualNumEnt) {
            os << ", graph2.getLocalRowCopy(...) reported different # entries"
               << endl;
            success = false;
            continue; // don't even bother with the rest
          }
        }
        //os << ", origLclInds: " << Teuchos::toString (origLclColInds);

        // The indices in both graphs don't need to be in the same
        // order; they just need to be the same indices.
        Tpetra::sort (origLclColInds, origLclColInds.extent(0));
        std::sort (oldLclColInds.begin (), oldLclColInds.end ());

        // Compare the two sets of indices.
        bool arraysSame = true;
        if ((size_t)oldLclColInds.size () != (size_t)origLclColInds.size ()) {
          arraysSame = false;
        } else {
          for (size_type k = 0; k < oldLclColInds.size (); ++k) {
            if (oldLclColInds[k] != origLclColInds[k]) {
              arraysSame = false;
            }
          }
        }
        if (! arraysSame) {
          success = false;
          os << ", WRONG!";
        }
        os << endl;
      }
    }

    comm->barrier ();
    os << endl;
    if (myRank == 2) {
      os << "Proc 2: checking global indices [10, 11, 12] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {10, 11, 12};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    else if (myRank == 3) {
      os << "Proc 3: checking global indices [15, 16, 14] "
        "(should be owned on this process)" << endl;
      const LO testGids[] = {15, 16, 14};
      const LO numTestGids = 3;
      for (LO k = 0; k < numTestGids; ++k) {
        const GO gid = testGids[k];
        const GO lid_old = curColMap->getLocalElement (gid);
        const GO lid_new = newColMap->getLocalElement (gid);

        os << "  gbl: " << gid
           << ", gbl->lcl_old: " << lid_old
           << ", gbl->lcl_new: " << lid_new
           << ", gbl->lcl_old->gbl: " << curColMap->getGlobalElement (lid_old)
           << ", gbl->lcl_new->gbl: " << newColMap->getGlobalElement (lid_new)
           << ", gbl->lcl_old->gbl->lcl_new: "
           << newColMap->getLocalElement (curColMap->getGlobalElement (lid_old))
           << ", gbl->lcl_new->gbl->lcl_old: "
           << curColMap->getLocalElement (newColMap->getGlobalElement (lid_new))
           << endl;
      }
    }
    comm->barrier ();

    if (false) {
      comm->barrier ();
      RCP<Teuchos::FancyOStream> errStream =
        Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));

      if (myRank == 0) {
        cerr << "Original column Map:" << endl;
      }
      curColMap->describe (*errStream, Teuchos::VERB_EXTREME);

      if (myRank == 0) {
        cerr << endl << "New column Map:" << endl;
      }
      newColMap->describe (*errStream, Teuchos::VERB_EXTREME);
      comm->barrier ();
    }

    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << os.str ();
      }
      comm->barrier (); // let output complete
      comm->barrier ();
      comm->barrier ();
    }

    gblSuccess = 0;
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


//
// INSTANTIATIONS
//

// Tests to build and run in both debug and release modes.  We will
// instantiate them over all enabled local ordinal (LO), global
// ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphReindexColumns, ColMapOnlySortingOn, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphReindexColumns, ColMapAndImportSortingOff, LO, GO, NODE )


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

}


