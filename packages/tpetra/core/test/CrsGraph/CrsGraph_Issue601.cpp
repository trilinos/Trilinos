// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp" // REDUCE_MIN, reduceAll
#include "Tpetra_TestingUtilities.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  // See Github Issue #601.  This test is really only useful if run with > 2 MPI processes.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, GlobalAssembleOverlappingRowMap, LO, GO, NT )
  {
    typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
    typedef Tpetra::Map<LO, GO, NT> map_type;
    int lclSuccess = 1;
    int gblSuccess = 0;

    out << "Test Tpetra::CrsGraph::globalAssemble with overlapping row Map "
      "(see Github Issue #601)" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    TEST_ASSERT( numProcs > 2 );
    if (! success) {
      out << "This test is not useful unless you run it with 3 or more MPI "
        "processes!" << endl;
      return;
    }
    out << "Create overlapping row Map" << endl;

    // Each process r gets 2 rows: r, and (r+1) % numProcs.
    // The resulting row Map is overlapping, because numProcs > 2.

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    Teuchos::Array<GO> myGblRowInds (2);
    myGblRowInds[0] = static_cast<GO> (myRank);
    myGblRowInds[1] = static_cast<GO> ((myRank + 1) % numProcs);
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (INV, myGblRowInds (), indexBase, comm));

    bool insertLocalEntryValues[] = { true, false };

    {

      Teuchos::OSTab tab2 (out);

      // Test both with no local entries before globalAssemble(), and
      // with some (one) local entries before globalAssemble().  This
      // exercises lazy allocation of entries, if that's what CrsGraph
      // does (it does this as of 09 Sep 2016, though I've been aiming
      // to get rid of lazy allocation for years).
      for (bool insertLocalEntry : insertLocalEntryValues) {
        out << "insertLocalEntry: " << (insertLocalEntry ? "true" : "false")
            << endl;
        Teuchos::OSTab tab3 (out);

        const size_t maxNumEntPerRow = static_cast<size_t> (insertLocalEntry ? 2 : 1);
        crs_graph_type G (rowMap, maxNumEntPerRow);

        const GO gblRow0 = static_cast<GO> (myRank);
        const GO gblRow1 = static_cast<GO> ((myRank + 1) % numProcs);

        if (insertLocalEntry) {
          // Insert the diagonal entry in each of the rows on this
          // process.  This results in duplication across processes.
          // globalAssemble() should do the right thing in that case.
          G.insertGlobalIndices (gblRow0, tuple<GO> (gblRow0));
          G.insertGlobalIndices (gblRow1, tuple<GO> (gblRow1));
        }
        // Insert an entry into row (myRank - 1) % numProcs.  Add
        // numProcs to (myRank - 1) before taking the mod, so it
        // doesn't turn out negative.  (C(++) doesn't promise that the
        // result of mod is nonnegative if the input is negative.)
        {
          const GO nonlocalGblRow =
            static_cast<GO> ((numProcs + myRank - 1) % numProcs);
          // Make sure the new column index is not a duplicate (on the
          // process that will receive it at globalAssemble).  This is
          // very much NOT a square graph (/ matrix).
          const GO gblCol = static_cast<GO> (numProcs) + nonlocalGblRow;
          TEST_NOTHROW( G.insertGlobalIndices (nonlocalGblRow, tuple<GO> (gblCol)) );
        }

        out << "Call G.globalAssemble()" << endl;
        // This is the moment of truth.
        TEST_NOTHROW( G.globalAssemble () );

        lclSuccess = success ? 1 : 0;
        gblSuccess = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEST_EQUALITY_CONST( gblSuccess, 1 );
        if (! success) {
          out << "Test FAILED on some process; returning early." << endl;
          return;
        }

        out << "Now make sure each of our process' rows are correct" << endl;

        // Test gblRow0
        {
          typename crs_graph_type::global_inds_host_view_type gblInds;
          G.getGlobalRowView (gblRow0, gblInds);

          const LO expectedNumEnt = static_cast<LO> (maxNumEntPerRow);
          TEST_EQUALITY( static_cast<LO> (gblInds.size ()), expectedNumEnt );
          if (static_cast<LO> (gblInds.size ()) == expectedNumEnt) {
            if (insertLocalEntry) {
              auto lclEntIter = std::find (gblInds.data(), 
                                           gblInds.data() + gblInds.extent(0),
                                           gblRow0);
              TEST_ASSERT( lclEntIter != gblInds.data() + gblInds.extent(0));
            }
            const GO gblCol0 = gblRow0 + static_cast<GO> (numProcs);
            auto nonlclEntIter = std::find (gblInds.data(), 
                                            gblInds.data() + gblInds.extent(0),
                                            gblCol0);
            TEST_ASSERT( nonlclEntIter != gblInds.data() + gblInds.extent(0));
          }
        }

        // Test gblRow1
        {
          typename crs_graph_type::global_inds_host_view_type gblInds;
          G.getGlobalRowView (gblRow1, gblInds);

          const LO expectedNumEnt = static_cast<LO> (maxNumEntPerRow);
          TEST_EQUALITY( static_cast<LO> (gblInds.size ()), expectedNumEnt );
          if (static_cast<LO> (gblInds.size ()) == expectedNumEnt) {
            if (insertLocalEntry) {
              auto lclEntIter = std::find (gblInds.data(), 
                                           gblInds.data() + gblInds.extent(0),
                                           gblRow1);
              TEST_ASSERT( lclEntIter != gblInds.data() + gblInds.extent(0) );
            }
            const GO gblCol1 = gblRow1 + static_cast<GO> (numProcs);
            auto nonlclEntIter = std::find (gblInds.data(), 
                                            gblInds.data() + gblInds.extent(0),
                                            gblCol1);
            TEST_ASSERT( nonlclEntIter != gblInds.data() + gblInds.extent(0) );
          }
        }

        // Test across processes.
        lclSuccess = success ? 1 : 0;
        gblSuccess = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEST_EQUALITY_CONST( gblSuccess, 1 );
        if (! success) {
          out << "Test FAILED on some process; returning early." << endl;
          return;
        }
      } // insertLocalEntry
    } 
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, GlobalAssembleOverlappingRowMap, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
