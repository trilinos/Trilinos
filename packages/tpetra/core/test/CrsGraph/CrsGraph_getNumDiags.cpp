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
#include "Tpetra_Details_getNumDiags.hpp"
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

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, getNumDiags, LO, GO, NT )
  {
    typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
    typedef Tpetra::Map<LO, GO, NT> map_type;
    int lclSuccess = 1;
    int gblSuccess = 0;

    out << "Test Tpetra::Details::{getLocalNumDiags, getGlobalNumDiags}" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = getDefaultComm ();
    //const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    out << "Create row Map" << endl;
    const LO lclNumRows = 5;
    const GO gblNumRows = static_cast<GO> (numProcs) * static_cast<GO> (lclNumRows);
    constexpr GO indexBase = 0;
    RCP<const map_type> rowMap (new map_type (gblNumRows, indexBase, comm));
    auto domMap = rowMap;
    auto ranMap = rowMap;

    // Diagonal graph, globally indexed until first fillComplete.
    {
      const GO actualLocalNumDiags = lclNumRows;
      const GO actualGlobalNumDiags = gblNumRows;

      crs_graph_type graph (rowMap, 1);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        const GO gblCol = gblRow;
        graph.insertGlobalIndices (gblRow, 1, &gblCol);
      }
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
      }
      graph.fillComplete (domMap, ranMap);
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
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

    // Diagonal graph, locally indexed until first fillComplete.
    {
      const GO actualLocalNumDiags = lclNumRows;
      const GO actualGlobalNumDiags = gblNumRows;

      RCP<const map_type> colMap = rowMap;
      crs_graph_type graph (rowMap, colMap, 1);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const LO lclCol = lclRow;
        graph.insertLocalIndices (lclRow, 1, &lclCol);
      }
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
      }
      graph.fillComplete (domMap, ranMap);
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
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

    // Locally indexed graph, deliberately constructed to have full
    // structural rank and one entry per row, but no diagonal entries.
    // The local row and column indices in each row are always the
    // same, but the global row and column indices in each row are
    // always different.
    {
      const GO actualLocalNumDiags {0};
      const GO actualGlobalNumDiags {0};

      const GO gblNumCols = gblNumRows;
      RCP<const map_type> colMap;
      {
        std::vector<GO> gblColInds (lclNumRows);
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const GO gblRow = rowMap->getGlobalElement (lclRow);
          const GO gblCol = (gblRow + GO (1)) % gblNumCols;
          gblColInds[lclRow] = gblCol; // always shifted over by one
        }
        colMap = rcp (new map_type (gblNumCols, gblColInds.data (),
                                    gblColInds.size (), indexBase, comm));
      }

      crs_graph_type graph (rowMap, colMap, 1);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        // Not actually a diagonal entry, since the global indices differ.
        const LO lclCol = lclRow;
        graph.insertLocalIndices (lclRow, 1, &lclCol);
      }
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
      }
      graph.fillComplete (domMap, ranMap);
      {
        const GO reportedLocalNumDiags =
          Tpetra::Details::getLocalNumDiags (graph);
        const GO reportedGlobalNumDiags =
          Tpetra::Details::getGlobalNumDiags (graph);
        TEST_EQUALITY_CONST( reportedLocalNumDiags, actualLocalNumDiags );
        TEST_EQUALITY_CONST( reportedGlobalNumDiags, actualGlobalNumDiags );
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
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, getNumDiags, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
