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
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"

namespace { // (anonymous)
  TEUCHOS_UNIT_TEST(CrsGraph, DetectNonlocals)
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using crs_graph_type = Tpetra::CrsGraph<>;
    using map_type = Tpetra::Map<>;
    using LO = crs_graph_type::local_ordinal_type;
    using GO = crs_graph_type::global_ordinal_type;

    const bool debug = Tpetra::Details::Behavior::debug("CrsGraph");
    if (! debug) {
      out << "This test only does anything in debug mode.  "
        "Try rerunning with the environment variable "
        "TPETRA_DEBUG=CrsGraph set." << endl;
      return;
    }
    out << "Test that Tpetra::CrsGraph::fillComplete reports an "
      "error if the \"No Nonlocal Changes\" parameter is true, but "
      "there were nonlocal inserts." << endl;

    auto comm = getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    TEST_EQUALITY( numProcs, 2 );

    const LO lclNumRows = 1;
    const GO gblNumRows = 2;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp(new map_type(gblNumRows, lclNumRows, indexBase, comm));

    for (const bool assertNoNonlocalChanges : {false, true}) {
      crs_graph_type G(rowMap, 1);

      // Attempt an insert into a nonowned row, but only on Process 1.
      // In debug mode, the processes should do an all-reduce to check
      // whether any process had a nonowned insert.
      if (myRank == 1) {
        const GO gblRow(0);
        const GO gblColInd(0);
        TEST_ASSERT( ! rowMap->isNodeGlobalElement(gblRow) );
        G.insertGlobalIndices(gblRow, LO(1), &gblColInd);
      }

      auto params = Teuchos::parameterList("CrsGraph::fillComplete");
      // If set to true, this is a lie; the point is that in debug
      // mode, fillComplete should check whether it is a lie.  Don't
      // use TEST_THROW because that requires a specific type of
      // exception; we don't care about the exception type here, as
      // long as it's a subclass of std::exception.
      if (assertNoNonlocalChanges) {
        params->set("No Nonlocal Changes", true);
      }
      bool threw = false;
      try {
        G.fillComplete(params);
      }
      catch (std::exception& e) {
        threw = true;
        out << "G.fillComplete() threw: " << e.what() << endl;
      }
      if (assertNoNonlocalChanges) {
        TEST_ASSERT( threw );
      }
      else {
        TEST_ASSERT( ! threw );
      }

      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 0;
      using Teuchos::outArg;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;
      reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
    }
  }
} // namespace (anonymous)
