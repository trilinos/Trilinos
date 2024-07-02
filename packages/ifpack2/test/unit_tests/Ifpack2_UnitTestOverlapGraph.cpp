// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestOverlapGraph.cpp

\brief Ifpack2 Unit Test for OverlapGraph class.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_OverlapGraph.hpp>

namespace {

using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::reduceAll;
using Teuchos::REDUCE_MIN;
using std::endl;
typedef tif_utest::Node Node;

// Macro used inside the unit test below.  It tests for global error,
// and if so, prints each process' error message and quits the test
// early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define IFPACK2OVERLAPGRAPH_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
  TEST_EQUALITY_CONST( gblSuccess, 1 ); \
  if (gblSuccess != 1) { \
    out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
    for (int p = 0; p < numProcs; ++p) { \
      if (myRank == p && lclSuccess != 1) { \
        std::cout << errStrm.str () << std::flush; \
      } \
      comm->barrier (); \
      comm->barrier (); \
      comm->barrier (); \
    } \
    return; \
  } \
} while (false)


//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2OverlapGraph, OverlapGraphTest0, LO, GO)
{
  typedef Tpetra::CrsGraph<LO,GO,Node> crs_graph_type;
  typedef Ifpack2::OverlapGraph<LO,GO,Node> overlap_graph_type;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::OverlapGraph unit test" << endl;
  RCP<const Comm<int> > comm = tif_utest::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // Create a test graph.
  const size_t num_rows_per_proc = 5;
  RCP<const crs_graph_type> crsgraph;
  try {
    crsgraph = tif_utest::create_tridiag_graph<LO,GO,Node> (num_rows_per_proc);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "create_tridiag_graph threw exception: " << e.what () << endl;
  }
  IFPACK2OVERLAPGRAPH_REPORT_GLOBAL_ERR( "Tpetra::CrsGraph creation" );

  TEST_EQUALITY( crsgraph->getMap()->getLocalNumElements(), num_rows_per_proc )

  const LO overlap_levels = 2;
  RCP<overlap_graph_type> overlapgraph;
  try {
    overlapgraph = rcp (new overlap_graph_type (crsgraph, overlap_levels));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::OverlapGraph constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2OVERLAPGRAPH_REPORT_GLOBAL_ERR( "Ifpack2::OverlapGraph constructor" );


  // Now test how many local rows there are in the overlapped-graph.
  // For a tri-diagonal input-graph with 5 local rows on each proc:
  // 'end' procs (procs 0 and numProcs-1) should have an extra row in the
  // overlapped graph for each level of overlap. Other procs should
  // have an extra 2 rows for each level of overlap.

  // Special case: if numProcs==1 then overlapgraph should have the same
  // number of rows as the input-graph.

  // NOTE: The Teuchos unit test macros actually only work on Process
  // 0 in MPI_COMM_WORLD.  This is because they work by printing a
  // failure message to 'out', but 'out' is set up only to print on
  // Process 0.

  if (numProcs == 1) {
    // This works fine, because Process 0 is doing the check in this case.
    TEST_EQUALITY(overlapgraph->getOverlapGraph().getMap()->getLocalNumElements(), num_rows_per_proc)
  }
  else { // numProcs > 1
    const size_t actualOverlap =
      overlapgraph->getOverlapGraph ().getMap ()->getLocalNumElements ();
    size_t expectedOverlap = 0;
    if (myRank == 0 || myRank == numProcs - 1) {
      expectedOverlap = num_rows_per_proc + overlap_levels;
    } else {
      expectedOverlap = num_rows_per_proc + overlap_levels * 2;
    }
    lclSuccess = (actualOverlap == expectedOverlap) ? 1 : 0;
    if (lclSuccess != 1) {
      errStrm << "Process " << myRank << ": Actual overlap " << actualOverlap
              << " != expected overlap " << expectedOverlap << "." << endl;
    }
    IFPACK2OVERLAPGRAPH_REPORT_GLOBAL_ERR( "Overlap is wrong on at least one process." );
  }
}

#define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2OverlapGraph, OverlapGraphTest0, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

IFPACK2_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

} // namespace (anonymous)

