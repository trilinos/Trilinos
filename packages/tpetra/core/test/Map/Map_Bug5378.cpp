// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_Map.hpp>

namespace { // (anonymous)

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;


////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Map, Bug5378_GoodGIDs, LO, GO)
{
  RCP<const Teuchos::Comm<int> > comm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(10,10,comm);
  Array<GO> lookup_gids(  tuple<GO>(1,3,5) );
  Array<int> expected_ids(  tuple<int>( 0,0,0) ); // MPI process ranks are int
  Array<LO> expected_lids( tuple<LO>( 1,3,5) );
  Array<int> nodeIDs( lookup_gids.size() ); // MPI process ranks are int
  Array<LO> nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::AllIDsPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Map, Bug5378_BadGIDs, LO, GO)
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(10,10,comm);
  Array<GO> lookup_gids(  tuple<GO>(1,10,5) );
  Array<int> expected_ids(  tuple<int>( 0,-1,0) ); // MPI process ranks are int
  Array<LO> expected_lids( tuple<LO>( 1,-1,5) );
  Array<int> nodeIDs( lookup_gids.size() ); // MPI process ranks are int
  Array<LO> nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Map, Bug5378_GoodGIDsNoLIDs, LO, GO)
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(10,10,comm);
  Array<GO> lookup_gids(  tuple<GO>(1,3,5) );
  Array<int> expected_ids(  tuple<int>( 0,0,0) ); // MPI process ranks are int
  Array<int> nodeIDs( lookup_gids.size() ); // MPI process ranks are int
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::AllIDsPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Map, Bug5378_BadGIDsNoLIDs, LO, GO)
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  //
  RCP<const Tpetra::Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(10,10,comm);
  Array<GO> lookup_gids(  tuple<GO>(1,10,5) );
  Array<int> expected_ids(  tuple<int>( 0,-1,0) ); // MPI process ranks are int
  Array<int> nodeIDs( lookup_gids.size() ); // MPI process ranks are int
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );

  // Process 0 is responsible for printing the "SUCCESS" / "PASSED"
  // message, so without the barrier, it's possible for the test to be
  // reported as passing, even if the other processes crashed or hung.
  comm->barrier ();
}

#define UNIT_TEST_GROUP(LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(Map, Bug5378_GoodGIDs, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(Map, Bug5378_BadGIDs, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(Map, Bug5378_GoodGIDsNoLIDs, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(Map, Bug5378_BadGIDsNoLIDs, LO, GO)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

} // namespace (anonymous)
