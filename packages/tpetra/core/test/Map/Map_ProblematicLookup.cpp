// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;

////
TEUCHOS_UNIT_TEST( Map, ProblematicLookup )
{
  using std::cerr;
  using std::endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6
  //          Processor 1: Global IDs =                    9 10 11 12 13 14 15
  //
  // Lookup of global IDs 7 8 should return IDNotFound
  //
  if (myRank == 0) {
    cerr << "Creating Map" << endl;
  }
  comm->barrier ();
  comm->barrier ();
  comm->barrier (); // Just to make sure output finishes.

  typedef Tpetra::Map<> map_type;
  typedef map_type::local_ordinal_type LO;
  typedef map_type::global_ordinal_type GO;
  typedef int rank_type;

  RCP<const map_type> map;
  if (myRank == 0) {
    Array<GO> gids (tuple<GO> (1));
    map = Tpetra::createNonContigMap<LO, GO> (gids ().getConst () , comm);
  }
  else {
    Array<GO> gids (tuple<GO> (3));
    map = Tpetra::createNonContigMap<LO, GO> (gids ().getConst (), comm);
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": created Map" << endl;
    cerr << os.str ();
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": calling getRemoteIndexList" << endl;
    cerr << os.str ();
  }

  Array<rank_type> processRanks (1);
  Tpetra::LookupStatus lookup = map->getRemoteIndexList (tuple<GO> (2), processRanks ());

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": getRemoteIndexList done" << endl;
    cerr << os.str ();
  }
  comm->barrier ();
  if (myRank == 0) {
    cerr << "getRemoteIndexList finished on all processes" << endl;
  }
  comm->barrier (); // Just to make sure output finishes.

  TEST_EQUALITY_CONST( map->isDistributed(), true )
  TEST_EQUALITY_CONST( map->isContiguous(), false )
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( processRanks(), tuple<rank_type>(-1) );
}

