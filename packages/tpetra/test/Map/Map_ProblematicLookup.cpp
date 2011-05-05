#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;

//// 
TEUCHOS_UNIT_TEST( Map, ProblematicLookup ) 
{
  Teuchos::oblackholestream blackhole;
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int MyPid = comm->getRank();
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6  
  //          Processor 1: Global IDs =                    9 10 11 12 13 14 15
  //
  // Lookup of global IDs 7 8 should return IDNotFound
  // 
  RCP<const Tpetra::Map<int,int> > map;
  if (MyPid == 0) {
    Array<int> gids( tuple<int>(1) );
    map = Tpetra::createNonContigMap<int>( gids().getConst() , comm );
  }
  else {
    Array<int> gids( tuple<int>(3) );
    map = Tpetra::createNonContigMap<int>( gids().getConst() , comm );
  }
  Array<int> nodeIDs( 1 );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( tuple<int>(2), nodeIDs() );
  TEST_EQUALITY_CONST( map->isDistributed(), true )
  TEST_EQUALITY_CONST( map->isContiguous(), false )
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), tuple<int>(-1) );
}
