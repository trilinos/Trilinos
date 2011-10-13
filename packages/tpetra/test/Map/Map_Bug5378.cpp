#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include "Tpetra_Map.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#endif

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;

//// 
TEUCHOS_UNIT_TEST( Map, Bug5378_GoodGIDs ) 
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  // 
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,3,5) );
  Array<int> expected_ids(  tuple<int>( 0,0,0) );
  Array<int> expected_lids( tuple<int>( 1,3,5) );
  Array<int> nodeIDs( lookup_gids.size() ), 
             nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::AllIDsPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );
}

//// 
TEUCHOS_UNIT_TEST( Map, Bug5378_BadGIDs ) 
{
  RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD));
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6 7 8 9
  //
  // Lookup of any valid global ID should identify with proc 0
  // 
  RCP<const Tpetra::Map<int,long> > map = Tpetra::createContigMap<int,long>(10,10,comm);
  Array<long> lookup_gids(  tuple<long>(1,10,5) );
  Array<int> expected_ids(  tuple<int>( 0,-1,0) );
  Array<int> expected_lids( tuple<int>( 1,-1,5) );
  Array<int> nodeIDs( lookup_gids.size() ), 
             nodeLIDs( lookup_gids.size() );
  Tpetra::LookupStatus lookup = map->getRemoteIndexList( lookup_gids(), nodeIDs(), nodeLIDs() );
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( nodeIDs(), expected_ids() );
  TEST_COMPARE_ARRAYS( nodeLIDs(), expected_lids() );
}
