
#ifdef ZOLTAN_ORDER

#include <EDT_CrsGraph_ZoltanOrder.h>

#include <EDT_CrsGraph_Transpose.h>

#include <Epetra_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

#include <Epetra_MpiComm.h>

#include <vector>
#include <set>

EpetraExt::CrsGraph_ZoltanOrder::NewTypePtr EpetraExt::CrsGraph_ZoltanOrder::operator()( EpetraExt::CrsGraph_ZoltanOrder::OriginalTypeRef original )
{
  int err;

  Zoltan_LoadBalance * LB = lb_;
  Epetra_ZoltanQuery * Query = 0;
  Epetra_CrsGraph * TransGraph = 0;
  if( !LB )
  {
    //Setup Load Balance Object
    float version;
    char * dummy = 0;
    LB = new Zoltan_LoadBalance( 0, &dummy, &version );
    err = LB->Create( dynamic_cast<const Epetra_MpiComm&>(original.Comm()).Comm() );
    LB->Set_Param( "ORDER_METHOD", "METIS" );
    LB->Set_Param( "ORDER_TYPE", "LOCAL" );

    //Setup Query Object
    TransGraph = CrsGraph_Transpose()( original );
    Query = new Epetra_ZoltanQuery( original, TransGraph, true );
    if( err == ZOLTAN_OK ) err = LB->Set_QueryObject( Query );

    if( err != ZOLTAN_OK )
    { cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; exit(0); }
  }

  //Generate Reorder
  int num_elements = original.RowMap().NumMyElements();
  int num_gid_entries, num_lid_entries;
  vector<ZOLTAN_ID_TYPE> global_ids(num_elements,0);
  vector<ZOLTAN_ID_TYPE> local_ids(num_elements,0);
  vector<int> rank(num_elements);
  vector<int> iperm(num_elements);

  original.Comm().Barrier();
  err = LB->Order( &num_gid_entries,
                   &num_lid_entries,
                   num_elements,
                   &global_ids[0],
                   &local_ids[0],
                   &rank[0],
                   &iperm[0] );
  original.Comm().Barrier();

#ifdef EDT_ZOLTANORDER_DEBUG
cout << "------------------------------\n";
cout << "#GIDs: " << num_gid_entries << endl;
cout << "#LIDs: " << num_lid_entries << endl;
cout << "GIDs and LIDs\n";
for( int i = 0; i < num_elements; ++i ) cout << global_ids[i] << " " << local_ids[i] << endl;
cout << "Rank and Perm\n";
for( int i = 0; i < num_elements; ++i ) cout << rank[i] << " " << iperm[i] << endl;
cout << "------------------------------\n";
#endif

  if( !lb_ ) { delete LB; delete Query; }

  if( TransGraph ) delete TransGraph;

  //Generate New Row Map
  vector<int> gids(num_elements);
  for( int i = 0; i < num_elements; ++i )
    gids[ (num_elements-1) - rank[i] ] = global_ids[i];
  Epetra_Map * NewRowMap = new Epetra_Map( original.RowMap().NumGlobalElements(),
                                           num_elements,
                                           &gids[0],
                                           original.RowMap().IndexBase(),
                                           original.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *NewRowMap, original.RowMap() );

  //Create New Graph
  Epetra_CrsGraph * NewGraph( new Epetra_CrsGraph( Copy, *NewRowMap, 0 ) );
  NewGraph->Import( original, Importer, Insert );
  NewGraph->TransformToLocal();
  
  return NewGraph;
}

#endif //ZOLTAN_ORDER
