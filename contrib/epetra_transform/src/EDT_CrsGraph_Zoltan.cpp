
#include <vector>
#include <set>

#include <EDT_CrsGraph_Zoltan.h>

#include <Epetra_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <Epetra_MpiComm.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_Zoltan::operator()( const Epetra_CrsGraph & original )
{
  int err;

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan_LoadBalance  LB( 0, &dummy, &version );
  err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(original.Comm()).Comm() );
  if( err == LB_OK ) err = LB.Set_Method( "PARMETIS" );
  char * pM = new char[partitionMethod_.size()];
  memcpy( pM, partitionMethod_.c_str(), partitionMethod_.size() );
  if( err == LB_OK ) err = LB.Set_Param( "PARMETIS_METHOD", pM );

  //Setup Query Object
  Epetra_ZoltanQuery Query( original );
  if( err == LB_OK ) err = LB.Set_QueryObject( &Query );

  if( err != LB_OK )
  { cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; exit(0); }

  //Generate Load Balance
  int changes;
  int num_gid_entries, num_lid_entries;
  int num_import;
  LB_ID_PTR import_global_ids, import_local_ids;
  int * import_procs;
  int num_export;
  LB_ID_PTR export_global_ids, export_local_ids;
  int * export_procs;

  original.Comm().Barrier();
  err = LB.Balance( &changes,
                    &num_gid_entries, &num_lid_entries,
                    &num_import, &import_global_ids, &import_local_ids, &import_procs,
                    &num_export, &export_global_ids, &export_local_ids, &export_procs );
  original.Comm().Barrier();

  //Generate New Element List
  int numMyElements = original.RowMap().NumMyElements();
  vector<int> elementList( numMyElements );
  original.RowMap().MyGlobalElements( &elementList[0] );

  int newNumMyElements = numMyElements - num_export + num_import;
  vector<int> newElementList( newNumMyElements );

  set<int> gidSet;
  for( int i = 0; i < num_export; ++i )
    gidSet.insert( export_global_ids[i] );

  //Add unmoved indices to new list
  int loc = 0;
  for( int i = 0; i < numMyElements; ++i )
    if( !gidSet.count( elementList[i] ) )
      newElementList[loc++] = elementList[i];
  
  //Add imports to end of list
  for( int i = 0; i < num_import; ++i )
    newElementList[loc+i] = import_global_ids[i];

  //Free Zoltan Data
  if( err == LB_OK )
    err = LB.Free_Data( &import_global_ids, &import_local_ids, &import_procs,
                        &export_global_ids, &export_local_ids, &export_procs );

  //Create Import Map
  Epetra_Map * ImportMap = new Epetra_Map( original.RowMap().NumGlobalElements(),
                                           newNumMyElements,
                                           &newElementList[0],
                                           original.RowMap().IndexBase(),
                                           original.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *ImportMap, original.RowMap() );

  //Create New Graph
  std::auto_ptr<Epetra_CrsGraph> NewGraph( new Epetra_CrsGraph( Copy, *ImportMap, 0 ) );
  NewGraph->Import( original, Importer, Insert );
  NewGraph->TransformToLocal();
  
  return NewGraph;
}

}
