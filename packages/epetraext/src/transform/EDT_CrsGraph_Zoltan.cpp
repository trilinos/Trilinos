
#include <EDT_CrsGraph_Zoltan.h>

#ifdef ZOLTAN_ORDER
#include <EDT_CrsGraph_ZoltanOrder.h>
#endif

#include <EDT_CrsGraph_Transpose.h>

#include <Epetra_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <Epetra_MpiComm.h>

#include <vector>
#include <set>

using std::vector;
using std::set;

namespace EpetraExt {
namespace Transform {

NewTypePtr CrsGraph_Zoltan::operator()( OriginalTypeRef original )
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
    if( err == ZOLTAN_OK ) err = LB->Set_Param( "LB_METHOD", "PARMETIS" );
    if( err == ZOLTAN_OK ) err = LB->Set_Param( "PARMETIS_METHOD", partitionMethod_ );

    //Setup Query Object
    Epetra_CrsGraph * TGraph = CrsGraph_Transpose()( original );
    TransGraph = TGraph.release();
    Query = new Epetra_ZoltanQuery( original, TransGraph );
    if( err == ZOLTAN_OK ) err = LB->Set_QueryObject( Query );

    if( err != ZOLTAN_OK )
    { cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; exit(0); }
  }

  //Generate Load Balance
  int changes;
  int num_gid_entries, num_lid_entries;
  int num_import;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  int * import_procs;
  int num_export;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;
  int * export_procs;

  original.Comm().Barrier();
  err = LB->Balance( &changes,
                     &num_gid_entries, &num_lid_entries,
                     &num_import, &import_global_ids, &import_local_ids, &import_procs,
                     &num_export, &export_global_ids, &export_local_ids, &export_procs );
  original.Comm().Barrier();

  if( TransGraph ) delete TransGraph;

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
  if( err == ZOLTAN_OK )
    err = LB->Free_Data( &import_global_ids, &import_local_ids, &import_procs,
                         &export_global_ids, &export_local_ids, &export_procs );

  if( !lb_ ) { delete LB; delete Query; }

  //Create Import Map
  Epetra_Map * ImportMap = new Epetra_Map( original.RowMap().NumGlobalElements(),
                                           newNumMyElements,
                                           &newElementList[0],
                                           original.RowMap().IndexBase(),
                                           original.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *ImportMap, original.RowMap() );

  //Create New Graph
  Epetra_CrsGraph * NewGraph( new Epetra_CrsGraph( Copy, *ImportMap, 0 ) );
  NewGraph->Import( original, Importer, Insert );
  NewGraph->TransformToLocal();

#ifdef ZOLTAN_ORDER
  if( reorder_ )
  {
    Epetra_CrsGraph * oldGraph = NewGraph;
    NewGraph = CrsGraph_ZoltanOrder()( *oldGraph );
    delete oldGraph;
    delete ImportMap;
  }
#endif
  
  return NewGraph;
}

} //namespace Transform
} //namespace EpetraExt

