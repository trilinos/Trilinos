
#include <EDT_CrsGraph_Zoltan.h>

#include <EDT_CrsGraph_Transpose.h>

#include <Epetra_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

#include <Epetra_MpiComm.h>

#include <vector>
#include <set>

using std::vector;
using std::set;

namespace EpetraExt {

CrsGraph_Zoltan::
~CrsGraph_Zoltan()
{
  if( newObj_ ) delete newObj_;

  if( NewRowMap_ ) delete NewRowMap_;
}

CrsGraph_Zoltan::NewTypeRef
CrsGraph_Zoltan::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan_LoadBalance LB( 0, &dummy, &version );
  err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(orig.Comm()).Comm() );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "LB_METHOD", "PARMETIS" );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "PARMETIS_METHOD", partitionMethod_ );

  //Setup Query Object
  CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( orig );
  Epetra_ZoltanQuery Query( orig, &TransGraph );
  if( err == ZOLTAN_OK ) err = LB.Set_QueryObject( &Query );

  if( err != ZOLTAN_OK )
  { cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; exit(0); }

  //Generate Load Balance
  int changes;
  int num_gid_entries, num_lid_entries;
  int num_import;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  int * import_procs;
  int num_export;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;
  int * export_procs;

//  orig.Comm().Barrier();
//  err = LB.Generate_Files( "zoltan_output" );
  orig.Comm().Barrier();
  err = LB.Balance( &changes,
                     &num_gid_entries, &num_lid_entries,
                     &num_import, &import_global_ids, &import_local_ids, &import_procs,
                     &num_export, &export_global_ids, &export_local_ids, &export_procs );
  LB.Evaluate( 1, 0, 0, 0, 0, 0, 0 );
  orig.Comm().Barrier();

  //Generate New Element List
  int numMyElements = orig.RowMap().NumMyElements();
  vector<int> elementList( numMyElements );
  orig.RowMap().MyGlobalElements( &elementList[0] );

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
    err = LB.Free_Data( &import_global_ids, &import_local_ids, &import_procs,
                         &export_global_ids, &export_local_ids, &export_procs );

  //Create Import Map
  NewRowMap_ = new Epetra_Map( orig.RowMap().NumGlobalElements(),
                                           newNumMyElements,
                                           &newElementList[0],
                                           orig.RowMap().IndexBase(),
                                           orig.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *NewRowMap_, orig.RowMap() );

  //Create New Graph
  Epetra_CrsGraph * NewGraph = new Epetra_CrsGraph( Copy, *NewRowMap_, 0 );
  NewGraph->Import( orig, Importer, Insert );
  NewGraph->TransformToLocal();

  newObj_ = NewGraph;

  return *NewGraph;
}

} // namespace EpetraExt

