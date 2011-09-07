// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "EpetraExt_config.h"

#include <EpetraExt_Zoltan_CrsGraph.h>

#include <EpetraExt_Transpose_CrsGraph.h>

#include <EpetraExt_ZoltanQuery.h>
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

EPETRAEXT_DEPRECATED
Zoltan_CrsGraph::
~Zoltan_CrsGraph()
{
  if( newObj_ ) delete newObj_;

  if( NewRowMap_ ) delete NewRowMap_;
}

EPETRAEXT_DEPRECATED
Zoltan_CrsGraph::NewTypeRef
Zoltan_CrsGraph::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan::LoadBalance LB( 0, &dummy, &version );
  err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(orig.Comm()).Comm() );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "LB_METHOD", "GRAPH" );
#ifdef HAVE_LIBPARMETIS
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "GRAPH_PACKAGE", "PARMETIS" );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "PARMETIS_METHOD", partitionMethod_ );
#endif

  //Setup Query Object
  CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( orig );
  ZoltanQuery Query( orig, &TransGraph );
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
  NewGraph->FillComplete();

  Zoltan::LoadBalance LB2( 0, &dummy, &version );
  err = LB2.Create( dynamic_cast<const Epetra_MpiComm&>(orig.Comm()).Comm() );
  if( err == ZOLTAN_OK ) err = LB2.Set_Param( "LB_METHOD", "GRAPH" );
#ifdef HAVE_LIBPARMETIS
  if( err == ZOLTAN_OK ) err = LB2.Set_Param( "GRAPH_PACKAGE", "PARMETIS" );
  if( err == ZOLTAN_OK ) err = LB2.Set_Param( "PARMETIS_METHOD", partitionMethod_ );
#endif
  CrsGraph_Transpose transTrans;
  Epetra_CrsGraph & trans2 = transTrans( *NewGraph );
  ZoltanQuery query( *NewGraph, &trans2 );
  if( err == ZOLTAN_OK ) err = LB2.Set_QueryObject( &query );
  //err = LB2.Balance( &changes,
  //                   &num_gid_entries, &num_lid_entries,
  //                   &num_import, &import_global_ids, &import_local_ids, &import_procs,
  //                   &num_export, &export_global_ids, &export_local_ids, &export_procs );
  LB2.Evaluate( 1, 0, 0, 0, 0, 0, 0 );

  newObj_ = NewGraph;

  return *NewGraph;
}

} // namespace EpetraExt

