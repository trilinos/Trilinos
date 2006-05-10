//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_Zoltan_Rebalance.hpp>

#ifdef HAVE_EPETRAEXT_ZOLTAN

#ifndef HAVE_MPI
#error "Isorropia_Zoltan requires MPI."
#endif

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>

#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#include <EpetraExt_Transpose_CrsGraph.h>

#include <EpetraExt_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

namespace Isorropia_Zoltan {

Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     Teuchos::ParameterList& paramlist)
{
  std::string bal_alg_str("Balancing algorithm");
  std::string bal_alg = paramlist.get(bal_alg_str, "none_specified");
  if (bal_alg == "hypergraph") {
    throw Isorropia::Exception("hypergraph partitioning not yet supported.");
  }

  std::string part_method_str("PARMETIS_METHOD");
  std::string part_method = paramlist.get(part_method_str, "PartKway");

  int err;

  Epetra_CrsGraph& nonconst_input = const_cast<Epetra_CrsGraph&>(input_graph);

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan::LoadBalance LB( 0, &dummy, &version );
  err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(nonconst_input.Comm()).Comm() );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "LB_METHOD", "PARMETIS" );
  if( err == ZOLTAN_OK ) err = LB.Set_Param( "PARMETIS_METHOD", part_method );

  //Setup Query Object
  EpetraExt::CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( nonconst_input );
  EpetraExt::ZoltanQuery Query( nonconst_input, &TransGraph );
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

  nonconst_input.Comm().Barrier();
  err = LB.Balance( &changes,
                     &num_gid_entries, &num_lid_entries,
                     &num_import, &import_global_ids, &import_local_ids, &import_procs,
                     &num_export, &export_global_ids, &export_local_ids, &export_procs );
  LB.Evaluate( 1, 0, 0, 0, 0, 0, 0 );
  nonconst_input.Comm().Barrier();

  //Generate New Element List
  int numMyElements = nonconst_input.RowMap().NumMyElements();
  std::vector<int> elementList( numMyElements );
  nonconst_input.RowMap().MyGlobalElements( &elementList[0] );

  int newNumMyElements = numMyElements - num_export + num_import;
  std::vector<int> newElementList( newNumMyElements );

  std::set<int> gidSet;
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
  Epetra_Map* NewRowMap_ = new Epetra_Map( nonconst_input.RowMap().NumGlobalElements(),
                                           newNumMyElements,
                                           &newElementList[0],
                                           nonconst_input.RowMap().IndexBase(),
                                           nonconst_input.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *NewRowMap_, nonconst_input.RowMap() );

  //Create New Graph
  Teuchos::RefCountPtr<Epetra_CrsGraph> bal_graph =
    Teuchos::rcp(new Epetra_CrsGraph( Copy, *NewRowMap_, 0 ));

  bal_graph->Import( nonconst_input, Importer, Insert );
  bal_graph->FillComplete();

//   Zoltan::LoadBalance LB2( 0, &dummy, &version );
//   err = LB2.Create( dynamic_cast<const Epetra_MpiComm&>(orig.Comm()).Comm() );
//   if( err == ZOLTAN_OK ) err = LB2.Set_Param( "LB_METHOD", "PARMETIS" );
//   if( err == ZOLTAN_OK ) err = LB2.Set_Param( "PARMETIS_METHOD", partitionMethod_ );
//   CrsGraph_Transpose transTrans;
//   Epetra_CrsGraph & trans2 = transTrans( *bal_graph );
//   ZoltanQuery query( *bal_graph, &trans2 );
//   if( err == ZOLTAN_OK ) err = LB2.Set_QueryObject( &query );
//   err = LB2.Balance( &changes,
//                      &num_gid_entries, &num_lid_entries,
//                      &num_import, &import_global_ids, &import_local_ids, &import_procs,
//                      &num_export, &export_global_ids, &export_local_ids, &export_procs );
//   LB2.Evaluate( 1, 0, 0, 0, 0, 0, 0 );

  return(bal_graph);
}

}//namespace Isorropia_Zoltan

#endif

