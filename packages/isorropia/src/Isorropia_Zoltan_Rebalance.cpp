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

#ifdef HAVE_ISORROPIA_ZOLTAN

#ifndef HAVE_MPI
#error "Isorropia_Zoltan requires MPI."
#endif

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>
#include <Isorropia_Rebalance.hpp>

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

#ifdef HAVE_EPETRAEXT
#include <EpetraExt_Transpose_CrsGraph.h>
#endif

#include <Isorropia_ZoltanQuery.h>
#include <IZoltan_LoadBalance.h>

namespace Isorropia_Zoltan {

Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     Teuchos::ParameterList& paramlist)
{
  Teuchos::RefCountPtr<Epetra_Map> bal_map =
    Isorropia_Zoltan::create_balanced_map(input_graph, paramlist);

  Teuchos::RefCountPtr<Epetra_CrsGraph> bal_graph =
    Isorropia::redistribute_rows(input_graph, *bal_map);

  bal_graph->FillComplete();

  return(bal_graph);
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
		     Teuchos::ParameterList& paramlist)
{
  Teuchos::RefCountPtr<Epetra_Map> bal_map =
    Isorropia_Zoltan::create_balanced_map(input_matrix.Graph(), paramlist);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> bal_matrix =
    Isorropia::redistribute_rows(input_matrix, *bal_map);

  bal_matrix->FillComplete();

  return(bal_matrix);
}

Teuchos::RefCountPtr<Epetra_Map>
create_balanced_map(const Epetra_CrsGraph& input_graph,
		    Teuchos::ParameterList& paramlist)
{
  Teuchos::RefCountPtr<Epetra_Map> bal_map;

  Epetra_CrsGraph& nonconst_input = const_cast<Epetra_CrsGraph&>(input_graph);

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan::LoadBalance LB( 0, &dummy, &version );
  int err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(nonconst_input.Comm()).Comm() );

  std::string lb_method_str("LB_METHOD");

  //check for the user-specified value of LB_METHOD, and use "GRAPH" if
  //none is specified.
  std::string lb_meth = paramlist.get(lb_method_str, "GRAPH");

  if ( err == ZOLTAN_OK ) {
    err = LB.Set_Param( lb_method_str.c_str(), lb_meth.c_str() );
  }

  if (lb_meth == "GRAPH" || lb_meth == "PARMETIS") {
    std::string par_method_str("PARMETIS_METHOD");
    std::string par_method = paramlist.get(par_method_str, "PartKway");

    if( err == ZOLTAN_OK ) err = LB.Set_Param( "PARMETIS_METHOD", par_method );
  }

  if (lb_meth == "HYPERGRAPH") {
    //tell the load-balance object to register the hypergraph query functions
    //instead of the regular graph query functions.
    LB.Set_Hypergraph();
  }

  //Setup Query Object
  EpetraExt::CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( nonconst_input );
  Isorropia::ZoltanQuery Query( nonconst_input, &TransGraph );
  if( err == ZOLTAN_OK ) {
    err = LB.Set_QueryObject( &Query );
  }
  else {
    cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; return(bal_map);
  }

  //Generate Load Balance
  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;
  int * import_procs, * export_procs;

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
  for( int i = 0; i < num_export; ++i ) {
    gidSet.insert( export_global_ids[i] );
  }

  //Add unmoved indices to new list
  int loc = 0;
  for( int i = 0; i < numMyElements; ++i ) {
    if( !gidSet.count( elementList[i] ) ) {
      newElementList[loc++] = elementList[i];
    }
  }
  
  //Add imports to end of list
  for( int i = 0; i < num_import; ++i ) {
    newElementList[loc+i] = import_global_ids[i];
  }

  //Free Zoltan Data
  if( err == ZOLTAN_OK ) {
    err = LB.Free_Data( &import_global_ids, &import_local_ids, &import_procs,
                         &export_global_ids, &export_local_ids, &export_procs );
  }

  //Create Import Map
  bal_map =
    Teuchos::rcp(new Epetra_Map( nonconst_input.RowMap().NumGlobalElements(),
				 newNumMyElements, &newElementList[0],
				 nonconst_input.RowMap().IndexBase(),
				 nonconst_input.RowMap().Comm() ));
  return( bal_map );
}

}//namespace Isorropia_Zoltan

#endif

