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

#include <Isorropia_Zoltan_Repartition.hpp>

#ifdef HAVE_ISORROPIA_ZOLTAN

#ifndef HAVE_MPI
#error "Isorropia_Zoltan requires MPI."
#endif

#include <Isorropia_Exception.hpp>

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
#include <EpetraExt_Transpose_RowMatrix.h>
#endif

#include <Isorropia_ZoltanQuery.h>
#include <IZoltan_LoadBalance.h>

namespace Isorropia_Zoltan {

int
repartition(const Epetra_CrsGraph& input_graph,
            Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports)
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
    cout << "Setup of Zoltan Load Balancing Objects FAILED!\n";
    return(err);
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
  nonconst_input.Comm().Barrier();

  //Generate New Element List
  int numMyElements = nonconst_input.RowMap().NumMyElements();
  std::vector<int> elementList( numMyElements );
  nonconst_input.RowMap().MyGlobalElements( &elementList[0] );

  int newNumMyElements = numMyElements - num_export + num_import;
  myNewElements.resize( newNumMyElements );

  for( int i = 0; i < num_export; ++i ) {
    exports[export_global_ids[i]] = export_procs[i];
  }

  for( int i = 0; i < num_import; ++i ) {
    imports[import_global_ids[i]] = import_procs[i];
  }

  //Add unmoved indices to new list
  int loc = 0;
  for( int i = 0; i < numMyElements; ++i ) {
    if( !exports.count( elementList[i] ) ) {
      myNewElements[loc++] = elementList[i];
    }
  }
  
  //Add imports to end of list
  for( int i = 0; i < num_import; ++i ) {
    myNewElements[loc+i] = import_global_ids[i];
  }

  //Free Zoltan Data
  if( err == ZOLTAN_OK ) {
    err = LB.Free_Data( &import_global_ids, &import_local_ids, &import_procs,
                         &export_global_ids, &export_local_ids, &export_procs );
  }

  return( 0 );
}

}//namespace Isorropia_Zoltan

#endif

