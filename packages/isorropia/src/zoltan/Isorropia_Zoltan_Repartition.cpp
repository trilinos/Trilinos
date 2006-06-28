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

void
set_zoltan_parameters(Zoltan::LoadBalance& LB,
		      const Teuchos::ParameterList& paramlist)
{
  Teuchos::ParameterList::ConstIterator
    iter = paramlist.begin(),
    iter_end = paramlist.end();

  for(; iter != iter_end; ++iter) {
    const std::string& name = iter->first;
    const Teuchos::ParameterEntry& param = iter->second;

    if (param.isType<std::string>()) {
      const std::string& value = Teuchos::getValue<std::string>(param);
      LB.Set_Param(name, value);
    }
    else if (param.isType<double>()) {
      double value = Teuchos::getValue<double>(param);
      std::ostringstream osstr;
      osstr << value;
      std::string str = osstr.str();
      LB.Set_Param(name, str);
    }
    else if (param.isType<int>()) {
      int value = Teuchos::getValue<int>(param);
      std::ostringstream osstr;
      osstr << value;
      std::string str = osstr.str();
      LB.Set_Param(name, str);
    }
  }
}

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

  //if LB_METHOD has not been specified, then set it to GRAPH.
  std::string lb_method_str("LB_METHOD");
  if (!paramlist.isParameter(lb_method_str)) {
    paramlist.set(lb_method_str, "GRAPH");
  }

  //check for the value of LB_METHOD (using "GRAPH" by default)
  std::string lb_meth = paramlist.get(lb_method_str, "GRAPH");

  if (lb_meth == "HYPERGRAPH") {
    //tell the load-balance object to register the hypergraph query functions
    //instead of the regular graph query functions.
    LB.Set_Hypergraph();
  }

  //go ahead and pass the parameters in paramlist to Zoltan
  set_zoltan_parameters(LB, paramlist);

  //Setup Query Object
  EpetraExt::CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( nonconst_input );

  Zoltan::QueryObject* queryObject =
    new Isorropia::ZoltanQuery(nonconst_input, &TransGraph);

  if( err == ZOLTAN_OK ) {
    err = LB.Set_QueryObject( queryObject );
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

  delete queryObject;

  return( 0 );
}

}//namespace Isorropia_Zoltan

#endif

