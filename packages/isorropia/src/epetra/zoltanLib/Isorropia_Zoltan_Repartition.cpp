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
#include <Epetra_RowMatrix.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_CrsGraph.h>
#endif

#include <QueryObject.hpp>

namespace Isorropia {

namespace Epetra {

namespace ZoltanLib {

static int me = 0;

//TODO do we need to add these functions for geometric partitioning? 
// Are these being used?

int
repartition(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
	    Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
            Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports)
{
  int inputType = QueryObject::hgraph_input_;
  std::string lb_method_str("LB_METHOD");
  if (paramlist.isParameter(lb_method_str)){
    std::string lb_meth = paramlist.get(lb_method_str, "HYPERGRAPH");
    if (lb_meth == "GRAPH"){
      inputType = QueryObject::graph_input_;
    }
  }

  Teuchos::RefCountPtr<QueryObject> queryObject = 
      Teuchos::rcp(new QueryObject(input_graph, costs, inputType));

  const Epetra_Comm &ecomm = input_graph->RowMap().Comm();
#ifdef HAVE_MPI
  const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);

  MPI_Comm mpicomm = empicomm.Comm();

  return( load_balance(mpicomm, paramlist, *queryObject,
		       myNewElements, exports, imports) );
#else
  return (1);
#endif
}

int
repartition(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
	    Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
            Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports)
{
  int inputType = QueryObject::hgraph_input_;
  std::string lb_method_str("LB_METHOD");
  if (paramlist.isParameter(lb_method_str)){
    std::string lb_meth = paramlist.get(lb_method_str, "HYPERGRAPH");
    if (lb_meth == "GRAPH"){
      inputType = QueryObject::graph_input_;
    }
  }

  Teuchos::RefCountPtr<QueryObject> queryObject = 
    Teuchos::rcp(new QueryObject(input_matrix, costs, inputType));

  const Epetra_Comm &ecomm = input_matrix->RowMatrixRowMap().Comm();
#ifdef HAVE_MPI
  const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);

  MPI_Comm mpicomm = empicomm.Comm();

  return( load_balance(mpicomm, paramlist, *queryObject,
		       myNewElements, exports, imports) );
#else
  return (1);
#endif
}

#ifdef HAVE_MPI
int
load_balance(MPI_Comm &comm,
	     Teuchos::ParameterList& paramlist,
	     QueryObject& queryObject,
	     std::vector<int>& myNewElements,
	     std::map<int,int>& exports,
	     std::map<int,int>& imports)
{
  float version;
  int argcTmp=0;
  char *argvTmp[1];

  // create a Zoltan problem

  argvTmp[0] = NULL;
  Zoltan_Initialize(argcTmp, argvTmp, &version);

  Zoltan *zz = new Zoltan(comm);

  if (zz == NULL){
    throw Isorropia::Exception("Error creating Zoltan object");
    return -1;    
  }

  // set problem parameters

  std::string dbg_level_str("DEBUG_LEVEL");
  if (!paramlist.isParameter(dbg_level_str)) {
    paramlist.set(dbg_level_str, "0");
  }

  std::string lb_approach_str("LB_APPROACH");
  if (!paramlist.isParameter(lb_approach_str)) {
    paramlist.set(lb_approach_str, "PARTITION");
  }

  std::string lb_method_str("LB_METHOD");
  std::string lb_meth = paramlist.get(lb_method_str, "HYPERGRAPH");

  if (!paramlist.isParameter(lb_method_str)) {
    paramlist.set(lb_method_str, lb_meth);  // set to HYPERGRAPH
  }

  // If no process set vertex weights, let Zoltan default to
  // unit weights for vertices

  if (queryObject.haveVertexWeights()) {
    if (!paramlist.isParameter("OBJ_WEIGHT_DIM")) {
      paramlist.set("OBJ_WEIGHT_DIM", "1");
    }
  }

  // If no process set graph or hypergraph edge weights, 
  // let Zoltan default to unit weights for edges

  if (queryObject.haveGraphEdgeWeights() ||
      queryObject.haveHypergraphEdgeWeights()) {
    if (!paramlist.isParameter("EDGE_WEIGHT_DIM")) {
      paramlist.set("EDGE_WEIGHT_DIM", "1");
    }
  }

  Teuchos::ParameterList::ConstIterator
    iter = paramlist.begin(),
    iter_end = paramlist.end();

  for(; iter != iter_end; ++iter) {
    const std::string& name = iter->first;
    const std::string& value = Teuchos::getValue<std::string>(iter->second);
    zz->Set_Param(name, value);
  }

  // Set the query functions
  // M.M.W. do we need to support hierarchical partitioning here?

  zz->Set_Num_Obj_Fn(QueryObject::Number_Objects, (void *)&queryObject);
  zz->Set_Obj_List_Fn(QueryObject::Object_List, (void *)&queryObject);

  if (lb_meth == "HYPERGRAPH"){
    zz->Set_HG_Size_CS_Fn(QueryObject::HG_Size_CS, (void *)&queryObject);
    zz->Set_HG_CS_Fn(QueryObject::HG_CS, (void *)&queryObject);
    zz->Set_HG_Size_Edge_Wts_Fn(QueryObject::HG_Size_Edge_Weights , 
                                (void *)&queryObject);
    zz->Set_HG_Edge_Wts_Fn(QueryObject::HG_Edge_Weights, (void *)&queryObject);
  }
  else{
    zz->Set_Num_Edges_Multi_Fn(QueryObject::Number_Edges_Multi, (void *)&queryObject);
    zz->Set_Edge_List_Multi_Fn(QueryObject::Edge_List_Multi, (void *)&queryObject);
  }



  //Generate Load Balance
  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;
  int * import_procs, * export_procs;
  int *import_to_part, *export_to_part;

  int err = zz->LB_Partition(changes, num_gid_entries, num_lid_entries,
   num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
   num_export, export_global_ids, export_local_ids, export_procs, export_to_part );

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing partitioning with Zoltan");
    return -1;
  }

  //Generate New Element List
  int numMyElements = queryObject.RowMap().NumMyElements();
  std::vector<int> elementList( numMyElements );
  queryObject.RowMap().MyGlobalElements( &elementList[0] );

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
  zz->LB_Free_Part(&import_global_ids, &import_local_ids, 
                     &import_procs, &import_to_part);
  zz->LB_Free_Part(&export_global_ids, &export_local_ids, 
                     &export_procs, &export_to_part);

  delete zz;

  return( 0 );
}
#endif

}//namespace ZoltanLib
}//namespace Epetra
}//namespace Isorropia

#endif

