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

#include <Isorropia_EpetraPartitioner.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_EpetraZoltanLib.hpp>
#endif
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::~Partitioner(){}

void Partitioner::
clearPartSizes()
{
  if (partGIDs){
    delete [] partGIDs;
    partGIDs = NULL;
  }
  if (partSizes){
    delete [] partSizes;
    partSizes = NULL;
  }
  numPartSizes = 0;
}

void Partitioner::
setPartSizes(int len, int *global_part_id, float *part_size)
{
  clearPartSizes();

  if (len < 1) return;

  numPartSizes = len;

  partGIDs = new int [len];
  memcpy(partGIDs, global_part_id, sizeof(int) * len);

  partSizes = new float [len];
  memcpy(partSizes, part_size, sizeof(float) * len);
}

void Partitioner::
partition(bool force_repartitioning)
{

  int input_type = Library::unspecified_input_;

  std::string partitioning_method_str("PARTITIONING METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("ZOLTAN");

  if (alreadyComputed() && !force_repartitioning)
    return;

#ifdef HAVE_ISORROPIA_ZOLTAN
  if (partitioning_method == "SIMPLE_LINEAR") {
    throw Isorropia::Exception("Partitioner::partition - Only Zoltan Partitionner is now supported.");
  }

  if (input_graph_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, costs_));
  else if (input_matrix_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_));
  else if (input_coords_.get() != 0){
    if (weights_.get()){
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_, weights_));
    }
    else{
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_));
    }
  }
  else{
    throw Isorropia::Exception("Partitioner::partition - no input object.");
  }

  lib_->numPartSizes = numPartSizes;
  lib_->partGIDs = partGIDs;
  lib_->partSizes = partSizes;

#endif /* HAVE_ISORROPIA_ZOLTAN */
  Teuchos::ParameterList sublist = paramlist_.sublist(zoltan);
  // TODO: Add "block" and "random" partitioning.

  if (input_coords_.get() != 0){
    if (partitioning_method == "UNSPECIFIED")
      sublist.set("LB_METHOD", "RCB");
    else
      sublist.set("LB_METHOD", partitioning_method);
    input_type = Library::geometric_input_;
  }
  else{
    if (partitioning_method == "GRAPH"){
      input_type = Library::graph_input_;
      sublist.set("LB_METHOD", "GRAPH");
    }
    else // if (lb_meth == "HYPERGRAPH")  // Hypergraph by default
      {
      input_type = Library::hgraph_input_;
      sublist.set("LB_METHOD", "HYPERGRAPH");
    }
  }


  if (paramlist_.isParameter("NUM PARTS")) {
    sublist.set("NUM_GLOBAL_PARTS", paramlist_.get<std::string>("NUM PARTS"));
  }
  if (paramlist_.isParameter("IMBALANCE TOL")) {
    sublist.set("IMBALANCE_TOL", paramlist_.get<std::string>("IMBALANCE TOL"));
  }
  if (paramlist_.isParameter("BALANCE OBJECTIVE")
      && paramlist_.get<std::string>("BALANCE OBJECTIVE") == "NONZEROS") {
    sublist.set("ADD_OBJ_WEIGHT", "NONZEROS");
  }

  lib_->input_type_ = input_type;
  lib_->repartition(sublist, properties_, exportsSize_, imports_);
  computeNumberOfProperties();
  operation_already_computed_ = true;
}

void Partitioner::
compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}

// Create a new RowMap 
Teuchos::RCP<Epetra_Map>
Partitioner::createNewMap()
{
  if (!alreadyComputed()) {
    partition();
  }

  //Generate New Element List
  int myPID = input_map_->Comm().MyPID();
  int numMyElements = input_map_->NumMyElements();
  std::vector<int> elementList( numMyElements );
  if (numMyElements > 0)
    input_map_->MyGlobalElements( &elementList[0] );
  else
    input_map_->MyGlobalElements(NULL);

  std::vector<int> myNewGID (numMyElements - exportsSize_);
  std::vector<int>::iterator newElemsIter;
  std::vector<int>::const_iterator elemsIter;

  for (elemsIter = properties_.begin(), newElemsIter= myNewGID.begin() ;
       elemsIter != properties_.end() ; elemsIter ++) {
    if ((*elemsIter) == myPID) {
      (*newElemsIter) = elementList[elemsIter - properties_.begin()];
      newElemsIter ++;
    }
  }
  //Add imports to end of list
  myNewGID.insert(myNewGID.end(), imports_.begin(), imports_.end());

  int *gidptr;
  if (myNewGID.size() > 0)
    gidptr = &myNewGID[0];
  else
    gidptr = NULL;

  Teuchos::RCP<Epetra_Map> target_map =
    Teuchos::rcp(new Epetra_Map(-1, myNewGID.size(), gidptr, 0, input_map_->Comm()));

  return(target_map);
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

