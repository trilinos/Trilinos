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
#include <Isorropia_EpetraInternalPartitioner.hpp>
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
  Operator (input_graph, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, weights, paramlist)
{
  if (compute_partitioning_now)
    partition(true);
}

Partitioner::~Partitioner(){}

void Partitioner::
partition(bool force_repartitioning)
{

  bool use_zoltan = false;
  int input_type = Library::unspecified_input_;
  Teuchos::ParameterList sublist = paramlist_;

  std::string partitioning_method_str("PARTITIONING_METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("ZOLTAN");

  if (alreadyComputed() && !force_repartitioning)
    return;

#ifdef HAVE_ISORROPIA_ZOLTAN
  if (partitioning_method != "SIMPLE_LINEAR") {
    use_zoltan = true;
  }

  if (use_zoltan) {
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
    sublist = (paramlist_.sublist(zoltan));
  }

#else /* HAVE_ISORROPIA_ZOLTAN */
  if (paramlist_.isSublist(zoltan)) {
    throw Isorropia::Exception("Zoltan requested, but Zoltan not enabled.");
  }
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (use_zoltan == false) {
    if (input_matrix_.get() != 0)
      lib_ = Teuchos::rcp(new InternalPartitioner(input_matrix_, costs_));
    else if (input_graph_.get() != 0)
      lib_ = Teuchos::rcp(new InternalPartitioner(input_graph_, costs_));
    else if (input_coords_.get() != 0)
      if (weights_.get()){
        lib_ = Teuchos::rcp(new InternalPartitioner(input_coords_, weights_));
      }
      else{
        lib_ = Teuchos::rcp(new InternalPartitioner(input_coords_));
      }
    else{
      throw Isorropia::Exception("Partitioner::partition - no input object.");
    }
  }

  if (input_coords_.get() != 0){
    input_type = Library::geometric_input_;
  }
  else{
    std::string lb_method_str("LB_METHOD");
    if (sublist.isParameter(lb_method_str)){
      std::string lb_meth = sublist.get(lb_method_str, "HYPERGRAPH");
      if (lb_meth == "GRAPH"){
        input_type = Library::graph_input_;
      }
      else if (lb_meth == "HYPERGRAPH"){
        input_type = Library::hgraph_input_;
      }
      else{
        throw Isorropia::Exception("Valid LB_METHOD parameters for input type are HYPERGRAPH or GRAPH");
      }
    }
    else{
      input_type = Library::hgraph_input_;
    }
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

bool Partitioner::partitioning_already_computed() const {
  return (alreadyComputed());
}

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
  input_map_->MyGlobalElements( &elementList[0] );

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

  Teuchos::RCP<Epetra_Map> target_map =
    Teuchos::rcp(new Epetra_Map(-1, myNewGID.size(), &myNewGID[0], 0, input_map_->Comm()));

  return(target_map);
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

