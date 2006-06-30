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

#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_Zoltan_Repartition.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_CostDescriber.hpp>

#include <Teuchos_RefCountPtr.hpp>

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

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA

Epetra::Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Epetra::Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
		  Teuchos::RefCountPtr<const Isorropia::CostDescriber> costs,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(costs),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Epetra::Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Epetra::Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
		  Teuchos::RefCountPtr<const Isorropia::CostDescriber> costs,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(costs),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Epetra::Partitioner::~Partitioner()
{
}

void Epetra::Partitioner::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

void Epetra::Partitioner::compute_partitioning(bool force_repartitioning)
{
  if (partitioning_already_computed_) {
    if (!force_repartitioning) {
      return;
    }
  }

  if (input_graph_.get() == 0 && input_matrix_.get() == 0) {
    std::string str1("Isorropia::Epetra::Partitioner::compute_partitioning ERROR: ");
    std::string str2("not holding valid input graph or matrix.");
    throw Isorropia::Exception(str1+str2);
  }

  int err = 0;

  std::string zoltan("Zoltan");
  if (paramlist_.isSublist(zoltan)) {
#ifdef HAVE_ISORROPIA_ZOLTAN

    Teuchos::ParameterList& sublist = paramlist_.sublist(zoltan);

    if (input_graph_.get() == 0) {
      err = Isorropia_Zoltan::repartition(input_matrix_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }
    else {
      err = Isorropia_Zoltan::repartition(input_graph_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }
#else
    throw Isorropia::Exception("Zoltan requested, but zoltan not enabled.");
#endif
  }
  else {
    if (input_graph_.get() != 0) {
      weights_ = Teuchos::rcp(Epetra::create_row_weights_nnz(*input_graph_));
    }
    else {
      weights_ = Teuchos::rcp(Epetra::create_row_weights_nnz(*input_matrix_));
    }

    err = Isorropia::Epetra::repartition(*input_map_,
                                               *weights_,
                                               myNewElements_,
                                               exports_, imports_);
  }

  if (err != 0) {
    throw Isorropia::Exception("error in repartitioning");
  }

  partitioning_already_computed_ = true;
}

bool Epetra::Partitioner::partitioning_already_computed() const
{
  return partitioning_already_computed_;
}

int Epetra::Partitioner::newPartitionNumber(int myElem) const
{
  std::map<int,int>::const_iterator iter = exports_.find(myElem);
  if (iter != exports_.end()) {
    return(iter->second);
  }

  return( input_graph_->RowMap().Comm().MyPID() );
}

int Epetra::Partitioner::numElemsInPartition(int partition) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("Epetra::Partitioner::numElemsInPartition not implemented for non-local partitions.");
  }

  return(myNewElements_.size());
}

void
Epetra::Partitioner::elemsInPartition(int partition, int* elementList, int len) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("error in Epetra_Map::MyGlobalElements");
  }

  unsigned length = len;
  if (myNewElements_.size() < length) length = myNewElements_.size();

  for(unsigned i=0; i<length; ++i) {
    elementList[i] = myNewElements_[i];
  }
}

#endif //HAVE_EPETRA

}//namespace Isorropia

