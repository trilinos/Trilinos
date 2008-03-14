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
#include <Isorropia_EpetraCostDescriber.hpp>

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

namespace Epetra {

Partitioner::
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

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
		  Teuchos::RefCountPtr<CostDescriber> costs,
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

Partitioner::
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

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
		  Teuchos::RefCountPtr<CostDescriber> costs,
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

Partitioner::~Partitioner()
{
}

void Partitioner::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

void Partitioner::compute_partitioning(bool force_repartitioning)
{
  int err = 0;
  std::string str1("Isorropia::Partitioner::compute_partitioning ");
  std::string str2;

  if (partitioning_already_computed_) {
    if (!force_repartitioning) {
      return;
    }
  }

  if (input_graph_.get() == 0 && input_matrix_.get() == 0) {
    str2 = "ERROR: not holding valid input graph OR matrix.";
    throw Isorropia::Exception(str1+str2);
  }

  const Epetra_Comm &comm = input_map_->Comm();
  int localProc = comm.MyPID();

  //if Isorropia was configured with Zoltan support, then we will use
  //Zoltan unless the user specified "PARTITIONING_METHOD" = "SIMPLE_LINEAR".

  bool use_zoltan = false;

#ifdef HAVE_ISORROPIA_ZOLTAN
  std::string partitioning_method_str("PARTITIONING_METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");
  if (partitioning_method != "SIMPLE_LINEAR") {
    use_zoltan = true;
  }

  if (use_zoltan){

    // Check that we have a valid Zoltan problem.
    //
    // Zoltan parameter LB_METHOD:
    //
    //   GRAPH:
    //     matrix must be square and symmetric
    //     if the graph has self edges (nonzeros in diagonal) we 
    //       don't flag this as an error, but we will omit them 
    //       in the Zoltan queries
    //
    //   HYPERGRAPH:
    //     matrix can be any shape, need not be symmetric
    //     this is the default if LB_METHOD is not set
     
  
    std::string lb_method_str("LB_METHOD");
    if (paramlist_.isParameter(lb_method_str)){
      std::string lb_meth = paramlist_.get(lb_method_str, "HYPERGRAPH");
      if (lb_meth == "GRAPH"){
        bool square = false;
        bool symmetric = false;
        if (input_graph_.get() != 0){
          if (input_graph_->NumGlobalRows() == input_graph_->NumGlobalCols()){
            square = true;
            // TODO - is there a quick way to figure out if the graph is
            // symmetric?  I can't see a way to do it.  For now we let
            // Zoltan figure this out.
            symmetric = true;
          }
        }
        else{
          if (input_matrix_->NumGlobalRows() == input_matrix_->NumGlobalCols()){
            square = true;
            // TODO - is there a quick way to figure out if the matrix is
            // symmetric?  I can't see a way to do it.  For now we let
            // Zoltan figure this out.
            symmetric = true;
          }
        } 
        if (!square){
          str2 = "LB_METHOD=GRAPH, matrix or graph must be square";
          throw Isorropia::Exception(str1+str2);
        }
        if (!symmetric){
          str2 = "LB_METHOD=GRAPH, matrix or graph must be symmetric";
          throw Isorropia::Exception(str1+str2);
        }
      }
    }
    // If vertex/edge costs have been set, do a global operation to find
    // out how many weights were given.  Some processes may provide no
    // weights - they need to be informed that weights are being provided
    // by the application.  Do some sanity checks.
  
    err = 0;
    int gerr = 0;
    int base = input_map_->IndexBase();
  
    int numMyVWeights = 0;
    int numMyGWeights = 0;
    int numMyHGWeights = 0;
    int globalNumCols = 0;
    int myNZ = 0;
    int globalNZ = 0;
    int mySelfEdges = 0;
    int globalSelfEdges = 0;

    int myRows = input_map_->NumMyElements();
    int globalNumRows = input_map_->NumGlobalElements();

    if (input_graph_.get() == 0) {
      myNZ = input_matrix_->NumMyNonzeros();
      mySelfEdges = input_matrix_->NumMyDiagonals();
      globalNZ = input_matrix_->NumGlobalNonzeros();
      globalSelfEdges = input_matrix_->NumGlobalDiagonals();
      globalNumCols = input_matrix_->NumGlobalCols();
    }
    else{
      myNZ = input_graph_->NumMyNonzeros();
      mySelfEdges = input_graph_->NumMyDiagonals();
      globalNZ = input_graph_->NumGlobalNonzeros();
      globalSelfEdges = input_graph_->NumGlobalDiagonals();
      globalNumCols = input_graph_->NumGlobalCols();
    }
  
    if (costs_.get() != 0) {
  
      numMyVWeights = costs_->getNumVertices();
  
      if (costs_->haveGraphEdgeWeights()){
        for (int i=0; i<numMyVWeights; i++){
          int gid = input_map_->GID(i);
          if (gid >= base){
            numMyGWeights += costs_->getNumGraphEdges(gid);
          }
        }
      }
      numMyHGWeights = costs_->getNumHypergraphEdgeWeights();
  
      if ((numMyVWeights > 0) && (numMyVWeights != myRows)){
        str2 = "Number of my vertex weights != number of my rows";
        err = 1;
      }
      else if ((numMyGWeights > 0) && (numMyGWeights != (myNZ - mySelfEdges))){
        str2 = "Number of my graph edge weights != number of my nonzeros";
        err = 1;
      }
    }
    else{
      costs_ = Teuchos::rcp(new CostDescriber());
    }

    comm.SumAll(&err, &gerr ,1);
    if (gerr > 0){
      throw Isorropia::Exception(str1+str2);
    }
  
    int lval[4], gval[4];
    lval[0] = numMyVWeights;
    lval[1] = numMyGWeights;
    lval[2] = numMyHGWeights;
  
    comm.SumAll(lval, gval, 3);
  
    int numVWeights = gval[0];
    int numGWeights = gval[1];
    int numHGWeights = gval[2];
  
    if ((numVWeights > 0) && (numVWeights != globalNumRows)){
      str2 = "Number of vertex weights supplied by application != number of rows";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numGWeights > 0) && (numGWeights != (globalNZ - globalSelfEdges))){
      str2 = "Number of graph edge weights supplied by application != number of edges";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numHGWeights > 0) && (numHGWeights < globalNumCols)){
      str2 = "Number of hyperedge weights supplied by application < number of columns";
      throw Isorropia::Exception(str1+str2);
    }

    costs_->setNumGlobalVertexWeights(numVWeights);
    costs_->setNumGlobalGraphEdgeWeights(numGWeights);
    costs_->setNumGlobalHypergraphEdgeWeights(numHGWeights);

  }   // end if (use_zoltan)
#endif

  std::string zoltan("Zoltan");
  if (use_zoltan || paramlist_.isSublist(zoltan)) {
#ifdef HAVE_ISORROPIA_ZOLTAN

    Teuchos::ParameterList& sublist = paramlist_.sublist(zoltan);

    if (input_graph_.get() == 0) {
      err = Isorropia::Epetra::ZoltanLib::repartition(input_matrix_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }
    else {
      err = Isorropia::Epetra::ZoltanLib::repartition(input_graph_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }

    if (err != 0) {
      throw Isorropia::Exception("error in Isorropia_Zoltan::repartition");
    }
#else
    throw Isorropia::Exception("Zoltan requested, but Zoltan not enabled.");
#endif
  }
  else { //we'll use the built-in simple-linear partitioner.
    if (input_graph_.get() != 0) {
      weights_ = Teuchos::rcp(create_row_weights_nnz(*input_graph_));
    }
    else {
      weights_ = Teuchos::rcp(create_row_weights_nnz(*input_matrix_));
    }

    err = Isorropia::Epetra::repartition(*input_map_,
                                               *weights_,
                                               myNewElements_,
                                               exports_, imports_);

    if (err != 0) {
      throw Isorropia::Exception("error in simple linear repartitioning");
    }
  }

  partitioning_already_computed_ = true;
}

bool Partitioner::partitioning_already_computed() const
{
  return partitioning_already_computed_;
}

int Partitioner::newPartitionNumber(int myElem) const
{
  std::map<int,int>::const_iterator iter = exports_.find(myElem);
  if (iter != exports_.end()) {
    return(iter->second);
  }

  return( input_graph_->RowMap().Comm().MyPID() );
}

int Partitioner::numElemsInPartition(int partition) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("Partitioner::numElemsInPartition not implemented for non-local partitions.");
  }

  return(myNewElements_.size());
}

void
Partitioner::elemsInPartition(int partition, int* elementList, int len) const
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

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

