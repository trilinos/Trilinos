//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//************************************************************************
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 CostDescriber *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<CostDescriber>(costs,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 CostDescriber *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<CostDescriber>(costs,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_BlockMap> input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_map, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_BlockMap *input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_BlockMap>(input_map,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, coords, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
                         const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, coords, weights, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 CostDescriber *costs,
                         const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<CostDescriber>(costs,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
                         const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 CostDescriber *costs,
                         const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<CostDescriber>(costs,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

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

  if (input_graph_.get() != 0 && input_coords_.get() != 0)
  {
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_,costs_,input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_,input_coords_));
    }
  }
  else if (input_matrix_.get() != 0 && input_coords_.get() != 0)
  {
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_,costs_,input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_,input_coords_));
    }
  }
  else if (input_graph_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, costs_));
  else if (input_matrix_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_));
  else if (input_coords_.get() != 0)
  {
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_));
    }
  }
  else if (input_map_.get() != 0)
  {
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_map_));
  }
  else
  {
    throw Isorropia::Exception("Partitioner::partition - no input object.");
  }

  lib_->numPartSizes = numPartSizes;
  lib_->partGIDs = partGIDs;
  lib_->partSizes = partSizes;

#endif /* HAVE_ISORROPIA_ZOLTAN */
  Teuchos::ParameterList sublist = paramlist_.sublist(zoltan);
  // TODO: Add "block" and "random" partitioning.

  if (partitioning_method == "UNSPECIFIED" && sublist.isParameter("LB_METHOD")) 
  {
    throw Isorropia::Exception("Isorropia \"PARTITIONING METHOD\" as to be set\n"
			       "ZOLTAN/LB_METHOD is no longer supported.\n"
                               "See readme and release notes for details.");
  }

  if (input_coords_.get() != 0)
  {
    if (partitioning_method == "UNSPECIFIED")
    {
      sublist.set("LB_METHOD", "RCB");
      input_type = Library::geometric_input_;
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::graph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_graph_geometric_input_;
    }
    else // Default case: copy partitioning_method (e.g., RCB, RIB, HSFC)
    {
      sublist.set("LB_METHOD", partitioning_method);
      input_type = Library::geometric_input_;
    }
  }
  else if (input_graph_.get() != 0 || input_matrix_.get() != 0) // graph or matrix input
  {
    if (partitioning_method == "GRAPH")
    {
      input_type = Library::graph_input_;
      sublist.set("LB_METHOD", "GRAPH");
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_GRAPH")
    {
      input_type = Library::graph_input_;
      sublist.set("LB_METHOD", "HIER");
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH") // Can perhaps simplify this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_graph_input_;
    }
    
    else //Hypergraph by default
    {
      input_type = Library::hgraph_input_;
      sublist.set("LB_METHOD", "HYPERGRAPH");
    }
  }
  else  // BlockMap input
  {
    if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else // BLOCK by default
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
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

  std::string print_metrics_str("PRINT ZOLTAN METRICS");
  std::string print_metrics =
    paramlist_.get(print_metrics_str, "UNSPECIFIED");

  if (print_metrics == ("UNSPECIFIED")){
    printMetrics = 0;
  }
  else if (print_metrics == ("1")){
    printMetrics = 1;
  }
  else if (print_metrics == ("2")){
    printMetrics = 2;
  }
  else {
    printMetrics = 1;
  }

  lib_->input_type_ = input_type;
  lib_->repartition(sublist, properties_, exportsSize_, imports_);
  computeNumberOfProperties();
  operation_already_computed_ = true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Partitioner::
compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<Epetra_Map>
Partitioner::createNewMap()
{
  Epetra_Map *outputMap;

  createNewMap(outputMap);

  return( Teuchos::RCP<Epetra_Map>(outputMap) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
void
Partitioner::createNewMap(Epetra_Map * &outputMap)
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
    input_map_->MyGlobalElements((int*)NULL); // disambiguate int/long long

  int newGIDSize = numMyElements - exportsSize_;

  std::vector<int> myNewGID;

  if (newGIDSize > 0){
    myNewGID.resize(newGIDSize);
    std::vector<int>::iterator newElemsIter;
    std::vector<int>::const_iterator elemsIter;

    for (elemsIter = properties_.begin(), newElemsIter= myNewGID.begin() ;
         elemsIter != properties_.end() ; elemsIter ++) {
      if ((*elemsIter) == myPID) {
        (*newElemsIter) = elementList[elemsIter - properties_.begin()];
        newElemsIter ++;
      }
    }
  }
  //Add imports to end of list
  myNewGID.insert(myNewGID.end(), imports_.begin(), imports_.end());

  int *gidptr;
  if (myNewGID.size() > 0)
    gidptr = &myNewGID[0];
  else
    gidptr = NULL;

  outputMap = new Epetra_Map(-1, myNewGID.size(), gidptr, 0, input_map_->Comm());

  return;
}
////////////////////////////////////////////////////////////////////////////////


} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

