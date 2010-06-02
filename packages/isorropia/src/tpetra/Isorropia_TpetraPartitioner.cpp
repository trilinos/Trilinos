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

#include <Isorropia_TpetraPartitioner.hpp>

#ifdef HAVE_ISORROPIA_ZOLTAN

#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Exception.hpp>

#ifdef HAVE_ISORROPIA_TPETRA

#include <Isorropia_TpetraLibrary.hpp>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_ISORROPIA_TPETRA

namespace Tpetra {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::CrsGraph<int,int,Node> *input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> >(input_graph,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
			 Teuchos::RCP<CostDescriber<Node> > costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::CrsGraph<int,int,Node> *input_graph,
			 CostDescriber<Node>  *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> >(input_graph,false), 
            Teuchos::RCP<CostDescriber<Node> >(costs,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::RowMatrix<double,int,int,Node> *input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> >(input_matrix,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
			 Teuchos::RCP<CostDescriber<Node> > costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::RowMatrix<double,int,int,Node> *input_matrix,
			 CostDescriber<Node>  *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> >(input_matrix,false), 
            Teuchos::RCP<CostDescriber<Node> >(costs,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::MultiVector<double,int,int,Node> *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> >(coords,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > coords,
                         Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::MultiVector<double,int,int,Node> *coords,
			       const ::Tpetra::MultiVector<double,int,int,Node> *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> >(coords,false), 
            Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> >(weights,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(Teuchos::RCP<const ::Tpetra::Map<int,int,Node> > input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_map, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class Node>
Partitioner<Node>::Partitioner(const ::Tpetra::Map<int,int,Node> *input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const ::Tpetra::Map<int,int,Node> >(input_map,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0)
{
  if (compute_partitioning_now)
    partition(true);
}
////////////////////////////////////////////////////////////////////////////////


template<class Node>
Partitioner<Node>::~Partitioner(){}

template<class Node>
void Partitioner<Node>::clearPartSizes()
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

template<class Node>
void Partitioner<Node>::setPartSizes(int len, int *global_part_id, float *part_size)
{
  clearPartSizes();

  if (len < 1) return;

  numPartSizes = len;

  partGIDs = new int [len];
  memcpy(partGIDs, global_part_id, sizeof(int) * len);

  partSizes = new float [len];
  memcpy(partSizes, part_size, sizeof(float) * len);
}

template<class Node>
void Partitioner<Node>::partition(bool force_repartitioning)
{

  int input_type = Library<Node>::unspecified_input_;

  std::string partitioning_method_str("PARTITIONING METHOD");
  std::string partitioning_method =
    this->paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("ZOLTAN");

  if (alreadyComputed() && !force_repartitioning)
    return;

#ifdef HAVE_ISORROPIA_ZOLTAN

  if (this->input_graph_.get() != 0 && this->input_coords_.get() != 0)
  {
    if (this->weights_.get())
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_graph_,this->costs_,this->input_coords_, this->weights_));
    }
    else
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_graph_,this->input_coords_));
    }
  }
  else if (this->input_matrix_.get() != 0 && this->input_coords_.get() != 0)
  {
    if (this->weights_.get())
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_matrix_,this->costs_,this->input_coords_, this->weights_));
    }
    else
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_matrix_,this->input_coords_));
    }
  }
  else if (this->input_graph_.get() != 0)
    this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_graph_, this->costs_));
  else if (this->input_matrix_.get() != 0)
    this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_matrix_, this->costs_));
  else if (this->input_coords_.get() != 0)
  {
    if (this->weights_.get())
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_coords_, this->weights_));
    }
    else
    {
      this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_coords_));
    }
  }
  else if (this->input_map_.get() != 0)
  {
    this->lib_ = Teuchos::rcp(new ZoltanLibClass<Node>(this->input_map_));
  }
  else
  {
    throw Isorropia::Exception("Partitioner::partition - no input object.");
  }

  this->lib_->numPartSizes = numPartSizes;
  this->lib_->partGIDs = partGIDs;
  this->lib_->partSizes = partSizes;

#endif /* HAVE_ISORROPIA_ZOLTAN */
  Teuchos::ParameterList sublist = this->paramlist_.sublist(zoltan);

  if (partitioning_method == "UNSPECIFIED" && sublist.isParameter("LB_METHOD")) 
  {
    throw Isorropia::Exception("Isorropia \"PARTITIONING METHOD\" as to be set\n"
			       "ZOLTAN/LB_METHOD is no longer supported.\n"
                               "See readme and release notes for details.");
  }

  if (this->input_coords_.get() != 0)
  {
    if (partitioning_method == "UNSPECIFIED")
    {
      sublist.set("LB_METHOD", "RCB");
      input_type = Library<Node>::geometric_input_;
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library<Node>::graph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library<Node>::hgraph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library<Node>::hgraph_graph_geometric_input_;
    }
    else
    {
      sublist.set("LB_METHOD", partitioning_method);
      input_type = Library<Node>::geometric_input_;
    }
  }
  else if (this->input_graph_.get() != 0 || this->input_matrix_.get() != 0) // graph or matrix input
  {
    if (partitioning_method == "GRAPH")
    {
      input_type = Library<Node>::graph_input_;
      sublist.set("LB_METHOD", "GRAPH");
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH") // Can perhaps simplify this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library<Node>::hgraph_graph_input_;
    }
    
    else //Hypergraph by default
    {
      input_type = Library<Node>::hgraph_input_;
      sublist.set("LB_METHOD", "HYPERGRAPH");
    }
  }
  else  // BlockMap input
  {
    if (partitioning_method == "CYCLIC")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else // BLOCK by default
    {
      input_type = Library<Node>::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
  }

  if (this->paramlist_.isParameter("NUM PARTS")) 
  {
    std::string aStr = (this->paramlist_).get("NUM PARTS");
    sublist.set("NUM_GLOBAL_PARTS", aStr);
    //sublist.set("NUM_GLOBAL_PARTS", this->paramlist_.get<std::string>("NUM PARTS"));
  }
  if (this->paramlist_.isParameter("IMBALANCE TOL")) 
  {
    std::string aStr = this->paramlist_.get("IMBALANCE TOL");
    sublist.set("IMBALANCE_TOL", aStr);
  }
  //if (this->paramlist_.isParameter("BALANCE OBJECTIVE")
  //    && this->paramlist_.get<std::string>("BALANCE OBJECTIVE") == "NONZEROS") 

  if (this->paramlist_.isParameter("BALANCE OBJECTIVE")
      && this->paramlist_.get("BALANCE OBJECTIVE") == "NONZEROS") 
  {
    sublist.set("ADD_OBJ_WEIGHT", "NONZEROS");
  }

  this->lib_->input_type_ = input_type;
  this->lib_->repartition(sublist, this->properties_, this->exportsSize_, this->imports_);
  this->computeNumberOfProperties();
  this->operation_already_computed_ = true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class Node>
void Partitioner<Node>::compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
template <class Node>
Teuchos::RCP< ::Tpetra::Map<int,int,Node> > Partitioner<Node>::createNewMap()
{
  ::Tpetra::Map<int,int,Node> *outputMap;

  createNewMap(outputMap);

  return( Teuchos::RCP< ::Tpetra::Map<int,int,Node> >(outputMap) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
template <class Node>
void Partitioner<Node>::createNewMap(::Tpetra::Map<int,int,Node> * &outputMap)
{
  if (!alreadyComputed()) {
    partition();
  }

  //Generate New Element List
  int myPID = this->input_map_->Comm().MyPID();
  int numMyElements = this->input_map_->NumMyElements();
  std::vector<int> elementList( numMyElements );
  if (numMyElements > 0)
    this->input_map_->MyGlobalElements( &elementList[0] );
  else
    this->input_map_->MyGlobalElements(NULL);

  int newGIDSize = numMyElements - this->exportsSize_;

  std::vector<int> myNewGID;

  if (newGIDSize > 0){
    myNewGID.resize(newGIDSize);
    std::vector<int>::iterator newElemsIter;
    std::vector<int>::const_iterator elemsIter;

    for (elemsIter = this->properties_.begin(), newElemsIter= myNewGID.begin() ;
         elemsIter != this->properties_.end() ; elemsIter ++) {
      if ((*elemsIter) == myPID) {
        (*newElemsIter) = elementList[elemsIter - this->properties_.begin()];
        newElemsIter ++;
      }
    }
  }
  //Add imports to end of list
  myNewGID.insert(myNewGID.end(), this->imports_.begin(), this->imports_.end());

  int *gidptr;
  if (myNewGID.size() > 0)
    gidptr = &myNewGID[0];
  else
    gidptr = NULL;

  outputMap = new ::Tpetra::Map<int,int,Node>(-1, myNewGID.size(), gidptr, 0, this->input_map_->Comm()); //MMW need to check

  return;
}
////////////////////////////////////////////////////////////////////////////////


} // namespace TPETRA

#endif //HAVE_TPETRA

}//namespace Isorropia

