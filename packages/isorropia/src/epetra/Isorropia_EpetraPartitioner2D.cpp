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

// Assumes we have Epetra and Zoltan support
#include <Isorropia_EpetraPartitioner2D.hpp>
#include <Isorropia_Zoltan_Repartition.hpp>
#include <Isorropia_EpetraZoltanLib.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

namespace Epetra {


  /** Constructor that accepts an Epetra_CrsGraph object, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}

  /** Constructor that accepts an Epetra_CrsGraph object and a CostDescriber, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param costs CostDescriber object which allows for user-specified
       weights of varying types to be provided to the partitioner.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}


  /**
     Constructor that accepts an Epetra_RowMatrix object, called by
       API function create_partitioner().

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}


  /**
     Constructor that accepts an Epetra_RowMatrix object and a
     CostDescriber, called by API function create_partitioner(). 

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param costs CostDescriber object which allows for user-specified
       weights of varying types to be provided to the partitioner.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}


  /** Destructor */
Partitioner2D::~Partitioner2D(){}

////////////////////////////////////////////////////////////////////////////////
  /** setParameters() is an internal Partitioner2D method which handles
      the parameters from a Teuchos::ParameterList object. 

      The input
      ParameterList object is copied into an internal ParameterList
      attribute, and no reference to the input object is held after
      this function returns. (Thus, the input paramlist object may be
      altered or destroyed as soon as this method returns.)<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.
   */

  /**  partition is a method that computes 
       a rebalanced partitioning for the data in the object
      that this class was constructed with.

      \param force_repartitioning Optional argument defaults to false. By
         default, compute_partitioning() only does anything the first time
         it is called, and subsequent repeated calls are no-ops. If the user's
         intent is to re-compute the partitioning (e.g., if parameters
         or other inputs have been changed), then setting this flag to
         true will force a new partitioning to be computed.
   */
////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::
partition(bool force_repartitioning)
{
  std::string partitioning_method_str("PARTITIONING METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  if(partitioning_method == "UNSPECIFIED")
  {
    throw Isorropia::Exception("PARTITIONING_METHOD parameter must be specified.");
  }

  if(partitioning_method != "HGRAPH2D_FINEGRAIN")
  {
    throw Isorropia::Exception("PARTITIONING_METHOD parameter must be HGRAPH2D_FINEGRAIN.");
  }

  if (alreadyComputed() && !force_repartitioning)
    return;

  // Determine whether graph input or matrix input is used
  if (input_graph_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, costs_, 
             Library::hgraph2d_finegrain_input_));
  else
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_, 
             Library::hgraph2d_finegrain_input_));

  std::string zoltan("ZOLTAN");
  Teuchos::ParameterList &sublist = paramlist_.sublist(zoltan);

  if (paramlist_.isParameter("NUM PARTS")) 
  {
    sublist.set("NUM_GLOBAL_PARTS", paramlist_.get<std::string>("NUM PARTS"));
  }
  if (paramlist_.isParameter("IMBALANCE TOL")) 
  {
    sublist.set("IMBALANCE_TOL", paramlist_.get<std::string>("IMBALANCE TOL"));
  }

  lib_->repartition(sublist, properties_, exportsSize_, imports_);
  computeNumberOfProperties();
  operation_already_computed_ = true;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::
compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}
////////////////////////////////////////////////////////////////////////////////

//  /** An internal method which determines whether the 
//      method compute_partitioning() has already been
//      called on this class instance.
//  */
//bool Partitioner2D::partitioning_already_computed() const 
//{
//  return (alreadyComputed());
//}

//  /** An internal method which returns the new partition ID for a given element that
//     resided locally in the old partitioning.
//  */
//int Partitioner2D::newPartitionNumber(int myElem) const
//{
//  return ((*this)[myElem]);
//}

  /** An internal method which returns the number of elements in a given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
int Partitioner2D::numElemsInPart(int part) const
{
  return (numElemsWithProperty(part));
}

////////////////////////////////////////////////////////////////////////////////
  /** An internal method which fills caller-allocated list (of length len) with the
      global element ids to be located in the given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::elemsInPart(int partition, int* elementList, int len) const 
{
  //MMW
  std::cout << "MMW::NEED to reimplement" << std::endl;

  return (elemsWithProperty(partition, elementList, len));
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int
Partitioner2D::createNewMaps(Teuchos::RCP<Epetra_Map> domainMap, 
			     Teuchos::RCP<Epetra_Map> rangeMap) 
{
  //MMW
  std::cout << "MMW::NEED to reimplement" << std::endl;

  /*
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
       elemsIter != properties_.end() ; elemsIter ++) 
  {
    if ((*elemsIter) == myPID) 
    {
      (*newElemsIter) = elementList[elemsIter - properties_.begin()];
      newElemsIter ++;
    }
  }
  //Add imports to end of list
  myNewGID.insert(myNewGID.end(), imports_.begin(), imports_.end());

  Teuchos::RCP<Epetra_Map> target_map =
    Teuchos::rcp(new Epetra_Map(-1, myNewGID.size(), &myNewGID[0], 0, input_map_->Comm()));

  return(target_map);

  */
  return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int
Partitioner2D::partitionVectors() 
{
  //MMW
  std::cout << "MMW::NEED to implement" << std::endl;
  return 0;
}
////////////////////////////////////////////////////////////////////////////////

} // namespace EPETRA

}//namespace Isorropia

