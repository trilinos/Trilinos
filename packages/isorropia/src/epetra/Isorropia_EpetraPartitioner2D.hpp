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

#ifndef _Isorropia_EpetraPartitioner2D_hpp_
#define _Isorropia_EpetraPartitioner2D_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraOperator.hpp>
#include <Isorropia_Partitioner2D.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace Isorropia {

namespace Epetra {
  class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/

class Partitioner2D : virtual public Isorropia::Partitioner2D, virtual public Isorropia::Epetra::Operator  {
public:
  /** Constructor that accepts an Epetra_CrsGraph object, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RefCountPtr is used here because a
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
  Partitioner2D(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  /** Constructor that accepts an Epetra_CrsGraph object and a CostDescriber, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RefCountPtr is used here because a
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
  Partitioner2D(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
              Teuchos::RefCountPtr<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  /**
     Constructor that accepts an Epetra_RowMatrix object, called by
       API function create_partitioner().

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RefCountPtr is used here because a
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
  Partitioner2D(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  /**
     Constructor that accepts an Epetra_RowMatrix object and a
     CostDescriber, called by API function create_partitioner(). 

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RefCountPtr is used here because a
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
  Partitioner2D(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
              Teuchos::RefCountPtr<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  /** Destructor */
  virtual ~Partitioner2D();

  /** setParameters() is an internal Partitioner method which handles
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
  void partition(bool force_repartitioning=false);

  virtual void compute(bool forceRecomputing=false);

  //  /** An internal method which determines whether the
  //      method compute_partitioning() has already been
  //      called on this class instance.
  //  */
  //  bool partitioning_already_computed() const __deprecated;

  //  /** An internal method which returns the new partition ID for a given element that
  //     resided locally in the old partitioning.
  //  */
  //  int newPartitionNumber(int myElem) const __deprecated;

  /** An internal method which returns the number of elements in a given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
  int numElemsInPartition(int partition) const;

  /** An internal method which fills caller-allocated list (of length len) with the
      global element ids to be located in the given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
  void elemsInPartition(int partition, int* elementList, int len) const;

  Teuchos::RefCountPtr<Epetra_Map> createNewMap();

  int createNewMaps(Teuchos::RefCountPtr<Epetra_Map> domainMap, 
		    Teuchos::RefCountPtr<Epetra_Map> rangeMap);

  int partitionVectors();


};//class Partitioner2D

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif


