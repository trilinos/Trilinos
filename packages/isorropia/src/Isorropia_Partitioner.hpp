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

#ifndef _Isorropia_Partitioner_hpp_
#define _Isorropia_Partitioner_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new partitioning and
  describing the layout of elements in the new partition (the parts).

  If the methods which describe the new partitioning (e.g.,
  operator [], elemsInPart()) are called before compute_partitioning()
  has been called, behavior is not well defined. Implementations will
  either return empty/erroneous data, or throw an exception. In most
  cases, implementations will probably call compute_partitioning()
  internally in a constructor or factory method, so this won't usually
  be an issue.
*/
class Partitioner : virtual public Operator {
public:

  /** Destructor */
  virtual ~Partitioner() {}

  /** Method which does the work of computing a new partitioning.
     Implementations of this interface will typically be constructed
     with an object or information describing the existing ('old')
     partitioning. This method computes a 'new' rebalanced
     partitioning for that input data.

     \param forceRepartitioning Optional argument defaults to false.
        Depending on the implementation, partitioning() should
        only perform a repartitioning the first time it is called, and
        subsequent repeated calls are no-ops. If the user's intent is
        to re-compute the partitioning (e.g., if parameters or other
        inputs have been changed), then setting this flag to true
        will force a new partitioning to be computed.

     \sa Isorropia::Operator::compute()
   */
  virtual void partition(bool forceRepartitioning=false) = 0;



  //@{ \name Usage not recommanded
  /** Return the new part ID for a given element that
      resided locally in the old partitioning.

      \deprecated A better way to have the same results is to use
      Isorropia::Operator::operator[]

      \param[in] myElem local ID of element we want to know in which part
      it belongs to.
      \return new part number

      \sa Isorropia::Operator::operator[]
   */
  virtual  int newPartNumber(int myElem) const = 0;
  //@}


  /** Return the number of LOCAL elements in a given part.

      \param[in] part the part ID we want to know the number of local
      elements.

      \return number of local elements that belongs to the
      given part.

      \sa Isorropia::Operator::numElemsWithProperty()
   */
  virtual  int numElemsInPart(int part) const = 0;


  /** Fill user-allocated list (of length len) with the
      local element ids to be located in the given part

      \param[in] part the part ID we consider

      \param[out] elementList array of elements that belongs to this
      part ID, must be allocated by user with size at least @c len

      \param[in] len maximum number of elements we can put in the
      array. Usually, may be the result of
      Isorropia::Partitioner::numElemsInPart().  .

      \sa Isorropia::Operator::elemsWithProperty()
   */
  virtual  void elemsInPart(int part,
			    int* elementList,
			    int len) const = 0;

  //@{ \name Deprecated methods

  /** The deprecated way to compute partitioning.

      \param[in] forceRepartitioning Optional parameter to recompute the
      partitioning even one has been already computed.

      \deprecated This method uses the old name and will be removed in
      next version, use Isorropia::Partitioner::partition() (or
      possibly Isorropia::Operator::compute() )

      \sa Isorropia::Partitioner::partition()
      \sa Isorropia::Operator::compute()
   */
  virtual __deprecated void compute_partitioning(bool forceRepartitioning=false) {
    return (partition(forceRepartitioning));
  }

  /** Query whether a partitioning as already been successfully computed.

      \deprecated This method uses the old name and will be removed in
      next version, use Isorropia::Operator::alreadyComputed()

      \sa Isorropia::Operator::alreadyComputed()
   */
  virtual __deprecated bool partitioning_already_computed() const = 0 ;


  /** Return the new partition ID for a given element that
      resided locally in the old partitioning.

      \deprecated This method uses the old name and will be removed in
      next version, use Isorropia::Partitioner::newPartNumber(), or
      even better Isorropia::Operator::operator[]()

      \sa Isorropia::Partitioner::newPartNumber()
      \sa Isorropia::Operator::operator[]()
   */
  virtual __deprecated int newPartitionNumber(int myElem) {
    return (newPartNumber(myElem));
  }


  /** Return the number of LOCAL elements in a given partition.

      \deprecated This method uses the old name and will be removed in
      next version, use Isorropia::Partitioner::numElemsInPart()

      \sa Isorropia::Partitioner::numElemsInPart()
   */
  virtual __deprecated int numElemsInPartition(int partition) {
    return (numElemsInPart(partition));
  }

  /** Fill user-allocated list (of length len) with the
      local element ids to be located in the given partition.

      \deprecated This method uses the old name and will be removed in
      next version, use Isorropia::Partitioner::elemsInPart()

     \sa Isorropia::Partitioner::elemsInPart()
   */
  virtual __deprecated void elemsInPartition(int part,
                                int* elementList,
                                int len) {
    return (elemsInPart(part, elementList, len));
  }

  //@}

};//class Partitioner

}//namespace Isorropia

#endif

