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


  /** Give access of the part assignments array that is owned by the current
      processor.

      \param[out] size Number of elements in the array.

      \param[out] array Pointer to the the part assignements array inside
                        the object.

      \remark This pointer is only significant if the object still exists.
      Otherwise, you must use \see Isorropia::Operator::extractPartsCopy()

      \sa Isorropia::Operator::extractPropertiesView()
   */
  virtual int extractPartsView(int& size,
			       const int*& array) const {
    return extractPropertiesView(size, array);
  }


  /** Copy a part of the part assignment array.

      \param[in] len of the array given by the user.

      \param[out] size Number of elements in the array.

      \param[out] array Array of part assignments. Allocated by the user with
                        a size of at least @c len elements.

      \remark Memory space which is not useful in the array is not
      initialized or used in this method.

      \sa Isorropia::Operator::extractPropertiesCopy()
   */
  virtual int extractPartsCopy(int len,
			       int& size,
			       int* array) const {
    return extractPropertiesCopy(len, size, array);
  }

};//class Partitioner

}//namespace Isorropia

#endif

