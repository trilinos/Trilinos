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

#ifndef _Isorropia_Operator_hpp_
#define _Isorropia_Operator_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new partitioning/coloring/
    ordering and exploiting their results.

  If the accessor methods are called before the computation of the
  result (by a method like compute()) has been called, behavior is not
  well defined.  Implementations will either return empty/erroneous
  data, or throw an exception. In most cases, implementations will
  probably call compute_partitioning() internally in a constructor or
  factory method, so this won't usually be an issue.
*/
class Operator {
public:

  /** Destructor */
  virtual ~Operator() {}

  /** Set parameters for the Operator instance. The contents of the
      input paramlist object are copied into an internal ParameterList
      attribute. Instances of this interface should not retain a
      reference to the input ParameterList after this method returns.

      \param[in] paramlist List of parameters that the user wants to use.
  */
  virtual void setParameters(const Teuchos::ParameterList& paramlist) = 0;

  /** Method which does the work of computing a new
     partitioning/coloring/ordering, depending on the child class
     used.

     \param forceRecomputing Optional argument defaults to false.
     Depending on the implementation, compute() should
     only perform a computation the first time it is called, and
     subsequent repeated calls are no-ops. If the user's intent is
     to re-compute the results (e.g., if parameters or other
     inputs have been changed), then setting this flag to true
     will force a new result to be computed.
   */
  virtual void compute(bool forceRecomputing=false) = 0;


  //@{ \name Accessors
  /** Query whether the computation has already been called.

      \return @c True if the computation has already been done, @c
      False otherwise.
   */
  virtual bool alreadyComputed() const = 0;


  /** Return the number of different values used for "properties".

      For example, the number of colors or the number of parts used
      for the overall graph/matrix.

      \return Global number of values for properties

      \remark Infact, it returns the upper bound of the interval of
      taken values. For example, for the colors "1,2,4"" , it will
      return "4"
   */
  virtual int numProperties() const = 0;

  /** Return the number of different values used for "properties"
      for this process only.

      \return Local number of values for properties

   */
  virtual int numLocalProperties() const = 0;


  /** Return the "property" for a given element that resided locally.

      \param[in] myElem the local ID of the element we want to know
      the property.

      \return property associated to the local element.
   */
  virtual const int& operator[](int myElem) const = 0;


  /** Return the number of \b LOCAL elements with the given property.

      \param[in] property Value of the property to consider.

      \return Number of \b local elems which have this property.
   */
  virtual int numElemsWithProperty(int property) const = 0;


  /** Fill user-allocated list (of length len) with the local element
      ids of the LOCAL elements with the given property.

      \param[in] property  Value of the property to consider.

      \param[out] elementList User allocated array (of size at least
      \c len) of local ID that have the asked property.

      \param[in] len Maximum lenght for the array. If \c len is
      greater than the result of numElemsWithProperty() for \c
      property, only the first and relevant elements are filled.

      \remark Memory space which is not useful in the array is not
      initialized or used in this method.
   */
  virtual void elemsWithProperty(int property,
				 int* elementList,
				 int len) const = 0;


  /** Give access of the property array that is owned by the current
      processor.

      \param[out] size Number of elements in the array.

      \param[out] array Pointer to the the properties array inside
                        the object.

      \remark This pointer is only significant if the object still exists.
      Otherwise, you must use \see Isorropia::Operator::extractPropertiesCopy().

      \sa Isorropia::Operator::extractPropertiesCopy()
   */
  virtual int extractPropertiesView(int& size,
				    const int*& array) const = 0;


  /** Copy a part of the property array.

      \param[in] len of the array given by the user.

      \param[out] size Number of elements in the array.

      \param[out] array Array of properties. Allocated by the user with
                        a size of at least @c len elements.

      \remark Memory space which is not useful in the array is not
      initialized or used in this method.

      \sa Isorropia::Operator::extractPropertiesView()
   */
  virtual int extractPropertiesCopy(int len,
				    int& size,
				    int* array) const = 0;

  //@}

};//class Operator

}//namespace Isorropia

#endif

