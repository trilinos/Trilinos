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

  If the accessors methods are called before the computation of the
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
      attribute. Instances of this interface should not retain a reference
      to the input ParameterList after this method returns.
  */
  virtual void setParameters(const Teuchos::ParameterList& paramlist) = 0;

  /** Method which does the work of computing a new partitioning.
     Implementations of this interface will typically be constructed
     with an object or information describing the existing ('old')
     partitioning. This method computes a 'new' rebalanced
     partitioning for that input data.

     \param forceRecomputing Optional argument defaults to false.
     Depending on the implementation, compute() should
     only perform a computation the first time it is called, and
     subsequent repeated calls are no-ops. If the user's intent is
     to re-compute the results (e.g., if parameters or other
     inputs have been changed), then setting this flag to true
     will force a new result to be computed.
   */
  virtual void compute(bool forceRecomputing=false) = 0;

  /** Query whether the computation has already been called.
   */
  virtual bool alreadyComputed() const = 0;


  /** Return the number of differtent "properties", like the number of
   *  colors used for the overall graph/matrix.
   */
  virtual int numProperties() const = 0;

  /** Return the "property" for a given element that resided locally
   *  in the old partitioning.
   */
  virtual const int& operator[](int myElem) const = 0;

  /** Return the number of elements with the given property
   */
  virtual int numElemsWithProperty(int property) const = 0;

  /** Fill user-allocated list (of length len) with the local element
   *  ids of the LOCAL elements with the given property.
   */
  virtual void elemsWithProperty(int property,
				 int* elementList,
				 int len) const = 0;

};//class Operator

}//namespace Isorropia

#endif

