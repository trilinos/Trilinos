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

#ifndef _Isorropia_Orderer_hpp_
#define _Isorropia_Orderer_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Operator.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new ordering and
  describing the layout of elements in the new order.

  If the methods which describe the new ordering (e.g., operator[],
  etc.) are called before order() has been called, behavior is not
  well defined. Implementations will either return empty/erroneous
  data, or throw an exception. In most cases, implementations will
  probably call order() internally in a constructor or factory method,
  so this won't usually be an issue.
*/
class Orderer : virtual public Operator {
public:

  /** Destructor */
  virtual ~Orderer() {}

  /** Method which computes a new ordering.

     \param forceOrdering Optional argument defaults to false.
        Depending on the implementation, order() should
        only perform a reordering the first time it is called, and
        subsequent repeated calls are no-ops. If the user's intent is
        to re-compute the ordering (e.g., if parameters or other
        inputs have been changed), then setting this flag to true
        will force a new ordering to be computed.
   */
  virtual void order(bool forceOrdering=false) = 0;

};//class Orderer

}//namespace Isorropia

#endif

