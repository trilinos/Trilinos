// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_COUNT_NAN_INF_HPP
#define RTOPPACK_ROP_COUNT_NAN_INF_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Reduction operator that counts the number of entries that are NaN
 * or Inf.
 */
RTOP_ROP_1_REDUCT_SCALAR( ROpCountNanInf,
  index_type, // Reduction object type
  REDUCT_TYPE_SUM // Reduction object reduction
  )
{
  reduct += ( ScalarTraits<Scalar>::isnaninf(v0) ? 1 : 0 );
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_COUNT_NAN_INF_HPP
