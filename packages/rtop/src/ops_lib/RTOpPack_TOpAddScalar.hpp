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

#ifndef RTOPPACK_TOP_ADD_SCALAR_HPP
#define RTOPPACK_TOP_ADD_SCALAR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpAddScalar. */
template<class Scalar>
class TOpAddScalarEleWiseTransformation
{
public:
  TOpAddScalarEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( Scalar &z0 ) const
    {
      z0 += alpha_;
    }
private:
  Scalar alpha_;
};


/** \brief Add a scalar to a vector transformation operator: <tt>z0[i] +=
 * alpha, i=0...n-1</tt>.
 */
template<class Scalar>
class TOpAddScalar
  : public TOp_0_1_Base<Scalar, TOpAddScalarEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_0_1_Base<Scalar, TOpAddScalarEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpAddScalar( const Scalar &alpha )
    : base_t(TOpAddScalarEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpAddScalar");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ADD_SCALAR_HPP
