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

#ifndef RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP
#define RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Element-wise transformation operator for TOpEleWiseProdUpdate. */
template<class Scalar>
class TOpEleWiseProdUpdateEleWiseTransformation
{
public:
  TOpEleWiseProdUpdateEleWiseTransformation( const Scalar &alpha )
    : alpha_(alpha)
    {}
  void operator()( const Scalar &v0, Scalar &z0 ) const
    {
      z0 *= alpha_ * v0;
    }
private:
  Scalar alpha_;
};


/** \brief Element-wise product update transformation operator: 
 * <tt>z0[i] *= alpha*v0[i], i=0...n-1</tt>.
 */
template<class Scalar>
class TOpEleWiseProdUpdate
  : public TOp_1_1_Base<Scalar, TOpEleWiseProdUpdateEleWiseTransformation<Scalar> >
{
public:
  typedef TOp_1_1_Base<Scalar, TOpEleWiseProdUpdateEleWiseTransformation<Scalar> > base_t;
  /** \brief . */
  TOpEleWiseProdUpdate( const Scalar &alpha )
    : base_t(TOpEleWiseProdUpdateEleWiseTransformation<Scalar>(alpha))
    {
      this->setOpNameBase("TOpEleWiseProdUpdate");
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ELE_WISE_PROD_UPDATE_HPP
