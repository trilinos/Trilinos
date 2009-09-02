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

#ifndef RTOPPACK_ROP_WEIGHTED_NORM2_HPP
#define RTOPPACK_ROP_WEIGHTED_NORM2_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpWeightedNorm2EleWiseReduction
{
public:
  void operator()( const Scalar &v0, const Scalar &v1, Scalar &reduct ) const
    {
      reduct += v0 * ScalarTraits<Scalar>::conjugate(v1)*v1;
    }
};


/** \brief Weighted Two (Euclidean) norm reduction operator: <tt>result =
 * sqrt( sum( v0[i]*conj(v1[i])*v1[i], i=0...n-1 ) )</tt>.
 */
template<class Scalar>
class ROpWeightedNorm2
  : public ROp_2_ScalarReduction<Scalar, Scalar,
      ROpWeightedNorm2EleWiseReduction<Scalar> >
{
public:
  /** \brief . */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief . */
  ROpWeightedNorm2()
    {
      this->setOpNameBase("ROpWeightedNorm2");
    }
  /** \brief . */
  typename ST::magnitudeType operator()(const ReductTarget& reduct_obj ) const
    { return ST::magnitude(ST::squareroot(this->getRawVal(reduct_obj))); }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_WEIGHTED_NORM2_HPP
