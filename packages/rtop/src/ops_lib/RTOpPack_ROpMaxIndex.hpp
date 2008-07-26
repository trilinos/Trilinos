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

#ifndef RTOPPACK_ROP_MAX_INDEX_HPP
#define RTOPPACK_ROP_MAX_INDEX_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpMaxIndexEleWiseReductionOp {
public:
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0,
    ScalarIndex<Scalar> &reduct) const
    {
      if(
        v0 > reduct.scalar
        ||
        (v0 == reduct.scalar && i < reduct.index)
        )
      {
        reduct = ScalarIndex<Scalar>(v0, i);
      }
    }
};


/** \brief. */
template<class Scalar>
class ROpMaxIndexReductObjReductionOp {
public:
  /** \brief . */
  void operator()(
    const ScalarIndex<Scalar>& in_reduct, ScalarIndex<Scalar>& inout_reduct
    ) const
    {
      if(
        in_reduct.scalar > inout_reduct.scalar
        ||
        (
          in_reduct.scalar == inout_reduct.scalar
          &&
          in_reduct.index < inout_reduct.index
          )
        )
      {
        inout_reduct = in_reduct;
      }
    }
};


/** \brief Returns the maximum element and its index:
 * <tt>result.scalar = x(k)</tt> and <tt>result.index = k</tt> such
 * that <tt>x(k) >= x(i)</tt> for <tt>i=0...n-1</tt> and <tt>k</tt> is
 * the minimum index to break ties.
 */
template<class Scalar>
class ROpMaxIndex
  : public ROp_1_CoordVariantScalarReduction<
      Scalar,
      ScalarIndex<Scalar>,
      ROpMaxIndexEleWiseReductionOp<Scalar>,
      ROpMaxIndexReductObjReductionOp<Scalar> >
{
public:
  /** \brief . */
  ROpMaxIndex()
    {
      this->setOpNameBase("ROpMaxIndex");
      this->initReductObjValue(
        ScalarIndex<Scalar>(-ScalarTraits<Scalar>::rmax(), -1));
    }
  /** \brief . */
  ScalarIndex<Scalar> operator()(const ReductTarget& reduct_obj) const
    { return this->getRawVal(reduct_obj); }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MAX_INDEX_HPP
