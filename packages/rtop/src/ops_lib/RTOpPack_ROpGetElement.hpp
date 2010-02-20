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

#ifndef RTOPPACK_ROP_GET_ELEMENT_HPP
#define RTOPPACK_ROP_GET_ELEMENT_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpGetElementEleWiseReductionOp {
public:
  /** \brief . */
  ROpGetElementEleWiseReductionOp(const Ordinal &global_i = -1)
    :global_i_(global_i)
    {}
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0, Scalar &reduct) const
    {
      if (i == global_i_) {
        reduct = v0;
      }
    }
private:
  Ordinal global_i_;
};


/** \brief Returns the value of the element at index <tt>global_i</tt>.
 *
 * Warning! If the element is not found, then 0 is returned!
 */
template<class Scalar>
class ROpGetElement
  : public ROp_1_CoordVariantScalarReduction<
        Scalar,
        Scalar,
        ROpGetElementEleWiseReductionOp<Scalar>
      >
{
public:
  /** \brief . */
  ROpGetElement(const Ordinal &global_i)
    {
      this->setOpNameBase("ROpGetElement");
      globalIndex(global_i);
      this->initReductObjValue(ScalarTraits<Scalar>::zero());
    }
  /** \brief . */
  void globalIndex(const Ordinal &global_i)
    { 
      this->setEleWiseReduction(ROpGetElementEleWiseReductionOp<Scalar>(global_i));
    }
  /** \brief . */
  Scalar operator()(const ReductTarget& reduct_obj ) const
    {
      return this->getRawVal(reduct_obj);
    }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_GET_ELEMENT_HPP
